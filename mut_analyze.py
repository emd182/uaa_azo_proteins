# This script is for analyzing mutant proteins with unnatural amino acids included
# Editor: Elliott Dolan, Khare Lab, Rutgers University, 2020

#!/usr/bin/python

#Library importing
import pyrosetta as pr
import argparse, pickle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import sys, os
from glob import glob
from pyrosetta.rosetta.protocols.constraint_generator import \
    AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.core.select.residue_selector import \
    NeighborhoodResidueSelector, ResidueIndexSelector, OrResidueSelector,\
    ResidueNameSelector, NotResidueSelector, InterGroupInterfaceByVectorSelector,\
    AndResidueSelector
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.protocols.enzdes import AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.minimization_packing import \
    MinMover, PackRotamersMover
from pyrosetta.rosetta.core.pack.task.operation import \
    IncludeCurrent, ExtraRotamers, OperateOnResidueSubset, \
    RestrictToRepackingRLT, PreventRepackingRLT
from pyrosetta.rosetta.protocols.grafting.simple_movers import \
    DeleteRegionMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.simple_metrics.metrics import \
    TotalEnergyMetric, RMSDMetric, InteractionEnergyMetric, \
    SasaMetric
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import \
    PerResidueEnergyMetric

#Script starting from here
def parse_args():
    info = """
        This script should take a protein file and a few other parameters
        and output a series of files with a point mutation and relaxed proteins.
        """
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-p', '--wt_pdb', type=str, required=True,
                        help='Comma separated set of the relaxed PDB files with\
                        and without the ligand. Running the mutation_design.py \
                        script with --residue_num 0 with and without the ligand \
                        will generate these files.')
    parser.add_argument('-dir', '--source_directory', type=str, required=True,
                        help='Directory where the files that will be compared to\
                        will be found. Should contain mutants from the program \
                        mut_select.py - containing proteins in both isomers and \
                        with and without the bound ligand.')
    parser.add_argument('-l', '--ligand_type', type=str, required=True,
                        choices=['protein','ligand'],
                        help='Ligand type for identification - either 3-letter\
                        ligand name.')
    parser.add_argument('-res', '--residue_set', type=str, 
                        required=any(x in ['protein'] for x in sys.argv),
                        help='Comma delimited residue list for a peptide sequence \
                        behaving as a ligand in the mutation detection script.')
    parser.add_argument('-lig', '--ligand_name', type=str,  
                        required=any(x in ['ligand'] for x in sys.argv),
                        help='Ligand params file for ligand identification. \
                        Ligand params file should have a 3 letter name with \
                        the same 3 letters for the lig name.')
    parser.add_argument('-uaa', '--unnatural', action='store_true', required=False,
                        help='Stating if unnaturals are used - will require params files')
    parser.add_argument('-enz', '--enzdes_cstfiles', type=str, required=False,
                        help='Add for implementation of constraints used in the \
                        design of the protein. Run without for unconstrained scores - \
                        comma delimited for yes/no lig cst files')
    parser.add_argument('-sym', '--symmdef_file', type=str, required=False,
                        help='Add to consider symmetry for these proteins - \
                        needed by first removing first chains and then re-symmetrizing\
                        the protein.')
    parser.add_argument('-cis', '--cis_params_file', type=str, \
                        required=any(x in ['--uaa', '--unnatural'] for x in sys.argv), \
                        help='Cis state params files for unnatural amino acid inclusions')
    parser.add_argument('-trans', '--trans_params_file', type=str, \
                        required=any(x in ['--uaa', '--unnatural'] for x in sys.argv),\
                        help='Trans state params files for unnatural amino acid inclusions')
    parser.add_argument('-pre', '--previously_run_files', type=str, required=False, \
                        default=None, help="Load data saved from a previous run.")
    parser.add_argument('-per', '--per_resi', action='store_true', default=False, \
                        help="Make graphs for each residue for all proteins.")
    parser.add_argument('-lim', '--value_limits', type=float, required=False, \
                        default=15.0, help="Add limits of total scores in the case the\
                        apo states are diverge too greatly from the wildtype.")
    parser.add_argument('-out', '--out_file', type=str, required=False, \
                        help='Outputs a pickled set of information for later analysis.\
                        Or to prevent having to repeat the analysis.')
    parser.add_argument('-od', '--out_directory', type=str, required=False, \
                        help='Output Directory fo all output files.')

    return parser.parse_args()

def print_out(outstring):
    print("emd182::"+str(outstring))

def if_list_key_not_in_dict(newkey, dictionary):
    if newkey not in dictionary.keys():
        dictionary[newkey] = []

def out_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        print_out("Creating directory : " + directory + " as it doesn't exist")
    else:
        print_out("Directory " + directory + " already exists")

def total_energy(pose, score_function, selection=None):
    """
    Calculates total energy of a pose using a TotalEnergyMetric.
    If a residue selector is provided, it will calculate the total energy
    of that selection rather than the total pose energy.
    """
    tem = TotalEnergyMetric()
    tem.set_scorefunction(score_function)

    if selection:
        tem.set_residue_selector(selection)
    return tem.calculate(pose)

def total_sasa(pose):
    sm = SasaMetric()
    return sm.calculate(pose)

def residue_selection(selection_one, selection_two):
    #For now, just an or selector to group up the residues.
    first = ResidueIndexSelector(selection_one)
    second = ResidueIndexSelector(selection_one)
    both = OrResidueSelector(first, second)
    return both

#Simple definition file to convert a residue selector applied to a pose
#Into a vector of residue ints.
def selector_to_vector(residue_selector, pose):
    rosetta_vector = pr.rosetta.utility.vector1_unsigned_long()
    for res_int in [int(x) for x in \
                get_residues_from_subset(residue_selector.apply(pose)) ]:
        rosetta_vector.append(res_int)
    return rosetta_vector

def file_interpreter_to_dict(pdb_file_name):
    print_out("Processing this pdbID:" + pdb_file_name)
    split = pdb_file_name.split('/')[-1].strip('.pdb').split('_')
    pdb_info = {}
    pdb_info['pdb'] = split[0]
    pdb_info['mut_res'] = split[1]
    pdb_info['uaa'] = split[2]
    pdb_info['cisortrans'] = split[2][-1]
    pdb_info['contain_ligand'] = bool('yes' in split[3])
    pdb_info['n'] = int(split[4])
    return pdb_info

def apply_match_constraints(pose, match_cst_file):
    AoRMcsts = AddOrRemoveMatchCsts()
    AoRMcsts.set_cst_action(pr.rosetta.protocols.enzdes.CstAction.ADD_NEW)
    AoRMcsts.cstfile(match_cst_file)
    AoRMcsts.apply(pose)

def dictionaried_poses(pdb_file_names, score_function, cst_files):
    pose_set = {}
    for pdb_file_name in pdb_file_names:
        pose_set[pdb_file_name] = pr.pose_from_pdb(pdb_file_name)
        score_function(pose_set[pdb_file_name])
        for cst in cst_files.keys():
            print(cst)
            if cst in pdb_file_name:
                apply_match_constraints( pose_set[pdb_file_name], cst_files[cst])
    return pose_set

def scoring_poses(pose_dict, info_dict, score_function, ifsct=None):
    scores = {}
    for pose in pose_dict.keys():
        current_scores = {}
        print_out("Scoring pose: " + pose)
        current_scores['total_' + info_dict[pose]['cisortrans'] + '_' + \
            str(info_dict[pose]['contain_ligand'])] \
            = total_energy(pose_dict[pose], score_function)
        current_scores['sasa_' + info_dict[pose]['cisortrans'] + '_' + \
            str(info_dict[pose]['contain_ligand'])] \
            = total_sasa(pose_dict[pose])
        current_scores['per_resi'] = {}
        #Per residue scoring
        prem = PerResidueEnergyMetric()
        prem.set_scorefunction(score_function)
        e_vals = prem.calculate(pose_dict[pose])
        current_scores['per_resi'] = [ e_vals[i] for i in \
                range(1, pose_dict[pose].total_residue() + 1) ]
        scores[pose] = current_scores
    return scores

def interaction_energy_addition(pose, lig_res_selector, sf, add_to_dict, cisortrans):
    not_lig = NotResidueSelector()
    not_lig.set_residue_selector(lig_res_selector)
    interact_metric = InteractionEnergyMetric()
    interact_metric.set_scorefunction(sf)
    interact_metric.set_residue_selectors(lig_res_selector, not_lig)
    add_to_dict['interE_'+cisortrans] = interact_metric.calculate(pose)

def manage_cst_files(cst_dict, cst_string):
    if ',' in cst_string:
        for cst in cst_string.split(','):
            if 'noligcst' in cst:
                cst_dict['nolig'] = cst
            else:
                cst_dict['yeslig'] = cst
    else:
        cst_dict['yeslig'] = cst_string

def make_histogram(data, bins, xlabel, ylabel, title):
    print_out("not written yet")
    y_values = eval(data)
    x = list(range(len(y_values)))
    plt.bar(x, y_values, width=1, color='blue')
    plt.xticks(x, x_ticks)

def add_to_legend(legend_list, colo, label):
    legend_list.append(mpatches.Patch(color=colo, label=label))

def make_scatterplot(xvalues, yvalues, xlabel, ylabel, title, filename=None,\
        colors='b', xylimits=[[-10,10],[-10,10]], markers='o', outdir='./',\
        tick_labels=None, lines=None):
    xvals = xvalues
    yvals = yvalues
    plt.scatter(xvals, yvals, s=5, c=colors, marker=markers)
    if lines:
        for line in lines:
            plt.plot(line[0],line[1], '-k', linewidth=1)
    axes = plt.gca()
    axes.set_xlim(xylimits[0])
    axes.set_ylim(xylimits[1])
    if tick_labels:
        for num in range(0, len(xvals)):
            plt.annotate(tick_labels[num], (xvals[num],yvals[num]), \
                textcoords="offset points",xytext=(0,5), ha='center',\
                fontsize='x-small')
    #plt.plot([-10,10],[-10,10], '-ok')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    fig = plt.gcf()
    fig.set_size_inches(5,5)
    #plt.tight_layout()
    plt.show()
    if not filename:
       filename = title.replace(' ', '_')
    out_file = '/'.join([outdir,filename + '.png'])
    fig.savefig(out_file, dpi=150)#, bbox_inches='tight')
    plt.close('all')
    print_out("writing png file: " + out_file)

def make_bar_plot(bar_values, bar_labels, xlabel, ylabel, title, \
            ymax=10.0, ymin=-10.0, control=0.0, upper_mask=None,  \
            lower_mask=None, out_dir='./', legend_info={'blue':'normalized score'}):
    normal = np.asarray([ x-control for x in bar_values ])
    x = np.asarray(list(range(len(normal))))
    legends = []
    if upper_mask is not None:
        plt.bar(x[upper_mask], normal[upper_mask], width=0.8, align='center', color='red')
        add_to_legend(legends, 'red', legend_info['red'])
    else:
        upper_mask = normal < -1000000.0
    if lower_mask is not None:
        plt.bar(x[lower_mask], normal[lower_mask], width=0.8, align='center', color='green')
        add_to_legend(legends, 'green', legend_info['green'])
    else:
        lower_mask = normal < -1000000.0
    
    nor_mask = np.logical_not(np.logical_or(upper_mask, lower_mask))
    plt.bar(x[nor_mask], normal[nor_mask], width=0.8, align='center', color='blue')
    add_to_legend(legends, 'blue', legend_info['blue'])
    plt.legend(handles=legends)
    plt.xticks(ticks=x, labels=bar_labels, rotation='vertical', fontsize='small')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.ylim([ymin,ymax])
    plt.tight_layout()
    fig = plt.gcf()
    plt.show()
    out_directory(out_dir)
    out_file ='/'.join([ out_dir, title.replace(' ', '_') +'.png'])
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close('all')
    print_out("Writing png file: " + out_file)

def plot_res_ddg(ref_data, mut_data, name=None):
    
    # Plot
    fig, axes = plt.subplots(2, 1, sharex=True)
    axes[0].plot(ref_data, linewidth=0.3, color='red')
    axes[0].plot(mut_data, linewidth=0.3, color='black')
    axes[0].set_ylabel('Residue Energy', fontsize=9)
    axes[1].plot(mut_data - ref_data, linewidth=0.3, color='black')
    axes[1].set_ylabel('Mutant - Reference', fontsize=9)

    # Format Figure
    fig.set_size_inches(7.5, 5.0)
    fig.set_dpi(300)
    bottom, top, left, right = [0.05, 0.95, 0.1, 1]
    vcenter = bottom + ( top - bottom ) /2
    hcenter = left + ( right - left ) /2

    fig.tight_layout()
    fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right)

#Label axes
    fig.text(0.001, vcenter, 'Residue energy', fontsize=10, rotation='vertical', \
                ha='left', va='center')
    fig.text(hcenter, 0.001, 'Residue Number', fontsize=10, ha='center', va='bottom')
    fig.text(hcenter, 0.001, 'Per-Residue Energies', fontsize=11, ha='center', va='top')

    if name:
        plt.savefig(name)
        print_out("Making graph for: " + name)
    else: plt.show()
    plt.close()

#Main function starts here.
###FUNCTION IDEAS###
##add a count to determine the number of residues that change their rotameric
#state when they relax.
##measure changes in energy between original binding residues to the ligand
# compared to the remainding residues that still form interactions with the ligand.

def main(args):
    
    params = []
    if args.ligand_type == 'ligand':
        params.append(args.ligand_name)
        lig_name = args.ligand_name.split('/')[-1].strip('.params')
    if args.unnatural:
        params.append(args.cis_params_file)
        params.append(args.trans_params_file)
    
    #setting up constraint file conditionals, as well as other init arguments
    cstfiles = {}
    if args.enzdes_cstfiles:
        manage_cst_files(cstfiles, args.enzdes_cstfiles)
    
    init_args = ' '.join(['-run:preserve_header', '-extra_res_fa'] + params)

    #init_args 'yeslig' and 'nolig' either contain or don't contain 
    # the constraint files for the pdbs.
    
    if not args.previously_run_files:
    #Starting up rosetta with the appropriate params files -
    #Need to incorporate constraints somehow
        print_out("emd182::Starting rosetta with the following parameters: " \
                    + init_args)
        pr.init(init_args)

        #Setting up score functions
        sf = pr.rosetta.core.scoring.ScoreFunction()
        sfcst = pr.rosetta.core.scoring.ScoreFunction()
        sf.add_weights_from_file('ref2015')
        sfcst.add_weights_from_file('ref2015_cst')
        
        #Looking in the directory to determine what to group:
        pdb_list = glob('/'.join([args.source_directory, '*.pdb']))
    
        #Briefly setting up the wt poses
        wt_pose = {}
        for pdb_file_name in args.wt_pdb.split(','):
            print(pdb_file_name)
            pdb_list.remove(pdb_file_name)
            if 'yeslig' in pdb_file_name:
                #pr.init(init_args['yeslig'])
                wt_pose['yeslig'] = pr.pose_from_pdb(pdb_file_name)
            elif 'nolig' in pdb_file_name:
                #pr.init(init_args['yeslig'])
                wt_pose['nolig'] = pr.pose_from_pdb(pdb_file_name)
            else:
                print('Error, submitted wild type pdbs do not have yes or no ligand' + \
                        'name. Please submit properly named pdbs.')
                sys.exit()

        for wt_key in wt_pose.keys():
            sf(wt_pose[wt_key])

        print_out("Loaded and scored the wt_pose")
    
        #Residue selection for later analysis
        if args.ligand_type == 'protein':
            ligand = ResidueIndexSelector(args.residue_set)
        elif args.ligand_type == 'ligand':
            ligand = ResidueNameSelector()
            ligand.set_residue_name3(lig_name)
    
        pdb_dict = {}
        #Dictionary of pdb info - name, resi, uaa, cis/trans, ligand, nstruct
        group_set = {} 
        #Dictionary of [residue number] : pdb file name for analysis
        
        print_out("Interpreting the list of pdbs in the directory: " + \
                    args.source_directory )
        #Sorting file names into groups and into dictionaries
        for pdb_file_name in pdb_list:
            pdb_dict[pdb_file_name] = file_interpreter_to_dict(pdb_file_name)
            resi = pdb_file_name.split('/')[-1].strip('.pdb').split('_')[1]
            if_list_key_not_in_dict(resi, group_set)
            group_set[resi].append(pdb_file_name)
        
        scored_set = {} #Dictionary of scores at each residue index
        print_out('Taking groups of proteins and scoring them one by one.')
        #Taking protein groups and scoring them one by one.
        for residue_key in group_set.keys():
            print_out("Now on residue number: "+str(residue_key))
            cst_key = residue_key + "_cst"
            #Place poses in a dictionary with the key designated as the pdb_file_name.
            #Each item in a dictionary is a scored pose object.
            #Creates a new dictionary each time to allow for score analysis each round.
            resi_poses = dictionaried_poses(group_set[residue_key], sf, cstfiles)
            scored_set[residue_key] = scoring_poses(resi_poses, pdb_dict, sf)
            if cstfiles:
                scored_set[cst_key] = scoring_poses(resi_poses, pdb_dict, sfcst)
            for filename in resi_poses.keys():
                if pdb_dict[filename]['contain_ligand']:
                    interaction_energy_addition(resi_poses[filename], ligand, \
                        sf, scored_set[residue_key][filename], \
                        pdb_dict[filename]['cisortrans'])
                    if cstfiles:
                        interaction_energy_addition(resi_poses[filename], ligand, \
                            sfcst, scored_set[cst_key][filename], \
                            pdb_dict[filename]['cisortrans'])

        #Generating wt_data for the two wt poses.
        wt_data = {}
        wt_data['total_True'] = total_energy(wt_pose['yeslig'], sf)
        wt_data['sasa_True'] = total_sasa(wt_pose['yeslig'])
        interaction_energy_addition(wt_pose['yeslig'], ligand, sf, wt_data, 'True')

        prem = PerResidueEnergyMetric()
        prem.set_scorefunction(sf)
        e_vals = prem.calculate(wt_pose['yeslig'])
        wt_data['per_res_True'] = [ e_vals[i] for i in \
                range(1, wt_pose['yeslig'].total_residue() + 1) ]

        wt_data['total_False'] = total_energy(wt_pose['nolig'], sf)
        wt_data['sasa_False'] = total_sasa(wt_pose['nolig'])
        e_vals = prem.calculate(wt_pose['nolig'])
        wt_data['per_res_False'] = [ e_vals[i] for i in \
                range(1, wt_pose['yeslig'].total_residue() + 1) ]
        
        if cstfiles:
            wt_data['cst_total_True'] = total_energy(wt_pose['yeslig'], sfcst)
            interaction_energy_addition(wt_pose['yeslig'], ligand, \
                        sfcst, wt_data, 'True')
            wt_data['cst_total_False'] = total_energy(wt_pose['nolig'], sf)

            prem.set_scorefunction(sfcst)
            e_vals = prem.calculate(wt_pose['yeslig'])
            wt_data['cst_per_res_True'] = [ e_vals[i] for i in \
                    range(1, wt_pose['yeslig'].total_residue() + 1) ]
            e_vals = prem.calculate(wt_pose['nolig'])
            wt_data['cst_per_res_False'] = [ e_vals[i] for i in \
                    range(1, wt_pose['nolig'].total_residue() + 1) ]
        
        if args.out_file:
            print_out("Storing all the data into a pickle file")
            with open(args.out_file, 'wb') as savedfile:
                pickle.dump([scored_set, pdb_dict, group_set, wt_data], savedfile)
            print_out("Pickle has stored files into " + str(args.out_file))

    if args.previously_run_files:
        print_out("Loading data pickled earlier from this code: " + args.previously_run_files)
        with open(args.previously_run_files, 'rb') as loading_file:
            scored_set, pdb_dict, group_set, wt_data = pickle.load(loading_file)
        print_out("Loaded data into scored_set, pdb_dict, and group_set")
   
    """
    At this point, object -- scored_set -- is a dictionary of dictionaries that 
    contains all of the scored values - total score, sasa, and interaction score
    Goes by -- scored_set[residue_key][filename]
                -- this dictionary contains multiple types of scores,
                    - sasa, total, and inter_E
                -- for each filename.
    The object -- pdb_dict[any pdb filename] contains parsed information about the pdbs
            -'pdb', 'mut_res', 'uaa', 'cis/trans', 'contain_ligand', and 'n'
    The object -- group_set -- contains one set of information:
            -'resi' for the residue index, which contains a list of pdb file names.
    """
   
    ###Data analysis time - 
    collated = {}
    value_keys = []
    #First will be making graphs of each pdb structure on a per residue basis
    #print(scored_set.keys())
    #print(wt_data.keys())
    out_directory(args.out_directory)
    if args.per_resi:
        for resikey in scored_set.keys():
            out_dir = '/'.join([args.out_directory, resikey])
            out_directory(out_dir)
            for prokey,scores in scored_set[resikey].items():
                pdb_to_png = prokey.split('/')[-1].replace('.pdb','.png')  
                png_name = '/'.join([out_dir, pdb_to_png])  
                if 'yeslig' in prokey:
                    ref = 'per_res_True'
                else:
                    ref = 'per_res_False'
                if 'cst' in resikey:
                    ref = 'cst_' + ref
                
                plot_res_ddg(np.array(wt_data[ref]), \
                            np.array(scored_set[resikey][prokey]['per_resi']), \
                            name=png_name) 

    #Removes the per_resi key to prevent confusion from going on in the next section
    for resikey in scored_set.keys():
        for prokey in scored_set[resikey].keys():
            scored_set[resikey][prokey].pop('per_resi', None)

    #First function parses the scored_sets of data and sorts into this format:
    #collated[residue_number][scoretype] = [list of scores]
    for resikey in scored_set:
        keysets = {}
        for filename in scored_set[resikey]:
            for k,v in scored_set[resikey][filename].items():
                if_list_key_not_in_dict(k, keysets)
                keysets[k].append(v)
                if k not in value_keys:
                    value_keys.append(k)
        collated[resikey] = keysets
    skip_keys = ['per_resi', 'cst_per_resi']
    averages = {}
    for residue,dictionary in collated.items():
        averages[residue] = {}
        for k, v in dictionary.items():
            if k in skip_keys:
                continue
            averages[residue][k] = sum(v)/len(v)
    
    print(averages.keys())
    sorted_residues = sorted( set( [int(x.split('_')[0][0:-1]) \
                                    for x in averages.keys() ] ) ) 
    sorted_cstkeys = [x for x in averages.keys() if 'cst' in x]
    sorted_cstkeys = sorted(sorted_cstkeys, key=lambda x: float(x.split('_')[0][:-1]))
    sorted_reskeys = [x for x in averages.keys() if 'cst' not in x]
    sorted_reskeys = sorted(sorted_reskeys, key=lambda x: float(x[:-1]))
    print(sorted_reskeys)
    print(len(sorted_reskeys))
    print(sorted_cstkeys)
    print(len(sorted_cstkeys))
    print(sorted_residues)
    print(len(sorted_residues))

    #should be a list of numbers
#Keytypes should be - total_C_True total_T_True sasa_C_True sasa_T_True
#Keytypes should be - total_C_False total_T_False sasa_C_False sasa_T_False
#Keytypes should be - interE_C interE_T
    out_file_name_pdb = args.wt_pdb.split("/")[-1].split('_')[0]
#First thing to do is organize data by residue

    #First set of graphs are normalized against the wt.
    #sorted_resi
    for key_type in keysets:
        keybreak = key_type.split('_')
        wtkey = ''.join([keybreak[0],"_",keybreak[-1]])
        if keybreak[0] == 'interE':
            wtkey = ''.join([keybreak[0],"_True"])
        xlabel='REU'
        if 'sasa' in key_type:
            continue
        title = ' '.join(["Normalized", str(key_type), 'across', out_file_name_pdb])
        make_bar_plot([averages[x][key_type] for x in sorted_reskeys], \
                sorted_residues, "Residue Index", xlabel, title, \
                control=wt_data[wtkey], out_dir=args.out_directory)
        if cstfiles:
            title = ' '.join(["Normalized", str(key_type), 'across', \
                                out_file_name_pdb, '- with csts'])
            make_bar_plot([averages[x][key_type] for x in sorted_cstkeys], \
                    sorted_residues, "Residue Index", xlabel, title, \
                    control=wt_data[wtkey], out_dir=args.out_directory)
            title = ' '.join(["Normalized", str(key_type), ' constraints for', \
                                out_file_name_pdb])
            make_bar_plot([averages[x+'_cst'][key_type] - averages[x][key_type] \
                    for x in sorted_reskeys], sorted_residues, "Residue Index", \
                    xlabel, title, control=wt_data[wtkey], out_dir=args.out_directory)
     
     #Second set of graphs are comparing Cis to Trans States
    for key in ['total']:
        xlabel='REU'
        if key == 'sasa':
            continue
        for tf in [False, True]:
            ckey = '_'.join([key,'C',str(tf)])
            tkey = '_'.join([key,'T',str(tf)])
            title = ' '.join(["Cis minus Trans", str(key), 'in', out_file_name_pdb, \
                    'with ligand',str(tf)])
            data_points = np.asarray([averages[x][ckey] - averages[x][tkey] \
                            for x in sorted_reskeys])
            cst_points = np.asarray([averages[x][ckey] - averages[x][tkey] \
                            for x in sorted_cstkeys])
            if not tf:
                high_mask = data_points >= 2.0
                low_mask = data_points <= -2.0
                make_bar_plot(data_points, sorted_residues, "Residue Index", xlabel,\
                title, out_dir=args.out_directory)
                if cstfiles:
                    title = ' '.join(["Cis minus Trans", str(key), 'in', \
                        out_file_name_pdb,'with ligand',str(tf),'with constraints'])
                    make_bar_plot(cst_points, sorted_residues, "Residue Index", \
                        xlabel, title, out_dir=args.out_directory)
                    
            else:
                print(ckey)
                print(tf)
                make_bar_plot(data_points, sorted_residues, "Residue Index", xlabel,\
                    title, out_dir=args.out_directory, upper_mask=high_mask, \
                    lower_mask=low_mask, legend_info={'red':'nolig cis preferred', \
                    'green':'nolig trans preferred','blue':'nolig no preference'})

                if cstfiles:
                    title = ' '.join(["Cis minus Trans", str(key), 'in', \
                        out_file_name_pdb,'with ligand',str(tf),'with constraints'])
                    make_bar_plot(cst_points, sorted_residues, "Residue Index", \
                        xlabel, title, out_dir=args.out_directory, \
                        upper_mask=high_mask, lower_mask=low_mask, \
                        legend_info={'red':'nolig cis preferred', \
                        'green':'nolig trans preferred','blue':'nolig no preference'})

    key = 'interE'
    xlabel='REU'
    ckey = '_'.join([key,'C'])
    tkey = '_'.join([key,'T'])
    title = ' '.join(["Cis vs Trans", str('interE'), 'in', \
                    out_file_name_pdb, 'with ligand',str(tf)])
    make_bar_plot([averages[x][ckey] - averages[x][tkey] for x in sorted_reskeys], \
                sorted_residues, "Residue Index", xlabel, title, \
                out_dir=args.out_directory)

    if cstfiles:
        title = ' '.join(["Cis vs Trans", str('interE'), 'in', \
                    out_file_name_pdb, 'with ligand',str(tf),'and csts'])
        make_bar_plot([averages[x][ckey] - averages[x][tkey] for x in sorted_cstkeys],\
                    sorted_residues, "Residue Index", xlabel, title, \
                    out_dir=args.out_directory)
    
    #parser.add_argument('-lim', '--value_limits', type=float, required=False, \
    #Finding good data for generating scatterplots of apo vs bound states
    
    print_out("Setting up data for the scatterplot")
    partial_keysets = ['_'.join(x.split('_')[0:2]) for x in keysets if len(x.split('_')) >= 2]
    print_out("Set of keysets")
    print(keysets.keys())
    print_out("WT data keys")
    print(wt_data.keys())
    for lines in keysets:
        if 'inter' in lines:
            partial_keysets.append(lines)
    
    partial_keysets = list(set(partial_keysets))
    for tf in ['True', 'False']:
        remove_residues = []
        for resi in sorted_reskeys:
#ass    igning values to necessary checks
            tot_T_f = averages[resi]['total_T_'+tf]
            tot_C_f = averages[resi]['total_C_'+tf]
            wt_f = wt_data['total_'+tf]
            if abs(tot_T_f - wt_f) > args.value_limits and \
               abs(tot_C_f - wt_f) > args.value_limits:
                remove_residues.append(resi)

#Check to remove models where both models are bad in apo state
    # sorted_cstkeys, sorted_reskeys, sorted_residues
    print_out('removing these residues:')
    print_out(remove_residues)
    for resi in set(remove_residues):
        sorted_residues.remove(int(resi[:-1]))
        sorted_cstkeys.remove(resi+'_cst')
        sorted_reskeys.remove(resi)
    print(sorted_residues)
    print(len(sorted_residues))
    print(sorted_reskeys)
    print(len(sorted_reskeys))
    print(sorted_cstkeys)
    print(len(sorted_cstkeys))

    x20=[-20,20]
    x25=[-5,20]
    x5=[-5,5]
    x0=[0,0]
    apobound = {}
    print(sorted_residues)
    print(partial_keysets)
    print(averages['22A'].keys())
    for key_type in partial_keysets:
        tf = [key_type.split('_')[0] + "_" + x for x in ['False', 'True']]
        if 'inter' in key_type or 'sasa' in key_type:
            continue
        print_out('wt key - '+ str(tf))
        key_false = key_type + '_False'
        key_true = key_type + '_True'
        apo = [averages[x][key_false] - wt_data[tf[0]] for x in sorted_reskeys]
        apobound[key_false] = apo
        bound = [averages[x][key_true] - wt_data[tf[1]] for x in sorted_reskeys]
        apobound[key_true] = bound
        title = key_type + " apo vs bound"            
        make_scatterplot(np.asarray(apo), np.asarray(bound), 'apo REU', 'bound REU', \
                    title, outdir=args.out_directory, tick_labels=sorted_reskeys,\
                    xylimits=[x25,x25], lines=[[x20,x20]])
        if cstfiles:
            apo = [averages[x][key_false] - wt_data[tf[0]] for x in sorted_cstkeys]
            apobound['cst_'+key_false] = apo
            bound = [averages[x][key_true] - wt_data[tf[1]] for x in sorted_cstkeys]
            apobound['cst_'+key_true] = bound
            title = key_type + " apo vs bound - with constraints"            
            make_scatterplot(np.asarray(apo), np.asarray(bound), 'apo REU', \
                    'bound REU', title, outdir=args.out_directory, \
                    tick_labels=sorted_reskeys, xylimits=[x25,x25], lines=[[x20,x20]])

    
    #apobound is normalized data apporpriate to the wt_data
    #major_trans = [apobound['total_T_True'][resi] - apobound['total_C_False'][resi] \
    #            for resi in range(0,len(apobound['total_T_True']))]
    #major_cis = [apobound['total_C_True'][resi] - apobound['total_T_False'][resi] \
    #            for resi in range(0,len(apobound['total_C_True']))]
    #title = 'Cis vs Trans states-- bound-apo energies'            
    spaces=' '*10 + 'vs' + ' '*10
    #make_scatterplot(np.asarray(major_trans), np.asarray(major_cis), \
    #    'Trans bound' + spaces + 'Cis apo', 'Cis bound' + spaces + 'Trans apo', title, \
    #    outdir=args.out_directory, tick_labels=sorted_reskeys, \
    #    xylimits=[x20,x20], lines=[[x20,x0],[x0,x20]])
    
    apo = [apobound['total_C_False'][resi] - apobound['total_T_False'][resi] \
                for resi in range(0,len(apobound['total_C_False']))]
    bound = [apobound['total_C_True'][resi] - apobound['total_T_True'][resi] \
                for resi in range(0,len(apobound['total_C_True']))]
    title = 'Cis minus Trans for apo vs bound'            
    make_scatterplot(np.asarray(apo), np.asarray(bound), 'Cis apo' + spaces + \
            'Trans apo', 'Cis bound' + spaces + 'Trans bound', title, \
            outdir=args.out_directory, tick_labels=sorted_residues, \
            xylimits=[x20,x20], lines=[[x20,x0]])
    title = 'Cis minus Trans for apo vs bound - zoom'            
    make_scatterplot(np.asarray(apo), np.asarray(bound), 'Cis apo' + spaces + \
            'Trans apo', 'Cis bound' + spaces + 'Trans bound', title, \
            outdir=args.out_directory, tick_labels=sorted_residues, \
            xylimits=[x5,x5], lines=[[x5,x0]])
    if cstfiles:
        apo = [apobound['cst_total_C_False'][resi] - \
                    apobound['cst_total_T_False'][resi] \
                    for resi in range(0,len(apobound['cst_total_C_False']))]
        bound = [apobound['cst_total_C_True'][resi] - \
                    apobound['cst_total_T_True'][resi] \
                    for resi in range(0,len(apobound['cst_total_C_True']))]
        title = 'Cis minus Trans for apo vs bound with csts'     
        make_scatterplot(np.asarray(apo), np.asarray(bound), 'Cis apo' + spaces + \
                'Trans apo', 'Cis bound' + spaces + 'Trans bound', title, \
                outdir=args.out_directory, tick_labels=sorted_residues, \
                xylimits=[x20,x20], lines=[[x20,x0]])
        title = 'Cis minus Trans for apo vs bound with csts - zoom'  
        make_scatterplot(np.asarray(apo), np.asarray(bound), 'Cis apo' + spaces + \
                'Trans apo', 'Cis bound' + spaces + 'Trans bound', title, \
                outdir=args.out_directory, tick_labels=sorted_residues, \
                xylimits=[x5,x5], lines=[[x5,x0]])
    
   
    
   
if __name__ == '__main__':
    args = parse_args()
    main(args)

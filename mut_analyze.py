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

#Script starting from here
def parse_args():
    info = """
        This script should take a protein file and a few other parameters
        and output a series of files with a point mutation and relaxed proteins.
        """
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-p', '--wt_pdb', type=str, required=True,
                        help='Relaxed starting PDB file which the input files \
                        are based upon. Should be relaxed with no alterations.')
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
    parser.add_argument('-enz', '--enzdes_cstfile', type=str, required=False,
                        help='Add for implementation of constraints used in the \
                        design of the protein. Run without for unconstrained scores')
    parser.add_argument('-cis', '--cis_params_file', type=str, \
                        required=any(x in ['--uaa', '--unnatural'] for x in sys.argv), \
                        help='Cis state params files for unnatural amino acid inclusions')
    parser.add_argument('-trans', '--trans_params_file', type=str, \
                        required=any(x in ['--uaa', '--unnatural'] for x in sys.argv),\
                        help='Trans state params files for unnatural amino acid inclusions')
    parser.add_argument('-pre', '--previously_run_files', type=str, required=False, \
                        default=None, help="Load data saved from a previous run.")
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

def dictionaried_poses(pdb_file_names, score_function):
    pose_set = {}
    for pdb_file_name in pdb_file_names:
        pose_set[pdb_file_name] = pr.pose_from_pdb(pdb_file_name)
        score_function(pose_set[pdb_file_name])
    return pose_set

def scoring_poses(pose_dict, info_dict, score_function):
    scores = {}
    for pose in pose_dict.keys():
        current_scores = {}
        current_scores['total_' + info_dict[pose]['cisortrans'] + '_' + \
            str(info_dict[pose]['contain_ligand'])] \
            = total_energy(pose_dict[pose], score_function)
        current_scores['sasa_' + info_dict[pose]['cisortrans'] + '_' + \
            str(info_dict[pose]['contain_ligand'])] \
            = total_sasa(pose_dict[pose])
        scores[pose] = current_scores
    return scores

def interaction_energy_addition(pose, lig_res_selector, sf, add_to_dict, cisortrans):
    not_lig = NotResidueSelector()
    not_lig.set_residue_selector(lig_res_selector)
    interact_metric = InteractionEnergyMetric()
    interact_metric.set_scorefunction(sf)
    interact_metric.set_residue_selectors(lig_res_selector, not_lig)
    add_to_dict['interE_'+cisortrans] = interact_metric.calculate(pose)

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
    out_file = outdir + filename + '.png'
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
    out_file = out_dir + title.replace(' ', '_') + '.png'
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close('all')
    print_out("Writing png file: " + out_file)

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
    
    init_args = ' '.join(['-extra_res_fa'] + params)
    

    if not args.previously_run_files:
    #Starting up rosetta with the appropriate params files -
    #-- need both unnatural aa and ligand file to be added (if the ligand isn't a protein).
        print_out("emd182::Starting rosetta with the following parameters: " + init_args)
        
        pr.init(init_args)
        #Setting up score function
        sf = pr.rosetta.core.scoring.ScoreFunction()
        sf.add_weights_from_file('ref2015')
        
        #Looking in the directory to determine what to group:
        pdb_list = glob('/'.join([args.source_directory, '*.pdb']))
    
        #Briefly setting up the wt poses
        wt_pose = {}
        for pdb_file_name in args.wt_pdb.split(','):
            print(pdb_file_name)
            pdb_list.remove(pdb_file_name)
            if 'yeslig' in pdb_file_name:
                wt_pose['yeslig'] = pr.pose_from_pdb(pdb_file_name)
            elif 'nolig' in pdb_file_name:
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
    
        pdb_dict = {} #Dictionary of pdb info - name, resi, uaa, cis/trans, ligand, nstruct
        group_set = {} #Dictionary of [residue number] : pdb file name for analysis
        
        print_out("Interpreting the list of pdbs in the directory: " + args.source_directory )
        #Sorting file names into groups and into dictionaries
        for pdb_file_name in pdb_list:
            pdb_dict[pdb_file_name] = file_interpreter_to_dict(pdb_file_name)
            resi = pdb_file_name.split('/')[-1].strip('.pdb').split('_')[1]
            if_list_key_not_in_dict(resi, group_set)
            group_set[resi].append(pdb_file_name)
        
        scored_set = {} #Dictionary of scores at each residue index
        print_out('Taking groups of proteins and scoring them one by one.')
        #Taking protein groups and scoring them one by one.
        #sys.exit()
        for residue_key in group_set.keys():
            print_out("Now on residue number: "+str(residue_key))
            #Place poses in a dictionary with the key designated as the pdb_file_name.
            #Each item in a dictionary is a scored pose object.
            #Creates a new dictionary each time to allow for score analysis each round.
            resi_poses = dictionaried_poses(group_set[residue_key], sf)
            scored_set[residue_key] = scoring_poses(resi_poses, pdb_dict, sf)
            for filename in resi_poses.keys():
                if pdb_dict[filename]['contain_ligand']:
                    #print_out('may take a shit if things break here')
                    interaction_energy_addition(resi_poses[filename], ligand, sf, \
                            scored_set[residue_key][filename], pdb_dict[filename]['cisortrans'])
        #Generating wt_data for the two wt poses.
        wt_data = {}
        wt_data['total_True'] = total_energy(wt_pose['yeslig'], sf)
        wt_data['sasa_True'] = total_sasa(wt_pose['yeslig'])
        interaction_energy_addition(wt_pose['yeslig'], ligand, sf, wt_data, 'True')
        wt_data['total_False'] = total_energy(wt_pose['nolig'], sf)
        wt_data['sasa_False'] = total_sasa(wt_pose['nolig'])

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
                -- this dictionary contains 3 types of scores,
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
    
    averages = {}
    for residue,dictionary in collated.items():
        averages[residue] = {}
        for k, v in dictionary.items():
            averages[residue][k] = sum(v)/len(v)
    
    sorted_resi_ints = sorted([int(x[0:-1]) for x in averages.keys()]) 
    sorted_residues = []
    for ints in sorted_resi_ints:
        for key in averages.keys():
            if ints == int(key[0:-1]):
                sorted_residues.append(key)

    #should be a list of numbers
#Keytypes should be - total_C_True total_T_True sasa_C_True sasa_T_True
#Keytypes should be - total_C_False total_T_False sasa_C_False sasa_T_False
#Keytypes should be - interE_C interE_T
    out_file_name_pdb = args.wt_pdb.split("/")[-1].split('_')[0]
    
    #First set of graphs are normalized against the wt.
    for key_type in keysets:
        keybreak = key_type.split('_')
        wtkey = ''.join([keybreak[0],"_",keybreak[-1]])
        if keybreak[0] == 'interE':
            wtkey = ''.join([keybreak[0],"_True"])
            #wtkey = ''.join([keybreak[0],"_F"])
        title = ' '.join(["Normalized", str(key_type), 'across', out_file_name_pdb])
        xlabel='REU'
        if 'sasa' in key_type:
            xlabel='SASA - A^2'
        make_bar_plot([averages[x][key_type] for x in sorted_residues], \
                sorted_residues, "Residue Index", xlabel, title, \
                control=wt_data[wtkey], out_dir=args.out_directory)
     
    #for x in sorted_residues:
    #    if len(averages[x].keys()) != 10:
    #        print(x)
    #sys.exit()
     #Second set of graphs are comparing Cis to Trans States
    for key in ['sasa', 'total']:
        xlabel='REU'
        if key == 'sasa':
            xlabel='SASA - A^2'
        for tf in [False, True]:
            ckey = '_'.join([key,'C',str(tf)])
            tkey = '_'.join([key,'T',str(tf)])
            title = ' '.join(["Cis minus Trans", str(key), 'in', out_file_name_pdb, \
                    'with ligand',str(tf)])
            data_points = np.asarray([averages[x][ckey] - averages[x][tkey] \
                            for x in sorted_residues])
            if not tf:
                high_mask = data_points >= 2.0
                low_mask = data_points <= -2.0
                make_bar_plot(data_points, sorted_residues, "Residue Index", xlabel, title, \
                    out_dir=args.out_directory)
            else:
                print(ckey)
                print(tf)
                make_bar_plot(data_points, sorted_residues, "Residue Index", xlabel, title, \
                    out_dir=args.out_directory, upper_mask=high_mask, lower_mask=low_mask, \
                    legend_info={'red':'nolig cis preferred', 'green':'nolig trans preferred',\
                            'blue':'nolig no preference'})

    key = 'interE'
    xlabel='REU'
    ckey = '_'.join([key,'C'])
    tkey = '_'.join([key,'T'])
    title = ' '.join(["Cis vs Trans", str('interE'), 'in', \
                    out_file_name_pdb, 'with ligand',str(tf)])
    make_bar_plot([averages[x][ckey] - averages[x][tkey] for x in sorted_residues], \
                sorted_residues, "Residue Index", xlabel, title, out_dir=args.out_directory)
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
        for resi in sorted_residues:
#ass    igning values to necessary checks
            tot_T_f = averages[resi]['total_T_'+tf]
            tot_C_f = averages[resi]['total_C_'+tf]
            wt_f = wt_data['total_'+tf]
            if abs(tot_T_f - wt_f) > args.value_limits and \
               abs(tot_C_f - wt_f) > args.value_limits:
                remove_residues.append(resi)
#Che    ck to remove models where both models are bad in apo state
        
        print_out('removing these residues:')
        print_out(remove_residues)
        for resi in set(remove_residues):
            sorted_residues.remove(resi)

    x20=[-20,20]
    x25=[-5,20]
    x5=[-5,5]
    x0=[0,0]
    apobound = {}
    print(sorted_residues)
    for key_type in partial_keysets:
        tf = [key_type.split('_')[0] + "_" + x for x in ['False', 'True']]
        if 'inter' in key_type or 'sasa' in key_type:
            continue
        print_out('wt key - '+ str(tf))
        key_false = key_type + '_False'
        key_true = key_type + '_True'
        apo = [averages[x][key_false] - wt_data[tf[0]] for x in sorted_residues]
        apobound[key_false] = apo
        bound = [averages[x][key_true] - wt_data[tf[1]] for x in sorted_residues]
        apobound[key_true] = bound
        title = key_type + " apo vs bound"            
        make_scatterplot(np.asarray(apo), np.asarray(bound), 'apo REU', 'bound REU', \
                    title, outdir=args.out_directory, tick_labels=sorted_residues,\
                    xylimits=[x25,x25], lines=[[x20,x20]])
    
    #apobound is normalized data apporpriate to the wt_data
    major_trans = [apobound['total_T_True'][resi] - apobound['total_C_False'][resi] \
                for resi in range(0,len(apobound['total_T_True']))]
    major_cis = [apobound['total_C_True'][resi] - apobound['total_T_False'][resi] \
                for resi in range(0,len(apobound['total_C_True']))]
    title = 'Cis vs Trans states-- bound-apo energies'            
    spaces=' '*10 + 'vs' + ' '*10
    make_scatterplot(np.asarray(major_trans), np.asarray(major_cis), \
        'Trans bound' + spaces + 'Cis apo', 'Cis bound' + spaces + 'Trans apo', title, \
        outdir=args.out_directory, tick_labels=sorted_residues, \
        xylimits=[x20,x20], lines=[[x20,x0],[x0,x20]])
    
    apo = [apobound['total_C_False'][resi] - apobound['total_T_False'][resi] \
                for resi in range(0,len(apobound['total_C_False']))]
    bound = [apobound['total_C_True'][resi] - apobound['total_T_True'][resi] \
                for resi in range(0,len(apobound['total_C_True']))]
    title = 'Cis minus Trans for apo versus bound states'            
    make_scatterplot(np.asarray(apo), np.asarray(bound), 'Cis apo' + spaces + 'Trans apo', \
            'Cis bound' + spaces + 'Trans bound', title, outdir=args.out_directory, \
            tick_labels=sorted_residues, xylimits=[x20,x20], lines=[[x20,x0]])
    title = 'Cis minus Trans for apo versus bound states - zoom'            
    make_scatterplot(np.asarray(apo), np.asarray(bound), 'Cis apo' + spaces + 'Trans apo', \
            'Cis bound' + spaces + 'Trans bound', title, outdir=args.out_directory, \
            tick_labels=sorted_residues, xylimits=[x5,x5], lines=[[x5,x0]])
    
   
if __name__ == '__main__':
    args = parse_args()
    main(args)

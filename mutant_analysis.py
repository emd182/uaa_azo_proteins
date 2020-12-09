# This script is for analyzing mutant proteins with unnatural amino acids included
# Editor: Elliott Dolan, Khare Lab, Rutgers University, 2020

#!/usr/bin/python

#Library importing
###BE SURE TO CLEAN OUT WHEN DONE###
import pyrosetta as pr
import multiprocessing as mp
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import sys, os, itertools, ntpath, math, argparse, pickle
from glob import glob
from pyrosetta.rosetta.protocols.constraint_generator import \
    AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.core.select.residue_selector import \
    NeighborhoodResidueSelector, ResidueIndexSelector, OrResidueSelector,\
    ResidueNameSelector, NotResidueSelector, InterGroupInterfaceByVectorSelector,\
    AndResidueSelector
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.protocols.enzdes import AddOrRemoveMatchCsts
from pyrosetta.rosetta.core.pack.task.operation import \
    IncludeCurrent, ExtraRotamers, OperateOnResidueSubset, \
    RestrictToRepackingRLT, PreventRepackingRLT
from pyrosetta.rosetta.protocols.grafting.simple_movers import \
    DeleteRegionMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.simple_metrics.metrics import \
    TotalEnergyMetric, RMSDMetric, InteractionEnergyMetric, \
    SasaMetric, SequenceMetric, SequenceSimilarityMetric
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import \
    PerResidueEnergyMetric
from pyrosetta.rosetta.core.scoring import ScoreType

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
    parser.add_argument('-dir', '--source_directory', type=str, required=True, \
                        help='Directory where the files that will be compared to\
                        will be found. Should contain mutants from the program \
                        mut_select.py - containing proteins in both isomers and \
                        with and without the bound ligand.')
    parser.add_argument('-part', '--partial_set', type=str, required=False, \
                        help='Instead of selecting a full set for analysis,\
                        select a partial set for analysis for parallel analysis.\
                        Select a component of the files for analysis - no/yeslig, \
                        residue id, uaa selection, etc')
    parser.add_argument('-l', '--ligand_type', type=str, required=True,
                        choices=['protein','ligand'],
                        help='Ligand type for identification - either 3-letter\
                        ligand name.')
    parser.add_argument('-res', '--residue_set', type=str, 
                        required=any(x in ['protein'] for x in sys.argv),
                        help='Comma delimited residue list for a peptide sequence \
                        behaving as a ligand in the mutation detection script.')
    parser.add_argument('-lig', '--ligand_params', type=str, default=None,
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
    parser.add_argument('-pkl', '--pickled_file', type=str, required=False, \
                        default=None, help="Load data saved from a single previous run.")
    parser.add_argument('-pkls', '--pickled_directory', type=str, required=False, \
                        default=None, help="Load data saved from multiple runs.")
    parser.add_argument('-per', '--per_resi_graph', action='store_true', default=False, \
                        help="Make graphs for each residue for all proteins.")
    parser.add_argument('-filt','--filter', action='store_true', default=False, \
                        help='Call --filter to begin filtering the proteins.')
    parser.add_argument('-lim', '--filter_limits', type=float, required=False, \
                        default=15.0, help="Add limits of total scores in the case the\
                        apo states are diverge too greatly from the wildtype.")
    parser.add_argument('-fval', '--filter_on_E', type=str, required=False, \
                        default='delta_total_E,nolig', help="Sets energy value to \
                        apply the above filter to. Default is delta_total_E. Add comma\
                        yeslig / nolig to apply filter to either bound or apo state")
    parser.add_argument('-solo', '--analysis_only', action='store_true', required=False,\
                        default=False, help='Include this flag to only run the analysis \
                        on the data and not attempt to analyze the data.')
    parser.add_argument('-out', '--out_file', type=str, required=False, \
                        help='Outputs a pickled set of information for later analysis.\
                        Or to prevent having to repeat the analysis.')
    parser.add_argument('-od', '--out_directory', type=str, required=False, \
                        help='Output Directory fo all output files.')

    return parser.parse_args()

"""
##############################################################################
Either general use or currently Unused
##############################################################################
"""

def make_histogram(data, bins, xlabel, ylabel, title):
    print_out("not written yet")
    y_values = eval(data)
    x = list(range(len(y_values)))
    plt.bar(x, y_values, width=1, color='blue')
    plt.xticks(x, x_ticks)

def add_to_legend(legend_list, colo, label):
    legend_list.append(mpatches.Patch(color=colo, label=label))

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

def interaction_energy_addition(pose, lig_res_selector, sf, add_to_dict, cisortrans):
    not_lig = NotResidueSelector()
    not_lig.set_residue_selector(lig_res_selector)
    interact_metric = InteractionEnergyMetric()
    interact_metric.set_scorefunction(sf)
    interact_metric.set_residue_selectors(lig_res_selector, not_lig)
    add_to_dict['interE_'+cisortrans] = interact_metric.calculate(pose)

"""
##############################################################################
Classes for group pdb analysis. Sorting for pdb must happen before.
##############################################################################
"""

class single_mutant_analysis():
    """
        Data storage object for a mutant model. The input pdbs need to be saved
        as the asymmetric model if symmetrized, and the constraint files need to
        be the correct set - no ligand constraints for no ligands, ligand for 
        ligand constraints. This class doesn't analyze the mutants based on 
        file name, so that understanding must be performed elsewhere.
    """
    def __init__(self, pdb_file_name, pdb_pose, score_function, \
                    wild_type_pose=None, cst_file_name=None, sym=None):
        #OK, apparently pickle hates pickling rosetta objects.
        #Writing a function to move all the information extracted from
        #These Classes into a Dataframe.

        self.pdb_name = ntpath.basename(pdb_file_name)
        self.symmetry = sym
        self.cst_file = cst_file_name
        self.sf_default = pr.rosetta.core.scoring.deep_copy(score_function)

        #Rosetta Utilities that shouldn't be remade every time
        #A new function is called:
        self.ris = ResidueIndexSelector()
        self.prEm = PerResidueEnergyMetric()
        
        #Input PDBS need to only have asymmetric unit and the 
        #symmdef file for this class
        #That includes the wildtype pdbs
        #Isolate the constraint energies from the bulk energies 
        self.mutant = pdb_pose
        self.wt = pr.Pose()
        self.wt.assign(wild_type_pose)

        #Energies and energy types
        self.total_energy = 0
        self.cstE = 0 
        self.per_res_total_E = {}
        self.per_res_total_cstE = {}
        
        self.delta_total_E = 0
        self.delta_cstE = 0
        self.delta_per_res_E = {}
        self.delta_per_res_cstE = {}

        #Mutation info before symmetry
        self.sequence = ''
        self.resi_order = {}
        self.mutated_residues = {}
        self.calculate_sequence()
        
        if self.symmetry:
            self.symmetric = False
            self.apply_symmetry()

        if self.cst_file:
            self.apply_match_constraints()
        
        self.calculate_energies()
        
        self.energies = {'total_energy': self.total_energy,\
                         'cstE': self.cstE,\
                         'per_res_total_E': self.per_res_total_E,\
                         'per_res_total_cstE': self.per_res_total_cstE,\
                         'delta_total_E': self.delta_total_E,\
                         'delta_cstE': self.delta_cstE,\
                         'delta_per_res_E': self.delta_per_res_E,\
                         'delta_per_res_cstE': self.delta_per_res_cstE, \
                         'resi_order':self.resi_order }

    def apply_match_constraints(self):
        """
        Simply adds loaded constraints to the pose and wt pose.
        """
        AoRMcsts = AddOrRemoveMatchCsts()
        AoRMcsts.set_cst_action(pr.rosetta.protocols.enzdes.CstAction.ADD_NEW)
        print(self.pdb_name)
        AoRMcsts.cstfile(self.cst_file)
        AoRMcsts.apply(self.mutant)
        AoRMcsts.apply(self.wt)
    
    def apply_symmetry(self):
        """
        Will attempt to create a symmetric pose from an asymmetric pose.
        Should only require a pose to be symmetric prior to being loaded,
            self.symmetric = True
        and then this function converts that regular pose into a symmetrized pose. 
        """
        pr.rosetta.core.pose.symmetry.make_symmetric_pose(self.mutant, self.symmetry) 
        assert pr.rosetta.core.pose.symmetry.is_symmetric(self.mutant) == True, \
            "Mutant pose of "+ self.pdb_name + " is not symmetric - Test 1"
        if self.wt is not None:
            pr.rosetta.core.pose.symmetry.make_symmetric_pose(self.wt, self.symmetry)
            assert pr.rosetta.core.pose.symmetry.is_symmetric(self.mutant) == True, \
                "Wildtype pose for "+ self.pdb_name + " is not symmetric - Test 1"
    
    def calculate_sequence(self):
        """
        Calculate sequence when compared to self - takes into account
        symmetry and only counts symmetric mutations once.
        *** If not considered symmetric, it will duplicate mutants ***
        """
        sm = SequenceMetric()
        self.sequence = sm.calculate(self.mutant)
        
        assert self.mutant.total_residue() == self.wt.total_residue()
        
        #checking between two sequences: wt and original design
        for i in range(1, self.mutant.total_residue() + 1):
            pdb_number = ''.join(self.wt.pdb_info().pose2pdb(i).split())
            self.resi_order[i] = pdb_number
            self.ris.set_index(pdb_number)
            sm.set_residue_selector(self.ris)
            final_aa = sm.calculate(self.mutant)
            orig_aa = sm.calculate(self.wt)
            
            if final_aa != orig_aa:
                pdb_number = self.wt.pdb_info().pose2pdb(i).split()[0]
                self.mutated_residues[pdb_number] = [orig_aa, final_aa]

    def calculate_energies(self):
        """
        Calculate all energies of the poses, compared to the wt.
        """
        self.sf_default(self.mutant)
        no_cst_sf = pr.rosetta.core.scoring.deep_copy(self.sf_default)
        self.set_constraint_scores_to_x(no_cst_sf, 0.0)
        #Unconstrained energies
        tem = TotalEnergyMetric()
        tem.set_scorefunction(no_cst_sf)
        self.total_energy = tem.calculate(self.mutant)
        self.per_res_total_E = self.calculate_per_res_energies( \
                                        self.mutant, no_cst_sf)
        if self.wt is not None:
            #Calculating total Energy Change compared to WT
            self.delta_total_E = self.total_energy - tem.calculate(self.wt)
            self.delta_per_res_E = self.calculate_delta_per_res_energies( \
                                self.per_res_total_E, no_cst_sf)

        #Now determine the constraint component of the poses
        if self.cst_file is not None:
            cst_only_sf = pr.rosetta.core.scoring.ScoreFunction()
            self.set_constraint_scores_to_x(cst_only_sf, 1.0)
            tem.set_scorefunction(cst_only_sf)
            self.cstE = tem.calculate(self.mutant) 
            self.per_res_total_cstE = self.calculate_per_res_energies( \
                                self.mutant, cst_only_sf)  
            if self.wt is not None:
                self.delta_cstE = self.cstE - tem.calculate(self.wt)
                self.delta_per_res_cstE = self.calculate_delta_per_res_energies( \
                                self.per_res_total_cstE, cst_only_sf)
        
    def calculate_per_res_energies(self, curr_pose, score_function):
        """ Calculate energies on a per-residue basis """
        self.prEm.set_scorefunction(score_function)
        per_res_dict = {}
        per_res_Es = self.prEm.calculate(curr_pose)
        for res in range(1, curr_pose.total_residue() + 1 ):
            resi = ''.join(curr_pose.pdb_info().pose2pdb(res).split())
            per_res_dict[resi] = per_res_Es[res]
        return per_res_dict
    
    def calculate_delta_per_res_energies(self, mut_per_res, score_function):
        wt_per_res = self.calculate_per_res_energies(self.wt, score_function)
        delta_per_res_dict = {}
        for key in wt_per_res.keys():
            delta_per_res_dict[key] = mut_per_res[key] - wt_per_res[key]
        #print(delta_per_res_dict)
        return delta_per_res_dict
            

    def set_constraint_scores_to_x(self, score_function, value):
        """ Sets constraint scores of atom_pair_constraint, angle constraint,
        dihedral_constraint, and coordinate_constraint to zero 
        """
        score_function.set_weight(ScoreType.atom_pair_constraint, value)
        score_function.set_weight(ScoreType.angle_constraint, value)
        score_function.set_weight(ScoreType.dihedral_constraint, value)
        score_function.set_weight(ScoreType.coordinate_constraint, value)

class point_mutant_set_analysis():
    """
    Takes a list of file names, makes them into a pose set and
    organizes those pose objects in order to make them managable.
    pdblist == list of pdbfiles -- score_function == initialized score function
    mutant_index == residue index for the main mutant residue
    """
    def __init__(self, pdblist, mutant_index, score_function, \
                        wt_pose=None, cstfile=None, sym=None):
        
        self.mut_resi = mutant_index
        self.pdb_list = pdblist
        self.total = len(pdblist)
        self.sf = score_function
        
        #Optional variables
        self.wt_pose = wt_pose
        self.cst_file = cstfile
        self.sym = sym

        self.poses = {}
        self.pose_setup()
        
        self.analyzed_poses = {}
        self.pose_analyze()
        
    def pose_setup(self):
        """ 
        Takes a pdblist and score function and generates a dictionary
        of scored and analyzed poses 
        """
        for pdb_file_name in self.pdb_list:
            self.poses[pdb_file_name] = pr.pose_from_pdb(pdb_file_name)
    
    def pose_analyze(self):
        """ Load the poses into the single_mutant_analysis class """
        for key, pose in self.poses.items():
            self.analyzed_poses[key] = single_mutant_analysis( key, pose, 
                        self.sf, wild_type_pose=self.wt_pose, \
                        cst_file_name=self.cst_file, sym=self.sym)                        

"""
##############################################################################
Takes a dataframe of the above point_mutant_set_analysis classes as an input 
and converts it into a Dataframe without rosetta objects.
##############################################################################
"""

def mutant_dataframe_collapse(input_DF):
    """
    Setting up a dataframe that will scroll through every 
    point_mutant_set_analysis class in the dataframe and dump the data into a
    dataframe such that it can be saved and more accessible.
    """
    DF_columns = [x for x in input_DF.columns]
    datalist = {col:[] for col in DF_columns}
    for key in DF_columns[:-1]:
        datalist[key] = [val for val in input_DF[key]]
    fcol = DF_columns[-1]
    for index, value in input_DF[fcol].iteritems():
        print('Converting ' + str(index))
        final_column = {}
        for pdb_key, class_data in value.analyzed_poses.items():
            final_column[pdb_key] = {key:class_data.energies[key] for \
                        key in class_data.energies.keys()}
        datalist[fcol].append(final_column)
    return pd.DataFrame(data=datalist, columns=DF_columns)

"""
##############################################################################
Analysis functions to extract data from a collapsed dataframe.
##############################################################################
"""

def single_E_set_calc(frame_part, min_max_avg='min', sele_E='total_energy', \
                                pdbkey=None):
    """ 
    Returns the energies of the analyzed poses in a portion of the input dataframe.
    Input dataframe part must only constitute the same set of pdb proteins
    min_max_avg - requested to determine the min, max, or avg
    sele_E - request specific types of data to be calculated.
    """
    totals = {}
    keys = [x for x in frame_part.keys()]
    testkey = keys[np.random.randint(len(keys))]
    assert sele_E in frame_part[testkey].keys() and\
            'per' not in sele_E, "Select an appropriate key"
    
    for pdbid, item in frame_part.items():
        if type(item) == list:
            totals[pdbid] = item[sele_E]
        else:
            totals[pdbid] = [item[sele_E]]

    if pdbkey:
        assert pdbkey in totals.keys(), "PDBkey is not valid."
        return totals[pdbkey], pdbkey
    
    E_set = pd.DataFrame.from_dict(totals, orient='index')

    assert min_max_avg in ['avg', 'min', 'max'], \
            "Select an appropriate return - min, max, or avg"
    
    if min_max_avg == 'avg':
        #energy will be the average, index will be the most 'average' index
        energy = E_set.mean()
        index = abs(E_set - energy).idxmin()

    elif min_max_avg == 'min':
        #energy will be the minimum, index will be the most 'minimum' index
        energy = E_set.min()
        index = E_set.idxmin()
    
    elif min_max_avg == 'max':
        #energy will be the maximum, index will be the most 'maximum' index
        energy = E_set.max()
        index = E_set.idxmax()
    
    else:
        print("You shouldn't be here")
        sys.exit()
    return energy[0], index[0] 

def return_per_res_data(frame_part, min_max_avg='min', sele_E='per_res_total_E',\
                                                pdbkey=None):
    """
    Returns per residue information depending on what's asked -
    Returns minimum, maximum, or average
    frame_part is a dictionary of all the data for 
    one set of residue/uaa/apo-bound models
    """
    corresponding_sets = {'per_res_total_E': 'total_energy',\
                          'per_res_total_cstE': 'cstE', \
                          'delta_per_res_E':'delta_total_E', \
                          'delta_per_res_cstE':'delta_cstE'}

    assert sele_E in corresponding_sets.keys(), "Select an appropriate key and set"
    
    assert min_max_avg in ['min','max','avg'], "Select an appropriate minmaxavg set"

    if pdbkey:
        return frame_part[pdbkey][sele_E], pdbkey 
    
    energy_set, key = single_E_set_calc(frame_part, min_max_avg, \
                                            corresponding_sets[sele_E])

    return frame_part[key][sele_E], key 

def multi_E_set_calc(frame, min_max_avg='min', multi_E=['total_energy','cstE']):
    """
    Calls single_E_set_calc multiple times and returns the summed data.
    typical combos: total_energy and cstE ; delta_total_E and delta_cstE
    mma = min_max_avg
    """
    summed_E = 0.0
    keys = []
    assert isinstance(multi_E, list), \
                    "If not using a list, please use single_E_set_calc"
    for item in multi_E:
        single_E, curr_key = single_E_set_calc(frame, min_max_avg=min_max_avg, sele_E=item)
        if not keys:
            keys.append(curr_key) 
        summed_E += single_E
    return summed_E, keys[0]

def multi_per_res_data(frame_part, min_max_avg='min', multi_E=['per_res_total_E',\
                                                 'per_res_total_cstE']):
    """
    Calls return_per_res_data multiple times and returns the summed data.
    typical combos: total_energy and cstE ; delta_total_E and delta_cstE
    minmaxavg sets for minimum, maximum, and average
    """
    summed_E = {}
    keys = []
    assert isinstance(multi_E, list), \
                    "If not using a list, please use single_E_set_calc"
    for item in multi_E:
        per_res_E, key_ = return_per_res_data(frame_part, min_max_avg=min_max_avg, \
                                                sele_E=item)
        keys.append(key_)
        if summed_E:
            summed_E = {key:summed_E[key]+per_res_E[key] for key in per_res_E.keys()}
        else:
            summed_E = per_res_E

    return summed_E, keys[0]

def dataframe_extract(mut_rows, mma, datatypes, columns=['residue','point_mutant_set']):
    """ 
    Extracts the data depending on min_max_avg and datatypes input
    Returns data, residues extracted from, and directories -
    Sorts data and outputs the data as well.
    """
    extract = {}
    for idx in mut_rows.index.values.tolist():
        mut_key = mut_rows[columns[0]][idx]
        top_key = mut_key + '_top'
        if isinstance(datatypes, list):
            if 'per' in ''.join(datatypes):
                #print('running multi_per_res_data')
                extract[mut_key], extract[top_key] = multi_per_res_data( \
                        mut_rows[columns[1]][idx], min_max_avg=mma, multi_E=datatypes)
            else:
                #print('running multi_E_set_calc')
                extract[mut_key], extract[top_key] = multi_E_set_calc( \
                        mut_rows[columns[1]][idx], min_max_avg=mma, multi_E=datatypes)
        else:
            if 'per' not in datatypes:
                #print('running single_E_set_calc')
                extract[mut_key], extract[top_key] = single_E_set_calc( \
                        mut_rows[columns[1]][idx], min_max_avg=mma, sele_E=datatypes)
            else:
                #print('running return_per_res_data')
                extract[mut_key], extract[top_key] = return_per_res_data( \
                        mut_rows[columns[1]][idx], min_max_avg=mma, sele_E=datatypes)

    print("Analyzed Set: " + str(datatypes))
    sorted_labels = sorted([key for key in extract.keys() if '_' not in key],\
            key=lambda x: float(x[:-1]))
    out_data = [extract[key] for key in sorted_labels] 
    top_models = {key:[extract[key], extract[key.strip('_top')]] \
                    for key in extract.keys() if '_top' in key}
    return out_data, sorted_labels, top_models 

"""
##############################################################################
Functions for integrating the above Classes into indentificatiable DataFrame 
Objects to allow for further analysis and sorting.
##############################################################################
"""

def nameframe_pose_analyzer(pdbF, wildtype_poses, score_function, source_dir, \
                            cst_files=None, symm_def_file=None):
    """ Determines all the unique types of pdb_frame information,
        Then sending each of those new files into pdb analysis classes 
        Finally, outputs a pdbframe detailed with the following columns:
        residue, uaa, apo/bound, point_mutant_set
    """
    unique_set = {col:pdbF[col].unique() for col in pdbF.columns}
    cols = [str(key) for key in unique_set.keys()] + ['point_mutant_set']
    keys, values = zip(*unique_set.items())
    data_list = []
    for cbo in itertools.product(*values):
        print(wildtype_poses)
        if cst_files:
            cst_file = cst_files[cbo[2]]
        else:
            cst_file = cst_files
        pdb_list = pdbF.loc[ (pdbF[keys[0]] == cbo[0]) & \
                             (pdbF[keys[1]] == cbo[1]) & \
                             (pdbF[keys[2]] == cbo[2]) ]
        pdb_filenames = ['/'.join([source_dir, x]) for x in pdb_list.index]
        data_list.append([cbo[0], cbo[1], cbo[2], \
                point_mutant_set_analysis(pdb_filenames, cbo[0], score_function,\
                wt_pose=wildtype_poses[cbo[2]], cstfile=cst_file, \
                sym=symm_def_file) ] )
    pdbframe = pd.DataFrame(data_list, columns=cols) 
    return pdbframe

def DF_column_value_search(dataframe, searched_values):
    """ Searches through each column in a given DataFrame for a given value,
        Returning the column and the corresponding value as a dictionary """
    columns = []
    values = []
    for val in searched_values:
        for column in list(dataframe.columns):
            if val in dataframe[column].values:
                columns.append(column)
                values.append(val)
                break
    return columns, values

def row_select_from_column_filters(data_frame, data_search):
    """Takes a dataframe and a set of data values to search for.
        Using a set of data values, it searches each column of data for that value
        And filters out everything except for those values in that specific column.
        Returns the rows in the data_frame for the filtered values.
    """
    col_keys, val_keys = DF_column_value_search(data_frame, data_search)
    lookup = [(k,v) for k,v in zip(col_keys, val_keys)]
    filt = '{0[0]} == "{0[1]}"'.format
    request = ' & '.join(filt(q) for q in lookup)
    mut_rows = data_frame.query(request)
    return mut_rows

def filter_dataframe(unfiltered, filter_keys, filter_limits):
    """
    Takes an unfiltered dataframe and filters to the filter_limits on the
    filter_keys key. filter_key is comma deliminited.
    filter_key.split(',')[0] is the energy key.
    filter_key.split(',')[1:] is the other keys on the dataframe
    """
    col = ['residue','uaa','apobound', 'point_mutant_set']
    keys = filter_keys.split(',')
    bad_residues = []
    for resi in list(unfiltered['residue'].unique()):
        filter_ = keys[1:] + [resi]
        first_filt = row_select_from_column_filters(unfiltered, filter_)
        for idx in first_filt.index.values.tolist():
            bad_resi = 0
            for pdbs in first_filt[col[-1]][idx].keys():
                if abs(first_filt[col[-1]][idx][pdbs][keys[0]]) > filter_limits:
                    bad_resi += 1
            if bad_resi == len(first_filt[col[-1]][idx].keys()):
                bad_residues.append(resi)
    remove_residues = list( set( [x for x in bad_residues \
                                    if bad_residues.count(x) > 1 ] ) )
    if len(remove_residues) > 0:
        remove_indices = []
        for residue in remove_residues:
            filtered_frame = row_select_from_column_filters(unfiltered, [residue])
            remove_indices += filtered_frame.index.values.tolist()
        filtered = unfiltered.drop(index=remove_indices, axis=0, inplace=False)
        filt = filtered.reset_index(drop=True)
        return filt
    else:
        return unfiltered 

"""
##############################################################################
Making Per-Residue DDG Plots
##############################################################################
"""
def plot_res_ddg(ref_data, mut_data, titlen=None, filen=None):
        # Plot
    fig, axes = plt.subplots(2, 1, sharex=True)
    axes[0].plot(ref_data, linewidth=0.3, color='red')
    axes[0].plot(mut_data, linewidth=0.3, color='black')
    axes[0].set_ylabel('Residue Energy', fontsize=9)
    axes[1].plot(mut_data - ref_data, linewidth=0.3, color='black')
    axes[1].set_ylabel('Mutant - Reference', fontsize=9)
    
    if titlen != None:
        fig.suptitle(titlen, fontsize=16)

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

    if filen != None:
        plt.savefig(filen)
        print_out("Making graph for: " + filen)
    else: plt.show()
    plt.close()

def make_per_residue_graphs(frame, out_dir, cst_add=False):
    """
    Takes a dataframe of individual proteins, and makes per residue graphs
    for every model, saved with the residues name. 
            delta_per_res_dict[key] = mut_per_res[key] - wt_per_res[key]
    """
    for row, data in frame.iterrows():
        per_res_dir = '/'.join([out_dir,'per_res_'+data['residue']])
        out_directory(per_res_dir) 
        for pdb_key, pdb_dict in data['point_mutant_set'].items():
            pdb = pdb_key.split('/')[-1]
            title = pdb.replace('.pdb','') + ' per residue'
            filename = pdb.replace('.pdb','') + '_per_res'
            keysort = sorted([x for x in pdb_dict['per_res_total_E'].keys()],\
                key=lambda x: int(x[:-1]))
            total_Es = np.array([pdb_dict['per_res_total_E'][key] for key in keysort])
            delt_Es = np.array([pdb_dict['delta_per_res_E'][key] for key in keysort])
            wt_total_Es = total_Es - delt_Es
            if cst_add and pdb_dict['cstE'] != 0.0:
                print("Constraint energies != 0")
                print(pdb_dict['cstE'])
                print("Printing graphs with constraint energies added")
                title = title + ' with csts'
                filename = filename + '_csts'
                cstEs = np.array([pdb_dict['per_res_total_cstE'][k] for k in keysort])
                total_Es += cstEs
                delt_cstEs = np.array([pdb_dict['delta_per_res_cstE'][k] for k in keysort])
                wt_total_Es = total_Es - delt_cstEs - delt_Es
            out_file = per_res_dir + '/' + filename + '.png'
            plot_res_ddg(wt_total_Es, total_Es, titlen=title, filen=out_file)
            #sys.exit()

"""
##############################################################################
Making Bar Plots
##############################################################################
"""

def make_bar_plot(bar_values, bar_labels, xlabel, ylabel, title, csv_out=None, \
            ymax=10.0, ymin=-10.0, control=0.0, out_dir='./',\
            legend_info={'blue':'normalized score'}):
    
    normal = np.array(bar_values) - control
    x = np.asarray(list(range(len(normal))))
    legends = []
    plt.bar(x, normal, width=0.8, align='center', color='blue')
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
    if csv_out:
        frame = pd.DataFrame(csv_out)
        csv_dir = '/'.join([out_dir, 'csv_dir'])
        out_directory(csv_dir)
        csv_file ='/'.join([csv_dir, title.replace(' ', '_') +'.csv'])
        frame.to_csv(csv_file)
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close('all')
    print_out("Writing png file: " + out_file)
    #Add in code to dump to csv file.

def make_simple_bar_graphs(frame, dataset, uaa_info, analysis_dir, pdb_name):
    """
    Make general bar graphs for each set.
    bar graphs: delta_total_E, delta_cstE, cstE
    bar graphs for every residue set
    """
    ylabel = 'normalized REU'
    xlabel = 'Residue Index'
    
    mut_rows = row_select_from_column_filters(frame, dataset)
    
    for avgmin in ['avg', 'min']:
        for graph_type in ['delta_total_E', ['delta_total_E','delta_cstE']]:
            data, data_labels, top_models = dataframe_extract(mut_rows, avgmin, \
                    graph_type, columns=['residue','point_mutant_set'])
            data_name = str(graph_type)
            legend = {'blue':avgmin+' score'}
            title = ' '.join(['Normalized', pdb_name, ''.join(data_name), \
                    dataset[1],'with', dataset[0]+'-'+avgmin])  
            make_bar_plot(data, data_labels, xlabel, ylabel, title, \
                    csv_out=top_models, out_dir=analysis_dir, legend_info=legend) 

def make_x_minus_y_bar_graphs(frame, x_filt, y_filt, analysis_dir, pdb_name,\
            ylabel='Delta REU', xlabel='Residue Index', graph_sets=['delta_total_E',\
            ['delta_total_E','delta_cstE']], select_type=['min']):
    """
    Makes bar graphs where it filters one set from another and subtracts one 
    dataset from the other.
    """
    
    x_mut_rows = row_select_from_column_filters(frame, x_filt)
    y_mut_rows = row_select_from_column_filters(frame, y_filt)
        
    for avgmin in select_type:
        for graph_type in graph_sets:
            x_data, x_tick_names, xtop_models = dataframe_extract(x_mut_rows, avgmin,\
                graph_type, columns=['residue','point_mutant_set'])
            y_data, y_tick_names, ytop_models = dataframe_extract(y_mut_rows, avgmin,\
                graph_type, columns=['residue','point_mutant_set'])
            data_name = ' '.join(graph_type)
            legend = {'blue':'min score'}
            title = ' '.join(['Normalized', pdb_name, x_filt[1], 'vs', y_filt[1], \
                    y_filt[0]+':'+avgmin])  
            make_bar_plot(np.array(x_data), x_tick_names, xlabel, ylabel, title, \
                            out_dir=analysis_dir, legend_info=legend, \
                            control=np.array(y_data)) 

"""
##############################################################################
Making Scatter Plots
##############################################################################
"""

def base_ceil(x, base=2.5):
    return base * math.ceil(float(x)/base)
def base_floor(x, base=2.5):
    return base * math.floor(float(x)/base)

def scatter_edges(xvals, yvals, zoom=False, round=5.0):
    """
    Calculated Scatterplot edges from the xvals and yvals
    Rounds upwards or downwards to the 'round' variable to keep consistency
    zoom == True will zoom into the values within only 1 standard deviation
    """
    if zoom:
        xstd, ystd = [np.std(xvals), np.std(yvals)]
        xmean, ymean = [np.mean(xvals), np.mean(yvals)]
        return [base_floor(xmean-xstd), base_ceil(xmean+xstd)], \
               [base_floor(ymean-ystd), base_ceil(ymean+ystd)]

    else:
        xmax, xmin = [max(xvals), min(xvals)]
        ymax, ymin = [max(yvals), min(yvals)]
        xymax = max([base_ceil(xmax), base_ceil(ymax)])
        xymin = min([base_floor(xmin), base_floor(ymin)])
        return [[xymin, xymax],[xymin, xymax]]

def make_scatterplot(xvalues, yvalues, xlabel, ylabel, title, filename=None,\
        colors='b', xylimits=None, markers='o', outdir='./',\
        tick_labels=None, lines=None, xcontrol=0.0, ycontrol=0.0, zoom=False):
    xvals = xvalues - xcontrol
    yvals = yvalues - ycontrol
    if xylimits:
        xlimits, ylimits = [[-10,10],[-10,10]]
    else:
        xlimits, ylimits = scatter_edges(xvals, yvals, zoom)
    plt.scatter(xvals, yvals, s=5, c=colors, marker=markers)
    if lines:
        for line in lines:
            print(line)
            plt.plot(line[0],line[1], '-k', linewidth=1)
    axes = plt.gca()
    axes.set_xlim(xlimits)
    axes.set_ylim(ylimits)
    if tick_labels:
        for num in range(0, len(xvals)):
            plt.annotate(tick_labels[num], (xvals[num],yvals[num]), \
                textcoords="offset points",xytext=(0,5), ha='center',\
                fontsize='x-small')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    fig = plt.gcf()
    fig.set_size_inches(5,5)
    plt.show()
    if not filename:
       filename = title.replace(' ', '_')
    out_file = '/'.join([outdir,filename + '.png'])
    fig.savefig(out_file, dpi=150)#, bbox_inches='tight')
    plt.close('all')
    print_out("writing png file: " + out_file)
        
def make_xla_vs_ylb_scatter_plots(frame, x_filt, y_filt, analysis_dir, \
            x_control_filt=0.0, y_control_filt=0.0, ylabel='Delta REU', \
            xlabel='Delta REU', title='Apo vs Bound', mma='min', \
            graph_type='delta_total_E', columns=['residue','point_mutant_set'],\
            zoom_graph=False, gen_line=[[[-100,100],[-100,100]]]):
    """
    Makes bar graphs where it filters one set from another and subtracts one 
    dataset from the other.
    """
    x_mut_rows = row_select_from_column_filters(frame, x_filt)
    y_mut_rows = row_select_from_column_filters(frame, y_filt)
    x_data, x_tick_names, xtop_models = dataframe_extract(x_mut_rows, mma,\
        graph_type)
    y_data, y_tick_names, ytop_models = dataframe_extract(y_mut_rows, mma,\
        graph_type)
    if x_control_filt != 0.0:
        x_control_rows = row_select_from_column_filters(frame, x_control_filt)
        x_control_data, x_control_names, xtop_control_models = dataframe_extract( \
                x_control_rows, mma, graph_type)
        x_control_filt = np.array(x_control_data)
    if y_control_filt != 0.0:
        y_control_rows = row_select_from_column_filters(frame, y_control_filt)
        y_control_data, y_control_names, ytop_control_models = dataframe_extract( \
                y_control_rows, mma, graph_type)
        y_control_filt = np.array(y_control_data)
    
    make_scatterplot(np.array(x_data), np.array(y_data), xlabel, ylabel, \
            title, outdir=analysis_dir, tick_labels=x_tick_names, \
            lines=gen_line, xcontrol=x_control_filt, ycontrol=y_control_filt,
            zoom=zoom_graph)

"""
##############################################################################
Reading in and executing PDBs, inputs, and arguments. Before running the 
decoys into the analysis classes
##############################################################################
"""

def init_arg_parsing(ligand_type, unnatural_params, ligand_params):
    """ Taking input arguments for params files and implementing them for
    Properly running pyrosetta. """
    args = []
    print(ligand_params)
    if ligand_params != None:
        assert type(ligand_params) == str, "need a params string"
        args.append(ligand_params)
    
    for uaa in unnatural_params:
        args.append(uaa)
    
    return ' '.join(['-run:preserve_header','-extra_res_fa'] + args)

def manage_cst_files(cst_dict, cst_string):
    if ',' in cst_string:
        for cst in cst_string.split(','):
            if 'noligcst' in cst:
                cst_dict['nolig'] = cst
            else:
                cst_dict['yeslig'] = cst
    else:
        cst_dict['yeslig'] = cst_string

def symm_cst_args_manage(cst_file_string, symmdef_file_string):
   #""" Sets up additional optional arguments """
    cstfiles = {}
    if args.enzdes_cstfiles:
        manage_cst_files(cstfiles, args.enzdes_cstfiles)
    symdef = ''
    if args.symmdef_file:
        symdef = args.symmdef_file
    return cstfiles, symdef

def wild_type_poses(wt_pdb_str, all_pdbs, score_function):
    wt_pose = {}
    for pdb_file_name in wt_pdb_str.split(','):
        if pdb_file_name in all_pdbs:
            all_pdbs.remove(pdb_file_name)
        if 'yeslig' in pdb_file_name:
            wt_pose['yeslig'] = pr.pose_from_pdb(pdb_file_name)
            score_function(wt_pose['yeslig'])
        elif 'nolig' in pdb_file_name:
            wt_pose['nolig'] = pr.pose_from_pdb(pdb_file_name)
            score_function(wt_pose['nolig'])
        else:
            print('Error, submitted wild type pdbs do not have yes or no ligand' + \
                    'name. Please submit properly named pdbs.')
            sys.exit()
    return wt_pose

def ligand_interpreter(lig_type, lig_info):
    """ Takes the ligand type and residues / params information and generates
        A residue selector for that type. """ 
    if args.ligand_type == 'protein':
        ligand = ResidueIndexSelector(lig_info[0])
    elif args.ligand_type == 'ligand':
        ligand = ResidueNameSelector()
        ligand.set_residue_name3(lig_info[1].split('/')[-1].strip('.params'))
    return ligand
    
def pdb_filename_data_interpreter(pdb_list):
    """ Takes in a list of pdbs containing 3 bits of information in each filename:
        UAA, mutant residue, bound or unbound ligand """
    col = ['residue','uaa','apobound']
    nameframe = pd.DataFrame(columns=col)
    for filename in pdb_list:
        pdbname = filename.split('/')[-1]
        data = pdbname.split('_')
        nameframe.loc[pdbname] = [data[1], data[2], data[3]]
    return nameframe

def dumping_dataframe_to_pickle(out_frame, out_file, out_part=None, out_dir=None):
    """ Takes the out_frame and dumps it into a pickle file into the
        directory out_dir. Takes the out_file name as a template and 
        modifies it with out_part - if an out_part is called. """
    out = out_file
    if out_part:
        out = out_file.split('.')[0] + '_' + out_part + '.pickle'
    if out_dir:
        out_directory(out_dir)
        out = '/'.join([out_dir,out])
    print_out("Storing all the data into a pickle file")
    with open(out, 'wb') as savedfile:
        pickle.dump(out_frame, savedfile)
        print_out("Pickle has stored dataframe into " + str(out))

def pulling_dataframe_from_pickles(pickled_file=None, pickled_directory=None):
    """ Takes a pickled file or directory of pickled files as in put and
        loads the pickled information into one dataframe."""
    if pickled_directory:
        pkl_list = glob('/'.join([pickled_directory, '*.pickle']))
        print('pickle list')
    if pickled_file:
        pkl_list = [pickled_file]
    frame_list = []
    for pkl in pkl_list:
        with open(pkl, 'rb') as pkl_file:
            frame_list.append(pickle.load(pkl_file))
            print('Loading from pickle: ' + pkl)
    out_frame = pd.concat(frame_list)
    out = out_frame.reset_index(drop=True)
    return out
    #return pd.concat(frame_list)

"""
##############################################################################
Main Function
This function takes in a collection of unnatural amino acid point-mutant 
decoys of ligand-binding proteins (Output from mutation_design.py), organizes
them, and analyzes them in apo/bound states and cis/trans states. 
The analyzed data is then extracted into Dataframes that are pickled for fast
retreval. Analyzed data is output into bar graphs and scatter plots.
Additionally, the data is dumped into csv files for decoy referencing.
##############################################################################
"""
###Considerations:
#convert from yes/nolig to bound/apo
def main(args):
    """
    Main Function
    """
    init_args = init_arg_parsing(args.ligand_type, \
                    [args.cis_params_file, args.trans_params_file], \
                    args.ligand_params)
    
    if not (args.pickled_file or args.pickled_directory):
        #Starting up rosetta with the appropriate params files -
        print_out("emd182::Starting rosetta with the following parameters: " + init_args)
        pr.init(init_args)
        #sys.exit()

        #Setup data for constraints and symmdef files
        cstfiles, symdef = symm_cst_args_manage( args.enzdes_cstfiles, \
                                                args.symmdef_file)
        
        #Setting up score functions
        sf = pr.rosetta.core.scoring.ScoreFunction()
        sf.add_weights_from_file('ref2015')
        
        #Getting the full list of pdbs:
        if args.partial_set:
            pdb_list = glob('/'.join([args.source_directory, \
                                '*'+ args.partial_set+'*.pdb']))
            #wt_list = glob('/'.join([args.source_directory, '*min.pdb']))
            #pdb_list+=wt_list
        else:
            pdb_list = glob('/'.join([args.source_directory, '*.pdb']))
        
        pdb_list = list(set(pdb_list))
        print(pdb_list)
        
        wt_poses = wild_type_poses(args.wt_pdb, pdb_list, sf)
        print_out("Loaded and scored the wt_pose")
    
        #Residue selection for later analysis
        ligand = ligand_interpreter(args.ligand_type, \
                                    [args.residue_set, args.ligand_params])
        
        #Organizing the pdb files into 3 states:
        # Residue mutation, Unnatural, and apo/bound state
        print_out("Interpreting the list of pdbs in the directory: " + \
                    args.source_directory )
        mutant_nameframe = pdb_filename_data_interpreter(pdb_list)
        
        #Processing Data by unique groups - sets of pdbfilenames
        #Each row is a set of point mutant, uaa_type, and apo/bound state
        #Ending with the class for each collection of respective models
        mutant_dataframe = nameframe_pose_analyzer(mutant_nameframe, wt_poses,\
                        sf, args.source_directory, cst_files=cstfiles, \
                        symm_def_file=symdef)
        
        #As the mutant dataframe can't be saved with pyrosetta objects
        #They need to be collapsed into a more pickleable and perhaps
        #Faster to run DataFrame
        collapsed_dataframe = mutant_dataframe_collapse(mutant_dataframe)
        
        dumping_dataframe_to_pickle(collapsed_dataframe, args.out_file, \
                    out_part=args.partial_set, out_dir=args.out_directory)
        
        if args.analysis_only:
            print('Analyzing the data only, then dumping to a pickle file. \
                   Run again without the --analysis_only flag to output \
                   graphs and a full data analysis.')
            sys.exit()

    elif args.pickled_file or args.pickled_directory:
        collapsed_dataframe = pulling_dataframe_from_pickles( \
                pickled_file=args.pickled_file, pickled_directory=args.pickled_directory)
        print_out("Loaded data into collapsed_dataframe")
        print(collapsed_dataframe)
    
    else:
        print_out("You need to input either a dataset to analyze or a pickled \
                    File to analyze. If not, this script doesn't run.")
        sys.exit()
   
    """
    At this point, the DataFrame -- mutant_dataframe -- is a DataFrame containing
    single_mutant_classes of each residue, as well as access to other data..
    mutant_dataframe is organized as such: 
        - columns = ['residue','uaa','apo/bound', 'point_mutant_set']
    each row is unique for the combination of residue, uaa, apo/bound, to allow
    for various ways of sorting the data. 
        The main data containing object is each point_mutant_set, the class object
        that contains each respective pose, as well as methods to calculate 
        commonly requested information information, including averaging and 
        min/max values of the sets of the pdbs. 
    """
    ###Data analysis time - 
    #Defining commonly used labels
    pdb_name = ntpath.basename(args.wt_pdb.split(',')[0]).split('_')[0]
    apobound = ['nolig','yeslig']
    cis_trans = {'cis': args.cis_params_file.split('/')[-1].replace('.params',''), \
                 'trans': args.trans_params_file.split('/')[-1].replace('.params','')}
    cis_trans_uaa = [cis_trans[key] for key in cis_trans.keys()]
    
    graph_types = {'Total E':'delta_total_E',\
                   'Total+Cst E':['delta_total_E','delta_cstE']}
    if args.per_resi_graph:
        make_per_residue_graphs(collapsed_dataframe, args.out_directory, cst_add=False)
        make_per_residue_graphs(collapsed_dataframe, args.out_directory, cst_add=True)

    #Filtering on the values
    print("pre filtering")
    print(collapsed_dataframe)
    if args.filter:
        collapsed_dataframe = filter_dataframe(collapsed_dataframe, \
                    args.filter_on_E, args.filter_limits)
        print('post filtering')
        print(collapsed_dataframe)
    
    #First set of graphs are normalized against the wt.
    #This set is for bar graphs of normalized data.
    #columns = ['residue','uaa','apo/bound', 'point_mutant_set']
    for combo in itertools.product(*(apobound, cis_trans_uaa) ): 
        make_simple_bar_graphs(collapsed_dataframe, combo, cis_trans_uaa, \
                    args.out_directory, pdb_name)
    
    #Second set of graphs - cis minus trans states - one in apo state, one in bound
    for ab in apobound:
        make_x_minus_y_bar_graphs(collapsed_dataframe, [ab, cis_trans['cis']], \
                [ab, cis_trans['trans']], args.out_directory, pdb_name)
    
    #May want to add a conditional the check if constraint energies > 0
    min_max_avg = ['min']
    #Third set of graphs - Scatter Plots - bound vs apo delta scores
    graph_types = {'Total E':'delta_total_E',\
                   'Total+Cst E':['delta_total_E','delta_cstE']}
    for key, item in cis_trans.items():
        for gkey, gitem in graph_types.items(): #score variants
            for zoomish in [True, False]: #turn on and off zoom
                for ma in min_max_avg: #selecting which item type we have
                    curr_title = ' '.join(["Apo vs Bound", gkey, 'for', item, \
                        pdb_name+',','score='+ma, 'zoom='+str(zoomish)])
                    make_xla_vs_ylb_scatter_plots( collapsed_dataframe, \
                        ['nolig', item], ['yeslig', item], args.out_directory, \
                        xlabel='Apo REU', ylabel='Bound REU', mma=ma,\
                        zoom_graph=zoomish, title=curr_title, graph_type=gitem)

    #Fourth set of graphs - Scatter Plots for Cis minus Trans apo/bound variants
    d=10
    for gkey, gitem in graph_types.items(): #score variants
        for zoomish in [True, False]: #turn on and off zoom
            for ma in min_max_avg: #selecting which item type we have
                curr_title = ' '.join(["Cis vs Trans", gkey, 'for', 'AZOF', \
                    pdb_name+',','score-'+ma, 'zoom-'+str(zoomish)])
                make_xla_vs_ylb_scatter_plots( collapsed_dataframe,\
                    ['nolig', cis_trans['cis']], ['yeslig', cis_trans['cis']],\
                    args.out_directory, mma=ma, zoom_graph=zoomish, \
                    title=curr_title, graph_type=gitem,\
                    x_control_filt=['nolig', cis_trans['trans']],\
                    y_control_filt=['yeslig', cis_trans['trans']],\
                    xlabel='Cis'+' '*d+'Apo REU'+' '*d+'Trans',\
                    ylabel='Cis'+' '*d+'Bound REU'+' '*d+'Trans')
    #args.out_type = 'min,avg,max'
    #args.zoom = 'True,False'
    
if __name__ == '__main__':
    args = parse_args()
    main(args)

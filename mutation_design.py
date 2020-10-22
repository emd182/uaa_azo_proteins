# This script is for taking a selected residue and mutating it to a selection of residues around a ligand for mutation and design.
# Editor: Elliott Dolan, Khare Lab, Rutgers University, 2020
'''
   Usage:
         python step14_mutation.py <number after 'UM_' in the file name> <residue number to be mutated> <model number right before '.pdb' in the file name>
'''

#!/usr/bin/python

#Library importing
import pyrosetta as pr
import argparse
import sys
from os.path import isdir, join
from os import makedirs
import numpy as np
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
from pyrosetta.rosetta.protocols.enzdes import AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.simple_metrics.metrics import \
    TotalEnergyMetric, RMSDMetric, InteractionEnergyMetric

#Script starting from here
def parse_args():
    info = """
        This script should take a protein file and a few other parameters
        and output a series of files with a point mutation and relaxed proteins.
        """
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-p', '--input_pdb', type=str, required=True,
                        help='Input starting PDB file with included ligand to \
                        determine the residues around which to begin mutations\
                        / relaxes upon.')
    parser.add_argument('-num', '--residue_number', type=str, required=True,
                        help='Number of the residue for mutation to --unnatural\
                        denoted residue')
    parser.add_argument('-n', '--nstruct', type=int, required=False, default=5,
                        help='N-structs to run on the point mutant. Defaults to 5.')
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
    parser.add_argument('-r', '--radius', type=float, required=False, default=5.5,
                        help='Radius around the point mutant to minimize.')
    parser.add_argument('-uaa', '--unnatural', type=str, required=True,
                        help='The Unnatural amino acid params file for point \
                        mutants.')
    parser.add_argument('-rem', '--remove_ligand', action='store_true', \
                        default=False, help="Turn on to run without the ligand")
    parser.add_argument('-cstdeg', '--dihedral_constraint_degrees', type=str, \
                        default='0,180',required=False, \
                        help='Set degrees needed for the dihedral\
                        constraints - REQUIRED if more than 1 angle for dihedral')
    parser.add_argument('-rig_lig', '--rigid_ligand', action='store_true', \
                        default=False, help="Turn on to prevent ligand repacking")
    parser.add_argument('-sym', '--symmdef_file', type=str, default=False, \
                        help="Turn on symmetry with the added symmdef file")
    parser.add_argument('-d', '--design', default=False, \
                        help='Input to allow design at a specific radius from')
    parser.add_argument('-ditoms', '--dihedral_constraint_atoms', default=None, \
                        type=str, required=False, help='Include if there are atoms\
                        to constrain in the mutant, comma seperated by set of \
                        atoms and semi-colon separate the groups of dihedral atoms.')
    parser.add_argument('-enz', '--enzdes_constraint_file', default=None, \
                        type=str, required=False, help='Include if enzdes constraints\
                        need to be added to for the design process. Note that the \
                        constraints added here will not be scored in the final\
                        design selection') 
    parser.add_argument('-fr', '--fast_relax', default=False, action='store_true', \
                        help="Call to perform fast_relax instead of repack-min")
    parser.add_argument('-nocst', '--no_constraints', action='store_true', \
                        default=False, help='Input to turn off constraints')
    parser.add_argument('-o', '--out_suffix', type=str, required=False,
                        help='Name of the output suffix - will contain mutated \
                        position')
    parser.add_argument('-dir', '--out_directory', type=str, required=False,
                        help='Output Directory fo all output files.')

    return parser.parse_args()

def print_out(outstring):
    print("emd182::"+str(outstring))

def out_directory(directory):
    """ 
    Make an output directory if requested and doesn't exist 
    """
    print_out("Making the directory:"+str(directory))
    if directory:
        outdir = directory
        if not isdir(directory):
            makedirs(directory)
    else:
        outdir = ''
    return outdir

def coord_constrain_pose(pose):
   """
   Applies coordinate constraints if needed to selected residues
   """
   print_out("Applying Constraints")
   ccg = CoordinateConstraintGenerator()
   ccg.set_bounded_width(0.2)
   ccg.set_bounded(True)
   ccg.set_ca_only(True)
   ccg.set_sd(0.5)
   ac = AddConstraints()
   ac.add_generator(ccg)
   ac.apply(pose)

   return pose
####
#Rebuild for incorporation of design into the protocol - need to be able to
#include designable ranges - primarily from the intergroup selector
####

def task_factory_builder(repacking_residues='all', designing_residues=None, \
        extra_rotamers=[1,2]):
    #Setting this up for now, need to implement these conditions:
        #Repackable set, Designable Set, prevent ligand repacking
    #Building a Task Factory for repacking/design around the neighbors residues 
    #Task factory builder should accept: repackable residues, designable
    #residues, specifically preventing repacking residues, and prevent the rest.
    #tf = task_factory_builder(repack=repacking_neighbors) 

    #distance - distance for Neighborhood residue selector - need to remov
    print_out("Building task factory")
    #Start up Task factory, starting ex1 and ex2, and include current
    tf = TaskFactory()
    tf.push_back(IncludeCurrent())

    if extra_rotamers:
        for x in extra_rotamers:
            tf.push_back(ExtraRotamers(0, x, 1))    
    
    #Set up repack and prevent repack
    repack = RestrictToRepackingRLT()
    prevent = PreventRepackingRLT()
    prevent_everything_else = NotResidueSelector()

    residues_for_repack = None
    residues_for_design = None
    changing_residues = OrResidueSelector()

    if repacking_residues != 'all':
        print_out("Repacking these residues: " + str(repacking_residues))
        residues_for_repack = ResidueIndexSelector(','.join(repacking_residues))
        changing_residues.add_residue_selector(residues_for_repack)
        tf.push_back(OperateOnResidueSubset(repack, residues_for_repack))
        
    if designing_residues:
        residues_for_design = ResidueIndexSelector(designing_residues)
        changing_residues.add_residue_selector(residues_for_design)
    
    prevent_everything_else.set_residue_selector(changing_residues)

    if repack == 'all':
        tf.push_back(OperateOnResidueSubset(repack, prevent_everything_else))
    else:
        tf.push_back(OperateOnResidueSubset(prevent, prevent_everything_else))
    return tf

def build_move_map(chi=False, bb=False, jump=False):
    """
    Building movemap - by default, everything is Off
    """
    mm = pr.MoveMap()
    print_out("chi: " + str(chi))
    print_out("bb: " + str(bb))
    print_out("jump: " + str(jump))
    mm.set_chi(chi)
    mm.set_bb(bb)
    mm.set_jump(jump)
    return mm

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

def fast_relax_mutant(pose, task_factory, move_map, score_function, decoys=1):
    """
    Runs Fast Relax on the mutant pose instead of just repack an minimize.
    Can modify to take varying amounts of decoys.
    """
    fast_relax = FastRelax(1)
    fast_relax.set_scorefxn(score_function)
    fast_relax.set_task_factory(task_factory)
    fast_relax.set_movemap(move_map)
    score_function(pose)
    packer = task_factory.create_task_and_apply_taskoperations(pose)
    #print_out("packer task")
    #print_out(packer)
    traj_idx = 0 
    lowest_energy = 1000.0
    print_out("emd182::Running Fast Relax")
    while traj_idx < decoys:
        print_out("Round Number: " + str(traj_idx))
        pose_copy = pr.Pose()
        pose_copy.assign(pose)
        print_out("emd182:: pose total residues: " + str(pose.total_residue()) )
        print_out("emd182:: pose total residues: " + str(pose_copy.total_residue()) )
        fast_relax.apply(pose_copy)
        decoy_energy = total_energy(pose_copy, score_function)
        if traj_idx == '0':
            mutated_pose = pose_copy
            lowest_energy = decoy_energy
        elif decoy_energy < lowest_energy:
            mutated_pose = pose_copy
            lowest_energy = decoy_energy
        traj_idx += 1
    
    return mutated_pose

def repack_and_minimize_mutant(pose, task_factory, move_map, score_function,  rounds=3):
    #Copying the pose
    ram_pose = pr.Pose(pose)

    #preparing repack and applying
    prm = PackRotamersMover()
    prm.score_function(score_function)
    prm.task_factory(task_factory)
    
    #preparing minimize and applying
    min_mover = MinMover()
    min_mover.movemap(move_map)
    min_mover.score_function(score_function)
    
    print_out("Checking Packertask")
    packer = task_factory.create_task_and_apply_taskoperations(pose)
    print_out("packer task")
    print_out(packer)
    for rnd in range(rounds):
        print_out("round " + str(rnd+1) + " of repack and min")
        prm.apply(ram_pose)
        min_mover.apply(ram_pose)
    
    return ram_pose

def residue_selection(selection_one, selection_two):
    #For now, just an or selector to group up the residues.
    first = ResidueIndexSelector(selection_one)
    second = ResidueIndexSelector(selection_one)
    both = OrResidueSelector(first, second)
    return both

def ligand_neighbor_selection(lig, radius, pose, include_ligand=False):
    rad = float(radius)
    not_ligand = NotResidueSelector()
    not_ligand.set_residue_selector(lig)

    lig_distant_neighbors = InterGroupInterfaceByVectorSelector()
    lig_distant_neighbors.group1_selector(lig)
    lig_distant_neighbors.group2_selector(not_ligand)
    lig_distant_neighbors.cb_dist_cut(2.0*rad)
    lig_distant_neighbors.nearby_atom_cut(rad)

    lig_contacts = pr.rosetta.core.select.residue_selector.CloseContactResidueSelector()
    lig_contacts.central_residue_group_selector(lig)
    lig_contacts.threshold(rad)

    lig_neighbors = OrResidueSelector()
    lig_neighbors.add_residue_selector(lig_distant_neighbors)
    lig_neighbors.add_residue_selector(lig_contacts)

    if include_ligand:
        selection = AndResidueSelector()
        selection.add_residue_selector(not_ligand)
    else:
        selection = OrResidueSelector()
        selection.add_residue_selector(lig)

    selection.add_residue_selector(lig_neighbors)
    return get_residues_from_subset(selection.apply(pose))

#Simple definition file to convert a residue selector applied to a pose
#Into a vector of residue ints.
def selector_to_vector(residue_selector, pose):
    residue_vector = []
    resids = get_residues_from_subset(residue_selector.apply(pose))
    for res_int in resids:
        res = pose.pdb_info().pose2pdb(res_int).split()
        residue_vector.append(''.join(res))
    return residue_vector
    
#Adding match constraints to the pose
def apply_match_constraints(pose, match_cst_file):
    AoRMcsts = AddOrRemoveMatchCsts()
    AoRMcsts.set_cst_action(pr.rosetta.protocols.enzdes.CstAction.ADD_NEW)
    AoRMcsts.cstfile(match_cst_file)
    AoRMcsts.apply(pose)

def apply_dihedral_constraint(atom_str_list, residue_index, pose, \
            cst_degrees='0,180'):
    atoms = [] #pr.rosetta.utility.vector1_core_id_AtomID()
    index = get_residues_from_subset(residue_index.apply(pose))
    print(index)
    for atomstr in atom_str_list:
        new_atom = pr.rosetta.core.id.NamedAtomID(atomstr, int(index.back()))
        atoms.append(pr.rosetta.core.pose.named_atom_id_to_atom_id(new_atom, pose))
    ambiguous = pr.rosetta.core.scoring.constraints.AmbiguousConstraint()
    for degree_skip in [float(x) for x in cst_degrees.split(',')]: 
        dihedral_func = pr.rosetta.core.scoring.func.CircularHarmonicFunc( \
                np.radians(degree_skip), np.radians(3.0))
        dihedral_cst = pr.rosetta.core.scoring.constraints.DihedralConstraint(\
                atoms[0], atoms[1], atoms[2], atoms[3], dihedral_func)
        ambiguous.add_individual_constraint(dihedral_cst)
    pose.add_constraint(ambiguous)

#Repurpose into LINK check.
def check_residue_against_constraints(pdbfilename, residuenumber):
    f = open(args.input_pdb, 'r')
    idx = residuenumber[0:-1]
    for line in f.readlines():
        if 'REMARK' in line:
            items = [x for x in line.split(' ') if x != '']
            residue_set = [int(items[6]), int(items[11])]
            print(items)
            sys.exit()
            print(residue_set)
            print(idx)
            if int(idx) in residue_set:
                print('The entered residue is a part of the constraints. \
                        There will be no good coming from this.\
                        So the script is ending here..')
                sys.exit()

#Main function starts here.
###FUNCTION IDEAS###
def main(args):
    
    params = [args.unnatural]
    uaa = params[0].split('/')[-1].strip(".params")
    if args.ligand_type == 'ligand':
        params.append(args.ligand_name)
        lig_name = args.ligand_name.split('/')[-1].strip('.params')
    
    init_args = ' '.join(['-run:preserve_header','-extra_res_fa'] + params)
    
#Adding enzdes constraint file if necessary
    if args.enzdes_constraint_file:
        init_args = init_args + " -enzdes::cstfile " \
                    + str(args.enzdes_constraint_file)

#Starting up rosetta with the appropriate params files -
#-- need both unnatural aa and ligand file to be added (if the ligand isn't a protein).
    print_out("emd182::Starting rosetta with the following parameters: " + init_args)
    pr.init(init_args)
    file_options = pr.rosetta.core.io.StructFileRepOptions()
    file_options.set_preserve_header(bool(True))
    
    pose = pr.pose_from_pdb(args.input_pdb)
    
    sf = pr.rosetta.core.scoring.ScoreFunction()
    sf.add_weights_from_file('ref2015_cst')
    #Scoring the pose
    sf(pose)
    
    if args.symmdef_file:
        sfsm = SetupForSymmetryMover(args.symmdef_file)
        sfsm.apply(pose)
    #sys.exit()

    #Residue selection
    if args.residue_number != '0':
        delta_resi = ResidueIndexSelector(args.residue_number)
        #check_residue_against_constraints(args.input_pdb, args.residue_number)
    
    #Determining if the interacting ligand is a protein or small molecule
    #and selecting the appropriate residues
    if args.ligand_type == 'protein':
        ligand = ResidueIndexSelector(args.residue_set)
    elif args.ligand_type == 'ligand':
        ligand = ResidueNameSelector()
        ligand.set_residue_name3(lig_name)
    
    #Adding the unnatural at the mutation site
    print_out("Loading Unnatural " + uaa + \
        " onto residue number " + args.residue_number )
    
    if args.residue_number in selector_to_vector(ligand, pose):
        print_out("Selected residue number IS a part of the ligand. Don't do that. \
                Chosen residue number: " + str(args.residue_number) )
    #Setting up Mutation function on the correct residue, so long as the
    #residue isn't #0
    if args.residue_number != '0':
        mutater = pr.rosetta.protocols.simple_moves.MutateResidue()
        mutating_residue = ResidueIndexSelector(args.residue_number)
        mutater.set_selector(mutating_residue)
        mutater.set_res_name(uaa)
        print_out("emd182::loading residue to mutate into: " + str(args.residue_number))


    lig_vector = selector_to_vector(ligand, pose)
    if args.residue_number != '0':
        if args.ligand_type == 'protein' and delta_resi in lig_vector:
            print_out("Selected residue for mutation is part of input selection: " + \
                    delta_resi + " is selected, but is contained in " \
                    + str(lig_vector))
            print_out("Exiting the python script")
            sys.exit()

    #Build the MoveMap for the mutations.
    move_map = build_move_map(True, True, True)
    
    #Determining if the ligand should be removed or not
    out_ligand_file = 'yeslig'
    if args.remove_ligand:
        args.rigid_ligand = False
        out_ligand_file = 'nolig'

    print_out("emd182::Start mutations and relaxation script " \
        + str(args.nstruct) + " times.")
    
    if args.residue_number == '0':
        args.nstruct = 1

    for struct in range(0,args.nstruct):
        print_out(struct)
        mutant_pose = pr.Pose(pose)
        sf(mutant_pose)
        
        #Making appropriate residue selections base on if the ligand will be
        #Removed or not, as well as if the ligand is rigid or not.
        residues_around_ligand = ligand_neighbor_selection(ligand, args.radius, \
                mutant_pose, bool(args.rigid_ligand) )
        if args.residue_number != '0':
            residues_around_mutant = ligand_neighbor_selection(mutating_residue, \
                    args.radius, mutant_pose, True)

        design_around_mutant = None
        #if design is turned on, will design around the mutant.
        #Can add more variability later
        if args.design and args.residue_number != '0':
            designing_residues = ligand_neighbor_selection(mutating_residue, \
                    args.design, mutant_pose, False)
        
        #Combining the repacking neighborhoods around the ligand and mutant
        if args.residue_number != '0':
            repacking_neighbors = residue_selection(residues_around_ligand, \
                    residues_around_mutant)
        else:
            repacking_neighbors = ResidueIndexSelector(residues_around_ligand)

        #Specifically converting residue selectors to vectors, and removing the 
        #Appropriate residue ids from the lists of repacking or designing residues
        repacking_resids = selector_to_vector(repacking_neighbors, mutant_pose)
        if args.rigid_ligand:
            lig_resids = selector_to_vector(ligand, mutant_pose)
            for res in lig_resids:
                if res in repacking_resids:
                    repacking_resids.remove(res)
            if args.design:
                designing_resids = selector_to_vector(designing_residues, mutant_pose)
                for res in lig_resids:
                    if res in designing_resids:
                        repacking_resids.remove(res)
        
        #If the remove-ligand is called, will remove the ligand from the mutant pose
        if args.remove_ligand:
            removing_resids = selector_to_vector(ligand, mutant_pose)
            for res in removing_resids:
                if res in repacking_resids:
                    repacking_resids.remove(res)
            remove_lig = DeleteRegionMover()
            remove_lig.set_residue_selector(ligand)
            remove_lig.apply(mutant_pose)
            
        tf = task_factory_builder(repacking_residues=repacking_resids, \
                designing_residues=design_around_mutant)
        
        #apply mutation
        if args.residue_number != '0':
            mutater.apply(mutant_pose)
            if args.dihedral_constraint_atoms:
                dihedral_atoms = args.dihedral_constraint_atoms.split(';')
                dihedral_values = args.dihedral_constraint_degrees.split(';')
                for dihedrals in range(len(dihedral_atoms)):
                    apply_dihedral_constraint(dihedral_atoms[dihedrals].split(','), \
                                delta_resi, mutant_pose, dihedral_values[dihedrals])
       
        #turn on match constraints if needed:
        if args.enzdes_constraint_file:
            if not args.remove_ligand:
                print_out('The code will break if the constraint file is attached \
                        to the ligand that is being removed')
            apply_match_constraints(mutant_pose, args.enzdes_constraint_file)
            print('just checking')

        #turn off constraints if requested - default is on
        if not args.no_constraints:
            mutant_pose = coord_constrain_pose(mutant_pose)
        
        #Repack or minimize if selected
        if args.fast_relax:
            print_out("Running With Fast Relax")
            new_pose = fast_relax_mutant(mutant_pose, tf, move_map, sf)
        else:
            print_out("Running repack and minimize")
            new_pose = repack_and_minimize_mutant(mutant_pose, tf, move_map, sf)
        
        #output the name of the file
        #Includes original pdb namne, residue number and uaa, nstruct number,
        #and if the ligand is included or not.
        base_pdb_filename = args.input_pdb.split('/')[-1].strip(".pdb")
        outname = '{}_{}_{}_{}_{}.pdb'.format(base_pdb_filename, \
                args.residue_number, uaa, out_ligand_file, str(struct))
        if args.residue_number == '0':
            outname = '{}_{}_min.pdb'.format(base_pdb_filename, \
                    out_ligand_file)
        if args.out_suffix:
            outname = outname.strip(".pdb") + args.out_suffix + ".pdb"
        print_out("Outputting protein file as : " + outname )

        out_dir = out_directory(args.out_directory)
        outname = '/'.join([out_dir, outname])
        print_out("Writing Outputting protein file as : " + outname )
        
        new_pose.dump_pdb(outname)

if __name__ == '__main__':
    args = parse_args()
    main(args)

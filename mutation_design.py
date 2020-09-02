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
    parser.add_argument('-num', '--residue_number', type=int, required=True,
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
    parser.add_argument('-rig_lig', '--rigid_ligand', action='store_true', \
                        default=False, help="Turn on to prevent ligand repacking")
    parser.add_argument('-d', '--design', default=False, \
                        help='Input to allow design at a specific radius from')
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
   ccg.set_bounded_width(0.1)
   ccg.set_bounded(True)
   ccg.set_sidechain(False)
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
        residues_for_repack = ResidueIndexSelector(repacking_residues)
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
    print("chi: " + str(chi))
    print("bb: " + str(bb))
    print("jump: " + str(jump))
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

def fast_relax_mutant(pose, task_factory, move_map, score_function):
    """
    Runs Fast Relax on the mutant pose instead of just repack an minimize.
    Can modify to take varying amounts of decoys.
    """
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(score_function)
    fast_relax.set_task_factory(task_factory)
    fast_relax.set_movemap(move_map)
    score_function(pose)
    packer = task_factory.create_task_and_apply_task_operations(pose)
    print_out("packer task")
    print(packer)

    traj_idx = 0 
    decoys = 1
    print("emd182::Running Fast Relax")
    while traj_idx < decoys:
        print("Round Number: " + str(traj_idx))
        pose_copy = pr.Pose()
        pose_copy.assign(pose)
        print("emd182:: pose total residues: " + str(pose.total_residue()) )
        print("emd182:: pose total residues: " + str(pose_copy.total_residue()) )
        fast_relax.apply(pose_copy)
        decoy_energy = total_energy(pose_copy, score_function)
        if traj_idx == 0:
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
    
    print("Checking Packertask")
    packer = task_factory.create_task_and_apply_taskoperations(pose)
    print_out("packer task")
    print(packer)

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

    lig_neighbors = InterGroupInterfaceByVectorSelector()
    lig_neighbors.group1_selector(lig)
    lig_neighbors.group2_selector(not_ligand)
    lig_neighbors.cb_dist_cut(2.0*rad)
    lig_neighbors.nearby_atom_cut(rad)
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
    rosetta_vector = pr.rosetta.utility.vector1_unsigned_long()
    for res_int in [int(x) for x in \
                get_residues_from_subset(residue_selector.apply(pose)) ]:
        rosetta_vector.append(res_int)
    return rosetta_vector
    

#Main function starts here.
###FUNCTION IDEAS###
##add a count to determine the number of residues that change their rotameric
#state when they relax.
##measure changes in energy between original binding residues to the ligand
# compared to the remainding residues that still form interactions with the ligand.

def main(args):
    
    params = [args.unnatural]
    uaa = params[0].split('/')[-1].strip(".params")
    if args.ligand_type == 'ligand':
        params.append(args.ligand_name)
        lig_name = args.ligand_name.split('/')[-1].strip('.params')
    
    init_args = ' '.join(['-extra_res_fa'] + params)
    
#Starting up rosetta with the appropriate params files -
#-- need both unnatural aa and ligand file to be added (if the ligand isn't a protein).
    #parser.add_argument('-rem', '--remove_ligand', action='store_true', \
    #parser.add_argument('-rig_lig', '--rigid_ligand', action='store_true', \

    print("emd182::Starting rosetta with the following parameters: " + init_args)
    pr.init(init_args)
    pose = pr.pose_from_pdb(args.input_pdb)
    
    sf = pr.rosetta.core.scoring.ScoreFunction()
    sf.add_weights_from_file('ref2015_cst')

    #Residue selection
    delta_resi = ResidueIndexSelector(args.residue_number)
    if args.ligand_type == 'protein':
        ligand = ResidueIndexSelector(args.residue_set)

    elif args.ligand_type == 'ligand':
    #3-letter ligand name
        ligand = ResidueNameSelector()
        ligand.set_residue_name3(lig_name)
    
    print("emd182::Loading Unnatural " + uaa + \
        " onto residue number " + str( args.residue_number ) )

    #Setting up Mutation
    mutater = pr.rosetta.protocols.simple_moves.MutateResidue()
    mutating_residue = ResidueIndexSelector(int(args.residue_number))
    mutater.set_selector(mutating_residue)
    mutater.set_res_name(uaa)
    print("emd182::loading residue to mutate into: " + str(args.residue_number))

    #Build the MoveMap for the mutations.
    move_map = build_move_map(True, True, True)
    
    if args.remove_ligand:
        args.rigid_ligand = False

    print("emd182::Start mutations and relaxation script " \
        + str(args.nstruct) + " times.")

    for x in range(0,args.nstruct):
        #load / reload pose 
        mutant_pose = pr.Pose(pose)
        sf(mutant_pose)
        
        #Making appropriate residue selections base on if the ligand will be
        #Removed or not, as well as if the ligand is rigid or not.
        residues_around_ligand = ligand_neighbor_selection(ligand, args.radius, \
                mutant_pose, bool(args.rigid_ligand) )
        residues_around_mutant = ligand_neighbor_selection(mutating_residue, \
                args.radius, mutant_pose, True)

        design_around_mutant = None
        #if design is turned on, will design around the mutant.
        #Can add more variability later
        if args.design:
            designing_residues = ligand_neighbor_selection(mutating_residue, \
                    args.design, mutant_pose, False)
        #If the remove-ligand is called, will remove the ligand from the mutant pose
        if args.remove_ligand:
            removing_residues = [int(x) for x in get_residues_from_subset(ligand.apply(mutant_pose)) ]
            print_out("Removing the following residues, from " + \
                    str(removing_residues[0]) + " to " + str(removing_residues[-1]))
            pr.rosetta.protocols.grafting.delete_region(mutant_pose, \
                    removing_residues[0], removing_residues[-1])
        
        #Combining the repacking neighborhoods around the ligand and mutant
        repacking_neighbors = residue_selection(residues_around_ligand, \
                residues_around_mutant)

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

            
        tf = task_factory_builder(repacking_residues=repacking_resids, \
                designing_residues=design_around_mutant)
        
        #apply mutation
        mutater.apply(mutant_pose)
        
        #turn off constraints if requested - default is on
        if not args.no_constraints:
            pose = coord_constrain_pose(mutant_pose)
        
        #Repack or minimize if selected
        if args.design:
            print_out("Running With Design")
            new_pose = fast_relax_mutant(mutant_pose, tf, move_map, sf)
        else:
            print_out("Running repack and minimize")
            new_pose = repack_and_minimize_mutant(mutant_pose, tf, move_map, sf)
        
        outname = '{}_{}{}_{}.pdb'.format(args.input_pdb.strip('.pdb'), uaa, args.residue_number, str(x))
        if args.out_suffix:
            outname = outname.strip(".pdb") + args.out_suffix + ".pdb"
        print("emd182::Outputting protein file as : " + outname )

        out_dir = out_directory(args.out_directory)
        outname = join(out_dir, outname)

        print("emd182::Writing Outputting protein file as : " + outname )
        
        new_pose.dump_pdb(outname)

if __name__ == '__main__':
    args = parse_args()
    main(args)
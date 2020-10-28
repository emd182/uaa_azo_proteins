import pyrosetta as pr
from pyrosetta.rosetta.protocols.minimization_packing import \
    PackRotamersMover
from pyrosetta.rosetta.core.pack.rotamer_set import \
    RotamerSetFactory
import sys

pr.init('-extra_res_fa inputs/AZC.params -ex1 -ex2')
pose = pr.pose_from_pdb('TEV_solo.pdb')

print("emd182::Building score function")
sf = pr.rosetta.core.scoring.ScoreFunction()
sf.add_weights_from_file('ref2015')

print("emd182::Setting up and making a mutation")
res_mut = 30 
mutater = pr.rosetta.protocols.simple_moves.MutateResidue()
mutater.set_target(res_mut)
mutater.set_res_name('AZC')
mutater.apply(pose)

print("emd182::Making Packertask and restricting to repacking")
packer_task = pr.standard_packer_task(pose)
packer_task.restrict_to_repacking()
packer_task.set_bump_check(False)
pack_mover = PackRotamersMover(sf, packer_task)

rsf = RotamerSetFactory()
rs = rsf.create_rotamer_set(pose)
rs.set_resid(res_mut)
sf(pose)
packer_graph = pr.rosetta.core.pack.create_packer_graph(pose, sf, packer_task)
rs.build_rotamers(pose, sf, packer_task, packer_graph)

print(rs.rotamer(1).name(), rs.num_rotamers())
short_pose = pr.rosetta.core.pose.Pose()
short_pose.detached_copy(pose)
for i in range(1, rs.num_rotamers()+1):
    short_pose.residue(res_mut).set_all_chi(rs.rotamer(i).chi())
    short_pose.dump_pdb('mutated_pos_'+str(i)+'.pdb')#, res_mut_loc)

sys.exit()
#task_pack = pr.rosetta.core.pack.task.TaskFactory.create_packer_task(pose)
#print("Task pack is :"+ str(task_pack))
#sf.setup_for_packing(pose, task_pack.designing_residues(), task_pack.designing_residues())
#pack_graph = pr.rosetta.core.pack.create_packer_graph(pose, sf, task_pack) 
#
#rotset = pr.rosetta.core.pack.rotamer_set.RotamerSet()
#rotset.set_resid(int(30))
##pack_task = pr.rosetta.core.pack.task.PackerTask()
##pack_task.add_behavior
##rotset.build_rotamers(pose, sf, task_pack, pack_graph)
#rotset.show()
#rotset.initialize_pose_for_rotset_creation(pose)
##rotset.prepare_sets_for_packing(pose, sf)
#
##rotset.num_rotamers(20)
##Stuck here - pose, score_function, pack.task.Packertask, packer_neighbor_graph
#sys.exit()
#true_filt = pr.rosetta.protocols.filters.TrueFilter()
#explosion=1
#jump_num=0
#clash_check=False
#solo_res=True
#include_current=True
#
#tr = pr.rosetta.protocols.protein_interface_design.movers.TryRotamers('30', sf, true_filt, explosion, jump_num, clash_check, solo_res, include_current)
#print("Built tryrotamers object")
#tr.apply(pose)
#pose.dump_pdb('test.pdb')

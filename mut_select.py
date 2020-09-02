# This script is for selection of residues around a ligand for mutation and design.
# Editor: Elliott Dolan, Khare Lab, Rutgers University, 2020

'''
   Usage:
         python step14_mutation.py <number after 'UM_' in the file name> <residue number to be mutated> <model number right before '.pdb' in the file name>
'''

#!/usr/bin/python

#Library importing
import pyrosetta as pr
from pyrosetta.rosetta.core.select.residue_selector import \
    ResidueIndexSelector, ResidueNameSelector, \
    NeighborhoodResidueSelector, InterGroupInterfaceByVectorSelector, \
    NotResidueSelector, AndResidueSelector
from pyrosetta.rosetta.core.select import get_residues_from_subset
import argparse
import sys
from rosetta.protocols.simple_moves import MutateResidue

#Script starting from here
def parse_args():
    info = """
        This script should take a protein file and a few other parameters
        and output a single file with a list of residues available for mutation.
        """
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-p', '--input_pdb', type=str, required=True,
                        help='Input starting PDB file with included ligand to \
                        determine the residues around which to begin mutations \
                        / relaxes upon.')
    parser.add_argument('-l', '--ligand_type', type=str, required=True,
                        choices=['protein','ligand'],
                        help='Ligand type for identification - either 3-letter \
                        ligand name.')
    parser.add_argument('-res', '--residue_set', type=str, 
                        required=any(x in ['protein'] for x in sys.argv),
                        help='Comma delimited residue list for a peptide sequence \
                        behaving as a ligand in the mutation detection script.')
    parser.add_argument('-lig', '--ligand_name', type=str,  
                        required=any(x in ['ligand'] for x in sys.argv),
                        help='Ligand params file for ligand identification. Ligand\
                        params file should have a 3 letter name with the same 3 \
                        letters for the lig name.')
    parser.add_argument('-r', '--radius', type=float, required=False, default=8.0,
                        help='Radius around the bound ligand to identify residues \
                        to mutate.')
    parser.add_argument('-o', '--out_file', type=str, required=False,
                        help='Name of the output file with the list of mutatable\
                        positions')

    return parser.parse_args()

#Main function starts here.
def print_notice(outstring):
    print("emd182::"+str(outstring))

def main(args):
    
    if args.ligand_type == 'ligand':
        init_args = ' '.join(['-extra_res_fa', args.ligand_name])
        lig_name = args.ligand_name.split('/')[-1].strip('.params')
    else:
        init_args = ''
    
    pr.init(init_args)

    sf = pr.get_fa_scorefxn()
    sf.add_weights_from_file('ref2015')

    pose = pr.pose_from_pdb(args.input_pdb)
    sf(pose)

    print_notice("Scaffold protein Loaded Successfully!")
    print_notice("Scaffold protein has"+str(pose.total_residue())+"residues.")
    #Importing list of residues if the ligand is a protein
    if args.ligand_type == 'protein':
        ligres = ResidueIndexSelector(args.residue_set)
    #Targeting the ligand if the ligand isn't protein
    elif args.ligand_type == 'ligand':
        ligres = ResidueNameSelector()
        ligres.set_residue_name3(lig_name)

    print_notice("Ligand found at resnum: " + str(get_residues_from_subset( ligres.apply(pose) ) ) )
    
    #Setting the proteins not in the ligand
    not_ligand = NotResidueSelector()
    not_ligand.set_residue_selector(ligres)
    
    #Setting the protein considered part of the ligand
    ligand_interface = InterGroupInterfaceByVectorSelector()
    ligand_interface.group1_selector(ligres)
    ligand_interface.group2_selector(not_ligand)
    ligand_interface.cb_dist_cut(2.0*float(args.radius))
    ligand_interface.nearby_atom_cut(float(args.radius))

    #Proper Calculating - making sure only the non-ligand residues are selected
    non_lig_residues = AndResidueSelector()
    non_lig_residues.add_residue_selector(ligand_interface)
    non_lig_residues.add_residue_selector(not_ligand)
    
    #Collecting the residues from the subset
    neighbor_residues = get_residues_from_subset(non_lig_residues.apply(pose) )
    print_notice("Ligand found, neighbor residues are: " + ', '.join([str(x) for x in neighbor_residues]))

    if args.out_file:
        out_name = args.out_file
    else:
        out_name = args.input_pdb.strip('.pdb') + '_lig.pos'
    
    f = open(out_name, "w")
    for x in [str(x) for x in neighbor_residues]:
        f.write(x+'\n')
    f.close
    print("emd182::Wrote ligand position residues of pdb file " + args.input_pdb + " to filename: " + out_name)

if __name__ == '__main__':
    args = parse_args()
    main(args)

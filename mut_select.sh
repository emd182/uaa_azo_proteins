#!/bin/bash
#python mut_select.py -p TEV_orig.pdb -l protein -res 301C,302C,303C,304C,305C,306C,307C,308C,309C,310C -r 8.0 -o TEV.pos
uaa="/home/emd182/uaa_small_mol"
pdb="${uaa}/pdbs/TEV_fr.pdb"
pos="${uaa}/TEV.pos"
lig="${uaa}/inputs/${2}.params"

python mut_select.py -p ${pdb} -l protein -res 301C,302C,303C,304C,305C,306C,307C,308C,309C,310C -o ${pos} -r 6.0
#python mut_select.py -p ${pdb} -l ligand -lig ${lig} -o ${pos} -r 6.0

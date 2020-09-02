#!/bin/bash
#1 - protein for loading
#for TEV protease:
#-res 301C,302C,303C,304C,305C,306C,307C,308C,309C,310C
# -l protein
TEV="301C,302C,303C,304C,305C,306C,307C,308C,309C,310C"
LIG="inputs/DOG.params"
CURR="/home/emd182/uaa_small_mol/pdbs"
scripts=/home/emd182/uaa_small_mol/scripts
PDB="TEV_orig"
POS="${PDB}.pos"
N=1
#python ${scripts}/mut_select.py -p ${PDB}.pdb -l protein -res ${TEV} -r 9.0 -o ${POS}
for uaa in AZC AZT; 
do
    while read a;
    do
        UAA="inputs/${uaa}.params"
        python ${scripts}/mutation_design.py -p ${PDB}.pdb -l protein -res ${TEV} -num ${a} -r 6.0 -uaa ${UAA} -dir ${PDB}_AZO -n ${N} -rem
        #nohup python ${scripts}/mutation_design.py -p ${PDB}.pdb -l protein -res ${TEV} -num ${a} -r 6.0 -uaa ${UAA} -dir ${PDB}_AZO -n ${N} > logs/${PDB}_${a}${uaa}.log &
        sleep 1
        exit
        running=$( top -bn 1 | grep ' R ' | wc -l )
        while [ ${running} -ge 40 ]
        do
            echo "sleeping ... for 120 seconds"
            sleep 120
            running=$( top -bn 1 | grep ' R ' | wc -l )
        done
    done<${POS}
done

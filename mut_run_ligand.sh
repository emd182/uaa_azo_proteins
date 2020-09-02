#!/bin/bash
#1 - protein for loading
#for TEV protease:
#-res 301C,302C,303C,304C,305C,306C,307C,308C,309C,310C
# -l protein
TEV="301C,302C,303C,304C,305C,306C,307C,308C,309C,310C"
CURR="/home/emd182/uaa_small_mol"
LIG="${CURR}/inputs/${2}.params"
pdb="${CURR}/pdbs/${1}_fr.pdb"
N="1"
#python mut_select.py -p ${pdb} -l ligand -lig ${LIG} -r 10.0 -o ${1}.pos
for uaa in AZC AZT; 
do
    while read a;
    do
        UAA="${CURR}/inputs/${uaa}.params"
        #nohup python mutation_design.py -p ${pdb} -l ligand -lig ${LIG} -num ${a} -r 12 -uaa ${UAA} -cst -dir ${1}_${uaa}_${a} -n ${N} > ${1}.${uaa}.${a}.log &
        python mutation_design.py -p ${pdb} -l ligand -lig ${LIG} -num ${a} -r 6.0 -uaa ${UAA} -dir ${CURR}/${1}_${uaa}_nodesign -cst -n ${N} 
        sleep 1
        exit
        running=$( top -bn 1 | grep ' R ' | wc -l )
        while [ ${running} -ge 50 ]
        do
            echo "sleeping ... for 120 seconds"
            sleep 120
            running=$( top -bn 1 | grep ' R ' | wc -l )
        done
    done<${CURR}/${1}.pos
done

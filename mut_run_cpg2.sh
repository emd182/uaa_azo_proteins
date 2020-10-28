#!/bin/bash
dir="/home/emd182/uaa_small_mol"
scripts=${dir}/scripts
pdbloc="${dir}/pdbs"
inputs="${dir}/inputs"
p="1cg2"
f="${p}_INPUT_frB"
PDB="${f}.pdb"
POS="${f}.pos"
sym="${inputs}/${p}_fr.symm"
LIG="${inputs}/MTX.params"
enzdes="${inputs}/${p}_znmtx.cst"
N=1
rad=5.5

#Returned pos file is listed as pose numbering - not pdb numbering.
#python ${scripts}/mut_select.py -p ${pdbloc}/${PDB} -l ligand -lig ${LIG} -sym ${sym} -r ${rad} -o ${POS}
#exit
for aa in AZO; # APF; 
do
    for uaa in T C; 
    do
        ncaa=${aa:0:2}${uaa}
        UAA="inputs/${ncaa}.params"
        while read a;
        do
            #python ${scripts}/mutation_design.py -p ${pdbloc}/${PDB} -l ligand -lig ${LIG} -num ${a} -r ${rad} -uaa ${UAA} -dir ${p}_${aa}_fr -n ${N} -sym ${sym} -enz ${enzdes} -fr -rig_lig
            python ${scripts}/mutation_design.py -p ${pdbloc}/${PDB} -l ligand -lig ${LIG} -num ${a} -r ${rad} -uaa ${UAA} -dir ${p}_${aa}_fr -n ${N} -rem -enz ${enzdes} -ditoms CE1,CZ,N2,N3 -cstdeg 0.0,180.0 -sym ${sym}
            exit
            nohup python ${scripts}/mutation_design.py -p ${pdbloc}/${PDB} -l ligand -lig ${LIG} -num ${a} -r ${rad} -uaa ${UAA} -dir ${f}_${aa}_fr -n ${N} -rig_lig -fr &
            nohup python ${scripts}/mutation_design.py -p ${pdbloc}/${PDB} -l ligand -lig ${LIG} -num ${a} -r ${rad} -uaa ${UAA} -dir ${f}_${aa}_fr -n ${N} -rem -fr &
            sleep 1
            running=$( top -bn 1 | grep ' R ' | wc -l )
            while [ ${running} -ge 50 ]
            do
                echo "sleeping ... for 10 seconds"
                sleep 20
                running=$( top -bn 1 | grep ' R ' | wc -l )
            done
        done<${POS}
    done
done


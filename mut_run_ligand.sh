#!/bin/bash
dir="/home/emd182/uaa_small_mol"
scripts=${dir}/scripts
pdbloc="${dir}/pdbs"
while (( $# > 1 ))
do
    case $1 in 
        -tp) template="$2";;
        -lig) ligname="$2";;
        *) break ;
    esac; shift 2
done

if ! [ -z "${template}" ]
    then
    template=${template}
fi
if ! [ -z "${ligname}" ]
    then
    ligname=${ligname}
fi

TEV="301C,302C,303C,304C,305C,306C,307C,308C,309C,310C"
LIG="${dir}/inputs/${ligname}.params"
f="${template}"
PDB="${f}.pdb"
POS="${f}.pos"
N=5
rad=5.0
#Returned pos file is listed as pose numbering - not pdb numbering.
python ${scripts}/mut_select.py -p ${pdbloc}/${PDB} -l ligand -lig ${LIG} -r ${rad} -o ${POS}
#exit
for aa in AZO APF;
do
    for uaa in T C; 
    do
        ncaa=${aa:0:2}${uaa}
        echo $ncaa
        #exit
        UAA="inputs/${ncaa}.params"
        while read a;
        do
            nohup python ${scripts}/mutation_design.py -p ${pdbloc}/${PDB} -l ligand -lig ${LIG} -num ${a} -r ${rad} -uaa ${UAA} -dir ${f}_${aa}_packmin -n ${N} -rig_lig &
            nohup python ${scripts}/mutation_design.py -p ${pdbloc}/${PDB} -l ligand -lig ${LIG} -num ${a} -r ${rad} -uaa ${UAA} -dir ${f}_${aa}_packmin -n ${N} -rem &
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


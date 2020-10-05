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
inp=${dir}/inputs
LIG="${inp}/${ligname}.params"
f="${template}"
PDB="${f}.pdb"
POS="${f}.pos"
N=5
proc_count=50
rad=5.0

#Returned pos file is listed as pose numbering - not pdb numbering.
#python ${scripts}/mut_select.py -p ${pdbloc}/${PDB} -l protein -res ${TEV} -r ${rad} -o ${POS}
#exit
for uaa in AZC AZT; # APF; 
do
    while read a;
    do
        #uaa=APT
        #a=148A
        #N=1
        UAA="inputs/${uaa}.params"
        python ${scripts}/mutation_design.py -p ${pdbloc}/${PDB} -l protein -res ${TEV} -num ${a} -r 6.0 -uaa ${UAA} -dir ${f}_AZO_fr -n 1 -fr -enz ${inp}/tev.cst -ditoms CE1,CZ,N2,N3 
        #python ${scripts}/mutation_design.py -p ${pdbloc}/${PDB} -l protein -res ${TEV} -num ${a} -r 6.0 -uaa ${UAA} -dir ${f}_AZO_fr -n 1 -fr -enz ${inp}/tev.cst -ditoms CE1,CZ,N2,N1 
        exit
        nohup python ${scripts}/mutation_design.py -p ${pdbloc}/${PDB} -l protein -res ${TEV} -num ${a} -r ${rad} -uaa ${UAA} -dir ${f}_AZO_fr -n ${N} -fr -enz ${inp}/tev.cst -ditoms CE1,CZ,N2,N3 &
        nohup python ${scripts}/mutation_design.py -p ${pdbloc}/${PDB} -l protein -res ${TEV} -num ${a} -r ${rad} -uaa ${UAA} -dir ${f}_AZO_fr -n ${N} -rem -fr -ditoms CE1,CZ,N2,N3 &
        sleep 1
        #exit
        running=$( top -bn 1 | grep ' R ' | wc -l )
        while [ ${running} -ge ${proc_count} ]
        do
            echo "sleeping ... for 10 seconds"
            sleep 10
            running=$( top -bn 1 | grep ' R ' | wc -l )
        done
    done<${POS}
done


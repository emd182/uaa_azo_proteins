#!/bin/bash
dir="/home/emd182/uaa_small_mol"
scripts=${dir}/scripts
#for TEV protease:
#-res 301C,302C,303C,304C,305C,306C,307C,308C,309C,310C
# -l protein
while (( $# > 1 ))
do
    case $1 in 
        -tp) template="$2";;
        -pos) proligdetect="$2";;
        *) break ;
    esac; shift 2
done

#if ! [ -z "${template}" ]
#    then
#    template=
#fi

TEV="301C,302C,303C,304C,305C,306C,307C,308C,309C,310C"
LIG="inputs/DOG.params"
PDB="TEVorig"
OUT="${PDB}_AZO_fr"
WT1="${OUT}/${PDB}_nolig_min.pdb"
WT2="${OUT}/${PDB}_yeslig_min.pdb"
INP="${dir}/inputs"
cis="${INP}/AZC.params"
trans="${INP}/AZT.params"

python ${scripts}/mut_analyze.py -p ${WT1},${WT2} -l protein -res ${TEV} \
        -uaa -cis ${cis} -trans ${trans} -pre ${OUT}.pickle \
        -dir ${OUT} -od ${OUT}_analysis/ 

#python ${scripts}/mut_analyze.py -p ${WT1},${WT2} -l protein -res ${TEV} \
#        -uaa -cis ${cis} -trans ${trans} -pre ${OUT}.pickle \
#        -dir ${OUT} -od ${OUT}_analysis/ 

#nohup python ${scripts}/mut_analyze.py -p ${WT1},${WT2} -l protein -res ${TEV} \
#        -uaa -cis ${cis} -trans ${trans} -out ${PDB}_AZOpm.pickle \
#        -dir ${PDB}_AZO_packmin -od ${PDB}_AZOpm_analysis/ &
#python ${scripts}/mut_analyze.py -p ${WT1},${WT2} -l protein \
#        -res ${TEV} -uaa -cis ${cis} -trans ${trans} \
#        -pre ${PDB}.pickle -dir ${PDB}_AZO -od ${PDB}_AZO_analysis/


#!/bin/bash
dir="/home/emd182/uaa_small_mol"
scripts=${dir}/scripts
#for TEV protease:
#-res 301C,302C,303C,304C,305C,306C,307C,308C,309C,310C
# -l protein

TEV="301C,302C,303C,304C,305C,306C,307C,308C,309C,310C"
LIG="inputs/DOG.params"
PDB="TEVorig"
OUT="${PDB}_AZO_fr"
WT1="${OUT}/${PDB}_nolig_min.pdb"
WT2="${OUT}/${PDB}_yeslig_min.pdb"
INP="${dir}/inputs"
cis="${INP}/AZC.params"
trans="${INP}/AZT.params"

#sed -i '/LINK/d' ${OUT}/*yeslig*pdb
echo "starting python script:"

python ${scripts}/mut_analyze.py -p ${WT1},${WT2} -l protein -res ${TEV} \
        -uaa -cis ${cis} -trans ${trans} -pre ${OUT}_all.pickle \
        -dir ${OUT} -od ${OUT}_nolimit_analysis/ -lim 150

#python ${scripts}/mut_analyze.py -p ${WT1},${WT2} -l protein -res ${TEV} \
#        -uaa -cis ${cis} -trans ${trans} -pre ${OUT}.pickle \
#        -dir ${OUT} -od ${OUT}_analysis/ 

#nohup python ${scripts}/mut_analyze.py -p ${WT1},${WT2} -l protein -res ${TEV} \
#        -uaa -cis ${cis} -trans ${trans} -out ${PDB}_AZOpm.pickle \
#        -dir ${PDB}_AZO_packmin -od ${PDB}_AZOpm_analysis/ &
#python ${scripts}/mut_analyze.py -p ${WT1},${WT2} -l protein \
#        -res ${TEV} -uaa -cis ${cis} -trans ${trans} \
#        -pre ${PDB}.pickle -dir ${PDB}_AZO -od ${PDB}_AZO_analysis/


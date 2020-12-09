#!/bin/bash
dir="/home/emd182/uaa_small_mol"
scripts=${dir}/scripts
#for TEV protease:
#-res 301C,302C,303C,304C,305C,306C,307C,308C,309C,310C
# -l protein

TEV="301B,302B,303B,304B,305B,306B,307B,308B,309B,310B"
PDB="TEVorig"
u="AZ"
ua="${u}O"
OUT="${PDB}_${ua}_fr"
WT1="${OUT}/${PDB}_nolig_min.pdb"
WT2="${OUT}/${PDB}_yeslig_min.pdb"
INP="${dir}/inputs"
cst="${INP}/tev.cst,${INP}/tev_noligcst.cst"
cis="${INP}/${u}C.params"
trans="${INP}/${u}T.params"
pickle="${OUT}_${ua}.pickle"

protein_info="-p ${WT1},${WT2} -l protein -res ${TEV}"
uaa_info="-uaa -cis ${cis} -trans ${trans} -enz ${cst}"
#data_info="-out ${pickle} -dir ${OUT} -od ${OUT}_analysis"
data_info="-pre ${pickle} -dir ${OUT} -od ${OUT}_analysis"
data_filter="-filt -lim 15 -fval delta_total_E,nolig"

#sed -i '/LINK/d' ${OUT}/*pdb
echo "starting python script:"

python ${scripts}/mutant_analysis.py ${protein_info} ${uaa_info} ${data_info} ${data_filter}

#nohup python ${scripts}/mut_analyze.py -p ${WT1},${WT2} -l protein -res ${TEV} \
#        -uaa -cis ${cis} -trans ${trans} -out ${PDB}_AZOpm.pickle \
#        -dir ${PDB}_AZO_packmin -od ${PDB}_AZOpm_analysis/ &
#python ${scripts}/mut_analyze.py -p ${WT1},${WT2} -l protein \
#        -res ${TEV} -uaa -cis ${cis} -trans ${trans} \
#        -pre ${PDB}.pickle -dir ${PDB}_AZO -od ${PDB}_AZO_analysis/


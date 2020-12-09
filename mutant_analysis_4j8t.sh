#!/bin/bash
dir="/home/emd182/uaa_small_mol"
scripts=${dir}/scripts
#for TEV protease:
#-res 301C,302C,303C,304C,305C,306C,307C,308C,309C,310C
# -l protein

#TEV="301B,302B,303B,304B,305B,306B,307B,308B,309B,310B"
LIG="${dir}/inputs/DOG.params"
PDB="4j8t"
u="AZ"
ua="${u}O"
OUT="${PDB}_${ua}"
WT1="${OUT}/${PDB}_nolig_min.pdb"
WT2="${OUT}/${PDB}_yeslig_min.pdb"
INP="${dir}/inputs"
#cst="${INP}/1cg2_zn.cst,${INP}/1cg2_zn_noligcst.cst"
cis="${INP}/${u}C.params"
trans="${INP}/${u}T.params"
pickle="${OUT}_${ua}.pickle"

#protein_info="-p ${WT1},${WT2} -l protein -res ${TEV}"
#data_info="-out ${pickle} -dir ${OUT} -od ${OUT}_analysis"

protein_info="-p ${WT1},${WT2} -l ligand -lig ${LIG}"
uaa_info="-uaa -cis ${cis} -trans ${trans}"
data_info="-dir ${OUT} -od ${OUT}_analysis"
#data_filter="-filt -lim 15 -fval delta_total_E,nolig"

#sed -i '/LINK/d' ${OUT}/*pdb
echo "starting python script:"

#while read pos;
#do
#  for ab in nolig yeslig;
#  do 
#    for ct in ${u}C ${u}T;
#    do
#      if [[ "${pos}" != "0" ]]; then
#        outinfo="--out_file ${pickle} --partial_set _${pos}_${ct}_${ab}"  
#        nohup python ${scripts}/mutant_analysis.py ${protein_info} ${uaa_info} ${data_info} ${outinfo} -solo &
#        sleep 1
#        running=$( top -bn 1 | grep ' R ' | wc -l )
#        while [ ${running} -ge 40 ]
#        do
#            echo "sleeping while analysis runs"
#            sleep 30
#            running=$( top -bn 1 | grep ' R ' | wc -l )
#        done
#      fi
#    done
#  done
#done<${PDB}.pos

#exit

outinfo="--pickled_directory ${OUT}_analysis" 
python ${scripts}/mutant_analysis.py ${protein_info} ${uaa_info} ${data_info} ${outinfo} --per_resi_graph

#nohup python ${scripts}/mut_analyze.py -p ${WT1},${WT2} -l protein -res ${TEV} \
#        -uaa -cis ${cis} -trans ${trans} -out ${PDB}_AZOpm.pickle \
#        -dir ${PDB}_AZO_packmin -od ${PDB}_AZOpm_analysis/ &
#python ${scripts}/mut_analyze.py -p ${WT1},${WT2} -l protein \
#        -res ${TEV} -uaa -cis ${cis} -trans ${trans} \
#        -pre ${PDB}.pickle -dir ${PDB}_AZO -od ${PDB}_AZO_analysis/


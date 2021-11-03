#!/bin/bash
###Directories
dir="/home/emd182/uaa_small_mol"
scripts="${dir}/scripts"
INP="${dir}/inputs"

###PDB Information
PDB="1cg2"
ex="_closed"
prep="${dir}/prepared/${PDB}${ex}"
cst="${prep}/${PDB}${ex}_mtz.cst,${prep}/${PDB}${ex}_mtz_noligcst.cst"
sym="${prep}/${PDB}${ex}_INPUT.symm"
POS="${prep}/${PDB}${ex}_INPUT.pos"

###Ligand information
LIG="${prep}/MTZ.params"

###UAA Information
u="AZ"
ua="${u}O"
cis="${INP}/${u}C.params"
trans="${INP}/${u}T.params"

###Wild type file locations
out="${PDB}${ex}_out/${ua}_fr"
WT1="${out}/${PDB}_nolig_min.pdb"
WT2="${out}/${PDB}_yeslig_min.pdb"

###Pickle Information
pickle="${PDB}${ex}_analysis"

###Cumulative info
protein_info="-p ${WT1},${WT2} -l ligand -lig ${LIG} -sym ${sym} -enz ${cst}"
uaa_info="-uaa ${ua} -cis ${cis} -trans ${trans}"
data_info="-dir ${out} -od ${out}_analysis"
data_filter="-filt -lim 15 -fval delta_total_E,nolig"

#sed -i '/LINK/d' ${out}/*pdb
echo "starting python script:"

#while read pos;
#do
#  for ab in nolig yeslig; ##Selects yeslig or nolig as 
#  do 
#    for ct in ${u}C ${u}T; ##Selecting the appropriate unnatural
#    do
#      if [[ "${pos}" != "0" ]]; then
#        outinfo="--out_file ${pickle} --partial_set ${pos}_${ct}_${ab}"  
#        nohup python ${scripts}/mutant_analysis.py ${protein_info} ${uaa_info} ${data_info} ${outinfo} -solo &
#        sleep 0.125
#        running=$( top -bn 1 | grep ' R ' | wc -l )
#        while [ ${running} -ge 60 ]
#        do
#            echo "sleeping while analysis runs"
#            sleep 10
#            running=$( top -bn 1 | grep ' R ' | wc -l )
#        done
#      fi
#    done
#  done
#done<${POS}
#
#exit

analysis="${out}_analysis"
outinfo="--pickled_directory ${analysis}" 
data_filter="-filt -lim 15 -fval delta_total_E,nolig"
python ${scripts}/mutant_analysis.py ${protein_info} ${uaa_info} ${data_info} ${outinfo} ${data_filter}
mkdir ${analysis}/filtered
mv ${analysis}/*png ${analysis}/filtered/
python ${scripts}/mutant_analysis.py ${protein_info} ${uaa_info} ${data_info} ${outinfo} # --per_resi_graph
mkdir ${analysis}/unfiltered
mv ${analysis}/*png ${analysis}/unfiltered/


#!/bin/bash
###Directories
dir="/home/emd182/uaa_small_mol"
scripts="${dir}/scripts"
inputs="${dir}/inputs"
params="${inputs}/params_files"

###PDB Information
pdb="1cg2"
pdb_loc="${inputs}/pdb_files/${pdb}"
sym_loc="${pdb_loc}_INPUT.symm"
pos_loc="${pdb_loc}_INPUT.pos"
cst_loc="${inputs}/cst_files/${pdb}_mtz"
cst="${cst_loc}.cst,${cst_loc}_noligcst.cst"

###Ligand information
lig_loc="${params}/MTZ.params"

###UAA Information
ua="FC"
uaa="F6AzoC"
cis="${params}/${ua}C.params"
trans="${params}/${ua}T.params"

###Wild type file locations
out="${pdb}_out/${uaa}_fr"
wt1="${out}/${pdb}_nolig_min.pdb"
wt2="${out}/${pdb}_yeslig_min.pdb"

###Pickle Information
pickle="${PDB}_analysis"

###Cumulative info
protein_info="-p ${wt1},${wt2} -l ligand -lig ${lig_loc} -sym ${sym_loc} -enz ${cst}"
uaa_info="-uaa ${uaa} -cis ${cis} -trans ${trans}"
data_info="-dir ${out} -od ${out}_analysis"
data_filter="-filt -lim 15 -fval delta_total_E,nolig"

#sed -i '/LINK/d' ${out}/*pdb
#exit
echo "starting python script:"

#while read num;
#do
#  for ab in nolig yeslig; ##Selects yeslig or nolig as 
#  do 
#    for ct in ${ua}C ${ua}T; ##Selecting the appropriate unnatural
#    do
#      if [[ "${num}" != "0" ]]; then
#        outinfo="--out_file ${pickle} --partial_set ${num}_${ct}_${ab}"  
#        nohup python ${scripts}/mutant_analysis.py ${protein_info} ${uaa_info} ${data_info} ${outinfo} -solo &
#        #python ${scripts}/mutant_analysis.py ${protein_info} ${uaa_info} ${data_info} ${outinfo} -solo
#        #exit
#        sleep 0.05
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
#done<${pos_loc}
#
#exit

analysis="${out}_analysis"
outinfo="--pickled_directory ${analysis}" 
data_filter="-filt -lim 15 -fval delta_total_E,nolig"
#python ${scripts}/mutant_analysis.py ${protein_info} ${uaa_info} ${data_info} ${outinfo} ${data_filter}
#mkdir ${analysis}/filtered
#mv ${analysis}/*png ${analysis}/filtered/
#python ${scripts}/mutant_analysis.py ${protein_info} ${uaa_info} ${data_info} ${outinfo}
#mkdir ${analysis}/unfiltered
#mv ${analysis}/*png ${analysis}/unfiltered/
python ${scripts}/mutant_analysis.py ${protein_info} ${uaa_info} ${data_info} ${outinfo} --per_resi_graph min
mkdir ${analysis}/unfiltered
mv ${analysis}/*png ${analysis}/unfiltered/


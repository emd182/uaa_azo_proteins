#!/bin/bash
dir="/home/emd182/uaa_small_mol"
scripts=${dir}/scripts

LIG="inputs/DOG.params"
INP="${dir}/inputs"
PDB="${1}"
#UAA=${2}
#DIR=${PDB}_${UAA}_fr
#WT1="${DIR}/${PDB}_nolig_min.pdb"
#WT2="${DIR}/${PDB}_yeslig_min.pdb"
#INP="${dir}/inputs"
#part=${UAA:0:2}
#
#cis="${INP}/${part}C.params"
#trans="${INP}/${part}T.params"
for pdb in 4j8tfr 4j9afr;
do
    PDB="${pdb}"
    for UAA in AZO APF; 
    do
        for relax in fr packmin;
        do
            DIR=${PDB}_${UAA}_${relax}
            WT1="${DIR}/${PDB}_nolig_min.pdb"
            WT2="${DIR}/${PDB}_yeslig_min.pdb"
            part=${UAA:0:2}
            cis="${INP}/${part}C.params"
            trans="${INP}/${part}T.params"
            analyze="${PDB}_${UAA}_${relax}"
            python ${scripts}/mut_analyze.py -p ${WT1},${WT2} \
                    -l ligand -lig ${LIG} \
                    -uaa -cis ${cis} -trans ${trans} \
                    -pre ${analyze}.pickle \
                    -dir ${analyze} \
                    -od ${analyze}_analysis/  
            #exit
        
        done
    done
done
#python ${scripts}/mut_analyze.py -p ${WT1},${WT2} -l ligand \
#        -lig ${LIG} -uaa -cis ${cis} -trans ${trans} \
#        -out ${PDB}_all.pickle -dir ${PDB}_AZO -od ${PDB}_AZO_analysis/


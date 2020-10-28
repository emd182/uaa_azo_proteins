#!/bin/bash
dir="/home/emd182/uaa_small_mol"
scripts=${dir}/scripts
pdbloc="${dir}/pdbs"
inp=${dir}/inputs

f="TEVorig"
PDB="${f}.pdb"
POS="${f}.pos"

N=1
proc_count=55
rad=5.5
outdir="${f}_APF_fr"

uat=APT
uac=APC
UAT="inputs/${uat}.params"
UAC="inputs/${uac}.params"
#LIG="${inp}/${ligname}.params"
TEV="301B,302B,303B,304B,305B,306B,307B,308B,309B,310B"

diatomAZO="CE1,CZ,N2,N3"
diatomAPF="CE1,CZ,N1,N2;N1,N2,C3,C4"
diatom=${diatomAPF}

degAZT="0.0,180.0"
degAZC="49.7,-130.4"
degAPC="50.7,-129.3;25.6,-154.4"
degAPT="0.0,180.0;0.0,180.0"

protein_info=" -p ${pdbloc}/${PDB} -l protein -res ${TEV} -enz ${inp}/tev.cst"
run_info="-r ${rad} -dir ${outdir} -n ${N} -fr"


#Returned pos file is listed as pose numbering - not pdb numbering.
python ${scripts}/mut_select.py -p ${pdbloc}/${PDB} -l protein -res ${TEV} -r ${rad} -o ${POS}
while read a;
do
    #a=178A
    #N=1
    #uaa=AZC
    #UAA="inputs/${uaa}.params"
    #python ${scripts}/mutation_design.py -p ${pdbloc}/${PDB} -l protein -res ${TEV} -num ${a} -r 6.0 -uaa ${UAA} -dir ${f}_AZO_fr -n 1 -fr -enz ${inp}/tev.cst -ditoms CE1,CZ,N2,N3 -cstdeg 49.7,-130.4 
    ###Unnatural Trans runs
    nohup python ${scripts}/mutation_design.py ${protein_info} ${run_info} \
            -num ${a} -uaa ${UAT} -ditoms ${diatom} -cstdeg ${degAPT} &
    nohup python ${scripts}/mutation_design.py ${protein_info} ${run_info} \
            -num ${a} -uaa ${UAT} -ditoms ${diatom} -cstdeg ${degAPT} -rem &
    
    ###Unnatural Cis Runs
    nohup python ${scripts}/mutation_design.py ${protein_info} ${run_info} \
            -num ${a} -uaa ${UAC} -ditoms ${diatom} -cstdeg ${degAPC} &
    nohup python ${scripts}/mutation_design.py ${protein_info} ${run_info} \
            -num ${a} -uaa ${UAC} -ditoms ${diatom} -cstdeg ${degAPC} -rem &
    
    #exit
    sleep 1
    running=$( top -bn 1 | grep ' R ' | wc -l )
    while [ ${running} -ge ${proc_count} ]
    do
        echo "sleeping ... for 10 seconds"
        sleep 10
        running=$( top -bn 1 | grep ' R ' | wc -l )
    done
done<${POS}


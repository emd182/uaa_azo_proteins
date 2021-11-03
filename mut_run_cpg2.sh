#!/bin/bash
dir="/home/emd182/uaa_small_mol"
scripts="${dir}/scripts"
inputs="${dir}/inputs"

p="1cg2"
prep="${dir}/prepared/${p}"
f="${p}_INPUT"
PDB="${f}.pdb"
POS="${prep}/${f}.pos"

sym="${prep}/${f}.symm"
LIG="${prep}/MTZ.params"
enzdes="${prep}/${p}_mtz.cst"

N=25
proc_count=50
rad=5.5

aa="AZ"
uat=${aa}T
uac=${aa}C
UAT="inputs/${uat}.params"
UAC="inputs/${uac}.params"
outdir="${dir}/${p}_out/AZO_fr"

diatomAZO="CE1,CZ,N2,N3"
diatomAPF="CE1,CZ,N1,N2;N1,N2,C3,C4"
diatom=${diatomAZO}

degAZT="0.0,180.0"
degAZC="49.7,-130.4"
degAPC="50.7,-129.3;25.6,-154.4"
degAPT="0.0,180.0;0.0,180.0"

degT=${degAZT}
degC=${degAZC}

pdbinfo="-p ${prep}/${PDB} -l ligand -lig ${LIG} -sym ${sym}"
runinfo="-r ${rad} -dir ${outdir} -n ${N} -fr"
cstinfo="-enz ${enzdes} -ditoms ${diatom}"
#Other options : -rig_lig, -rem, -fr

python ${scripts}/mut_select.py -p ${prep}/${PDB} -l ligand -lig ${LIG} -r ${rad} -o ${POS}
#exit
while read a; 
do

    #python ${scripts}/mutation_design.py ${pdbinfo} ${runinfo} ${cstinfo} \
    #        -cstdeg ${degAZT} -num ${a} -uaa ${UAT} 
    #python ${scripts}/mutation_design.py ${pdbinfo} ${runinfo} ${cstinfo} \
    #        -cstdeg ${degAZT} -num ${a} -uaa ${UAT} -rem B
    #exit
    ###Unnatural Trans runs
    nohup python ${scripts}/mutation_design.py ${pdbinfo} ${runinfo} ${cstinfo} \
            -cstdeg ${degAZT} -num ${a} -uaa ${UAT} &
    nohup python ${scripts}/mutation_design.py ${pdbinfo} ${runinfo} ${cstinfo} \
            -cstdeg ${degAZT} -num ${a} -uaa ${UAT} -rem B &
    
    ###Unnatural Cis Runs
    nohup python ${scripts}/mutation_design.py ${pdbinfo} ${runinfo} ${cstinfo} \
            -cstdeg ${degAZC} -num ${a} -uaa ${UAC} &
    nohup python ${scripts}/mutation_design.py ${pdbinfo} ${runinfo} ${cstinfo} \
            -cstdeg ${degAZC} -num ${a} -uaa ${UAC} -rem B &
    
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


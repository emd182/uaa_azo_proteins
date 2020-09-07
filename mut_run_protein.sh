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

TEV="301C,302C,303C,304C,305C,306C,307C,308C,309C,310C"
LIG="inputs/DOG.params"
PDB="TEVorig"
POS="${PDB}.pos"
txt="${PDB}".postxt
N=20
#nohup python ${scripts}/mut_select.py -p ${PDB}.pdb -l protein -res ${TEV} -r 9.0 -o ${POS} > out.check &
python ${scripts}/mut_select.py -p ${PDB}.pdb -l protein -res ${TEV} -r 12.0 -o ${POS} > ${txt}
if ! [ -z "${proligdetect}" ]
    then
        rm ${POS}
        frompos=$( grep "protocols.enzdes.Enzdes" ${txt} | grep numbering )
        echo $frompos
        for line in ${frompos}
        do
            echo $line
            if [ $( echo ${line} | grep -o '+' | wc -l ) -ge 2 ]; 
            then  
                echo "inside"
                newline=$line
                array=();
                while [[ $newline ]];
                do 
                    array+=("${newline%%'+'*}");
                    newline=${newline#*'+'};
                done
                declare -p array
                for x in "${array[@]}";
                do
                    echo ${x} >> ${POS}
                done
            fi
        done
fi
exit
for uaa in AZC AZT; 
do
    while read a;
    do
        UAA="inputs/${uaa}.params"
        #python ${scripts}/mutation_design.py -p ${PDB}.pdb -l protein -res ${TEV} \
        #        -num ${a} -r 6.0 -uaa ${UAA} -dir ${PDB}_AZO -n ${N} -rem
        nohup python ${scripts}/mutation_design.py -p ${PDB}.pdb -l protein -res ${TEV} \
                -num ${a} -r 6.0 -uaa ${UAA} -dir ${PDB}_AZO -n ${N} -rem &
        #exit
        #sleep 1
        #exit
        running=$( top -bn 1 | grep ' R ' | wc -l )
        while [ ${running} -ge 40 ]
        do
            echo "sleeping ... for 120 seconds"
            sleep 120
            running=$( top -bn 1 | grep ' R ' | wc -l )
        done
    done<${POS}
done

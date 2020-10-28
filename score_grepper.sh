#!/bin/bash

grep 'label ' *${1}*pdb | head -n1 > ${1}.scores
grep 'pose ' *${1}*.pdb >> ${1}.scores
sed -e 's/ /,/g' < ${1}.scores > ${1}.scores.csv

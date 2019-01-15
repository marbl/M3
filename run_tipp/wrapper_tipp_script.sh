#!/bin/bash

# source /fs/cbcb-lab/mpop/m3_taxonomy_workshop/test_user_profile.sh

TIPP="/fs/m3taxworkshop/bin/run_tipp.py"
REFPKG="/fs/m3taxworkshop/software/tipp/refpkg/RDP_2016_Bacteria.refpkg/"
UTILS_DIR="/fs/m3taxworkshop/software/tipp_utils/utils/"

A=100
P=1000


dat=$1
dat_base=`basename $dat`
dat_dir=`dirname $dat`
dat_rev=`echo ${dat_base%%.*}`
dat_rev_ext=`echo ${dat_base#*.}`
dat_name=`echo ${dat%%.*}`
dat_rev_name=$(pwd)/${dat_rev}_revcompl.${dat_rev_ext}
# echo $dat_rev_name
nam="$2_tipp_A${A}_P${P}"
nam_rev="${nam}_revcompl"
# echo $nam $nam_rev

#First get the reverse complement of the reads
echo "Reverse complementing fasta file - ${dat}"
python ${UTILS_DIR}/reverse_complement_fasta.py -i ${dat} -o ${dat_rev_name}
if [ $? -eq 0 ]; then
    echo OK
else
    echo FAIL
    echo "Error while reverse complementing fasta file"
    exit 
fi

echo "Running TIPP on ${dat}"

python ${TIPP} -t ${REFPKG}/pasta.taxonomy \
               -a ${REFPKG}/pasta.fasta \
               -r ${REFPKG}/RAxML_info.taxonomy \
               -tx ${REFPKG}/taxonomy.table \
               -txm ${REFPKG}/species.mapping \
               -A ${A} -P ${P} \
               -at 0.95 \
               -pt 0.95 \
               -f ${dat} \
               -o ${nam} -d $(pwd) \
               --tempdir $(pwd)/tmp \
               --cpu 2

if [ $? -eq 0 ]; then
    echo OK
else
    echo FAIL
    echo "Error while running TIPP on ${dat_r}"
    exit 
fi

echo "Running TIPP on ${dat_rev_name}"
python ${TIPP} -t ${REFPKG}/pasta.taxonomy \
               -a ${REFPKG}/pasta.fasta \
               -r ${REFPKG}/RAxML_info.taxonomy \
               -tx ${REFPKG}/taxonomy.table \
               -txm ${REFPKG}/species.mapping \
               -A ${A} -P ${P} \
               -at 0.95 \
               -pt 0.95 \
               -f ${dat_rev_name} \
               -o ${nam_rev} -d $(pwd) \
               --tempdir $(pwd)/tmp \
                --cpu 2

if [ $? -eq 0 ]; then
    echo OK
else
    echo FAIL
    echo "Error while running TIPP on ${dat_rev_name}"
    exit 
fi

# Combine the classification from both forward and reverse reads
python ${UTILS_DIR}/restructure_tipp_classification.py -i ${nam}_classification.txt -r ${nam_rev}_classification.txt -o ${nam}_merged
if [ $? -eq 0 ]; then
    echo OK
else
    echo FAIL
    echo "Error while merging TIPP outputs for forward and reverse reads"
    exit 
fi
python ${UTILS_DIR}/format_assignment_output.py -i ${nam}_merged.csv -o ${nam}_merged
if [ $? -eq 0 ]; then
    echo OK
else
    echo FAIL
    echo "Error while formatting the final taxonomic assignment output"
    exit 
fi


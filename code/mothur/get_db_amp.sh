#!/usr/bin/sh
DB_ROOT=${1}
DB_AMP_ROOT=${2}
DB_AMP_START=${3}
DB_AMP_END=${4}
MOTHUR=${5}

## Extracting Amplicon Target Region
DB_FILES="fasta=${DB_ROOT}.align, taxonomy=${DB_ROOT}.tax"
PCR_COORD="start=${DB_AMP_START}, end=${DB_AMP_END}"
${MOTHUR} "#pcr.seqs(${DB_FILES}, ${PCR_COORD}, keepdots=F, processors=4); unique.seqs(fasta=current)"


## cleanup
mv ${DB_ROOT}.pcr.unique.align ${DB_AMP_ROOT}.align
rm ${DB_ROOT}.pcr.align
rm ${DB_ROOT}.pcr.tax # not keeping taxonomy for amp region may want to use for tax
rm ${DB_ROOT}.pcr.names

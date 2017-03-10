#!/usr/bin/sh

TRAIN_ROOT=${1}
DB_AMP_ROOT=${2}
TRAIN_AMP_ROOT=${3}
TRAIN_AMP_START=${4}
TRAIN_AMP_END=${5}
MOTHUR=${6}


TRAIN_SEQ=${TRAIN_ROOT}.fasta
DB_ALIGN=${DB_AMP_ROOT}.align
TRAIN_TAX=${TRAIN_ROOT}.tax
COORDS="start=${TRAIN_AMP_START}, end=${TRAIN_AMP_END}"
${MOTHUR} "#align.seqs(fasta=${TRAIN_SEQ}, reference=${DB_ALIGN}, processors=4);\
        screen.seqs(fasta=current, taxonomy=${TRAIN_TAX}, ${COORDS});\
        degap.seqs(fasta=current)"

# ## Extracting Amplicon Target Region
# TRAIN_SEQ=${TRAIN_ROOT}.fasta
# TRAIN_TAX=${TRAIN_ROOT}.tax
# PCR_COORD="start=${TRAIN_AMP_START}, end=${TRAIN_AMP_END}"
# ./mothur "#pcr.seqs(fasta=${TRAIN_SEQ}, taxonomy=${TRAIN_TAX}, ${PCR_COORD}); degap.seqs(current)"

## Clean-up
mv ${TRAIN_ROOT}.good.ng.fasta ${TRAIN_AMP_ROOT}.fasta
mv ${TRAIN_ROOT}.good.tax ${TRAIN_AMP_ROOT}.tax
rm ${TRAIN_ROOT}.bad.accnos
rm ${TRAIN_ROOT}.flip.accnos


## Notes
# screens the training set for sequences that cover the ampicon region,
# might want to only include the amplicon region instead?

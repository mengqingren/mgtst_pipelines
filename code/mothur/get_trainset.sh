TRAIN_TAX=$1 #data/references/trainset10_082014.pds.tax
TRAIN_SEQ=$2 #data/references/trainset10_082014.pds.fasta
REF_ALIGN=$3 #data/refrences/silva.v34.align
POS_REF=$4 #data/references/start_stop.positions
OUTPUT=$5 #data/references


run_mothur (){
    START=$1
    END=$2
    REGION=$3

    echo "Running align.seqs"

    # Error with start and end - filters all sequences, not sure needed since aligning to truncated ref
    # ./mothur "#align.seqs(fasta=$TRAIN_SEQ, reference=$REF_ALIGN, flip=t, processors=8);
    #     screen.seqs(fasta=current, taxonomy=$TRAIN_TAX, start=$START, end=$END);
    #     degap.seqs(fasta=current)"

    ./mothur "#align.seqs(fasta=$TRAIN_SEQ, reference=$REF_ALIGN, flip=t, processors=8);
        screen.seqs(fasta=current, taxonomy=$TRAIN_TAX);
        degap.seqs(fasta=current);list.seqs(fasta=current)"


    mv $OUTPUT/trainset10_082014.pds.good.ng.fasta $OUTPUT/trainset10_082014.v34.fasta
    mv $OUTPUT/trainset10_082014.pds.good.tax $OUTPUT/trainset10_082014.v34.tax
    rm $OUTPUT/trainset10_082014.pds.align*
    rm $OUTPUT/trainset10_082014.pds.bad.accnos
    rm $OUTPUT/trainset10_082014.pds.flip.accnos
}

V34_START=$(grep "V34" $POS_REF | cut -f 2 -d " ")
V34_END=$(grep "V34" $POS_REF | cut -f 3 -d " ")
run_mothur $V34_START $V34_END v34

# ./mothur "#align.seqs(fasta=data/references/trainset10_082014.pds.fasta, reference=data/references/silva.v34.align, flip=t, processors=8);screen.seqs(fasta=current, taxonomy=data/references/trainset10_082014.pds.tax)"

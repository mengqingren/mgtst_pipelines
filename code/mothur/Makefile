#========================================================================================================
#
# VARIABLES
#
#========================================================================================================

## Directories
REF=reference
PROCESS=mothur/
FQ=data/
CODE=code/mothur/

## Amplicon Target Region
F=CCTACGGGNGGCWGCAG
R=GGACTACHVGGGTWTCTAAT
REGION=V34

## Reference Database
DB_GET=$(CODE)/get_silva.sh
DB_ROOT=$(REF)/silva.nr_v119

## Reference Amplicon Target Region
DB_AMP_START=6428
DB_AMP_END=23444
DB_AMP_ROOT=$(DB_ROOT).$(REGION)

## Trainset Database
TRAIN_GET=$(CODE)/get_trainset_v10_082014.sh
TRAIN_ROOT=$(REF)/trainset10_082014.pds

## Trainset Amplicon Target Region
TRAIN_AMP_START=2
TRAIN_AMP_END=17016
TRAIN_AMP_ROOT=$(REF)/trainset10_082014.$(REGION)
TRAIN_SEQ=$(TRAIN_AMP_ROOT).fasta
TRAIN_TAX=$(TRAIN_AMP_ROOT).tax

## Mothur
PROC=processors=8
# clustering requires more memory ~4Gb per thread
CLUSTER_PROC=processors=1
MOTHUR=bin/mothur

## Data Set Processing
# root name of mothur input file defining fq files and groups
PROCESS_ROOT=mgtst

## Mothur pipeline file names
FQ_FILES=$(PROCESS)/$(PROCESS_ROOT).files
ROOT=$(PROCESS)/$(PROCESS_ROOT)

# make contigs
CONTIG_FASTA=$(ROOT).trim.contigs.good.fasta
CONTIG_GROUP=$(ROOT).contigs.good.groups

# count contigs
CONTIG_NAMES=$(ROOT).trim.contigs.good.names
CONTIG_UNIQUE=$(ROOT).trim.contigs.good.unique.fasta
CONTIG_COUNT=$(ROOT).trim.contigs.good.count_table

# Align contigs to ref
ALIGN=$(ROOT).trim.contigs.good.unique.align
ALIGN_CLEAN=$(ROOT).trim.contigs.good.unique.good.filter.unique.fasta
ALIGN_COUNT=$(ROOT).trim.contigs.good.unique.good.filter.count_table

# Cluster
CLUSTER=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.fasta
CLUSTER_COUNT=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.count_table

# Chimeria removal
CHIMERA=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta
CHIMERA_COUNT=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table

# Classificaiton
CLASSIFY=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.V34.wang.taxonomy
CLUSTER_TAX=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.list

#========================================================================================================
#
# Pipeline
#
#========================================================================================================

######
#
# Get Reference Database
#
######

$(DB_ROOT).align :
	bash $(DB_GET) $(REF)

$(DB_AMP_ROOT).align : $(DB_ROOT).align
	bash $(CODE)/get_db_amp.sh $(DB_ROOT) $(DB_AMP_ROOT) $(DB_AMP_START) $(DB_AMP_END) $(MOTHUR)

.PHONY :  get.silva
get.silva :  $(DB_AMP_ROOT).align

######
#
# Get Trainset
#
######

$(TRAIN_ROOT).fasta :
	bash $(TRAIN_GET) $(REF)


#get the V34 region of the RDP training set
$(TRAIN_SEQ) : $(TRAIN_ROOT).fasta $(TRAIN_ROOT).tax $(DB_AMP_ROOT).align
	bash $(CODE)/get_train_amp.sh $(TRAIN_ROOT) $(DB_AMP_ROOT) $(TRAIN_AMP_ROOT)\
        $(TRAIN_AMP_START) $(TRAIN_AMP_END) $(MOTHUR)

.PHONY : get.trainset
get.trainset : $(TRAIN_ROOT).fasta $(TRAIN_SEQ)

######
#
# Generate Contigs
#
######

$(PROCESS) :
	mkdir $(PROCESS)

# Merging read pairs and removing low quality contigs
$(CONTIG_FASTA) : $(FQ_FILES) $(FQ) $(PROCESS)
	$(MOTHUR) "#make.contigs(file=$(FQ_FILES), inputdir=$(FQ)/, outputdir=$(PROCESS)/, $(PROC));\
        screen.seqs(fasta=current, group=current, maxambig=0, maxlength=500, maxhomop=8)"

$(CONTIG_UNIQUE) : $(CONTIG_FASTA)
	$(MOTHUR) "#unique.seqs(fasta=$(CONTIG_FASTA));\
         count.seqs(name=current, group=$(CONTIG_GROUP))"

.PHONY : make.contigs
make.contigs : $(CONTIG_FASTA) $(CONTIG_UNIQUE)

######
#
# Align Contigs to Ref Alignment
#
######

$(ALIGN) : $(CONTIG_UNIQUE) $(DB_AMP_ROOT).align
	$(MOTHUR) "#align.seqs(fasta=$(CONTIG_UNIQUE), reference=$(DB_AMP_ROOT).align, flip=t, $(PROC))"

$(ALIGN_CLEAN): $(ALIGN)
	$(MOTHUR) "#screen.seqs(fasta=$(ALIGN), count=$(CONTIG_COUNT), minlength=400, maxlength=500);\
        filter.seqs(fasta=current, vertical=T, trump=.);\
        unique.seqs(fasta=current, count=current)"

.PHONY : align.contigs
align.contigs :  $(ALIGN) $(ALIGN_CLEAN)

######
#
# Cluster Contigs
#
######

$(CLUSTER) : $(ALIGN_CLEAN)
	$(MOTHUR) "#pre.cluster(fasta=$(ALIGN_CLEAN), count=$(ALIGN_COUNT), diffs=5, $(PROC))"

$(CHIMERA) : $(CLUSTER)
	$(MOTHUR) "#chimera.uchime(fasta=$(CLUSTER), count=$(CLUSTER_COUNT), dereplicate=T, $(CLUSTER_PROC));\
        remove.seqs(fasta=current, accnos=current)"

$(CLASSIFY) : $(CHIMERA) $(TRAIN_SEQ)
	$(MOTHUR) "#classify.seqs(fasta=$(CHIMERA), count=$(CHIMERA_COUNT), reference=$(TRAIN_SEQ), taxonomy=$(TRAIN_TAX), cutoff=80, $(PROC))"

$(CLUSTER_TAX): $(CHIMERA) $(CLASSIFY)
	$(MOTHUR) "#cluster.split(fasta=$(CHIMERA), count=$(CHIMERA_COUNT), taxonomy=$(CLASSIFY), splitmethod=classify, taxlevel=4, $(CLUSTER_PROC));\
        make.shared(list=current,count=current, label=0.03);\
        classify.otu(list=current, count=current, taxonomy=current, label=0.03)"

.PHONY : cluster.contigs
cluster.contigs : $(CLUSTER) $(CHIMERA) $(CLASSIFY) $(CLUSTER_TAX)

############################################################################################################
#
#   Running Full Pipeline
#
############################################################################################################

all: cluster.contigs

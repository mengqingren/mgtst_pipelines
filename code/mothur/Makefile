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
DB_AMP_START=6388
DB_AMP_END=25316
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
MAKE_FQ_FILE=$(CODE)/make_mgtst_files.R
FQ_FILES=$(PROCESS)/$(PROCESS_ROOT).files
ROOT=$(PROCESS)/$(PROCESS_ROOT)

# make contigs
CONTIG_FASTA=$(ROOT).trim.contigs.good.unique.fasta
CONTIG_COUNT=$(ROOT).trim.contigs.good.count_table

# Align contigs to ref
ALIGN=$(ROOT).trim.contigs.good.unique.align
ALIGN_COUNT=$(ROOT).trim.contigs.good.good.count_table
ALIGN_CLEAN=$(ROOT).trim.contigs.good.unique.good.filter.unique.fasta
ALIGN_CLEAN_COUNT=$(ROOT).trim.contigs.good.unique.good.filter.count_table

# Unclustered
UNCLUSTERED_SEQ=$(ROOT).trim.contigs.good.unique.good.ng.good.fasta
UNCLUSTERED_COUNT=$(ROOT).trim.contigs.good.good.count_table
UNCLUSTERED_TAX=$(ROOT).trim.contigs.good.unique.ng.good.V34.wang.taxonomy

# Precluster
PRECLUSTER_FASTA=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.fasta
PRECLUSTER_COUNT=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.count_table

# Chimeria removal
CHIMERA_FASTA=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta
CHIMERA_COUNT=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table

# Classification
CLASSIFY=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.V34.wang.taxonomy

# Cluster
CLUSTER_DIST=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.dist
CLUSTER_LIST=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.list
CLUSTER_SHARED=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.shared
CLUSTER_TAX=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.cons.taxonomy


# Representative Seqs
REP_FASTA=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.rep.fasta
REP_DIST=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.list.0.03.rep.dist
REP_NG_FASTA=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.rep.ng.fasta
REP_TREE=$(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.rep.phylip.tre

############################################################################################################
#
#   Running Full Pipeline
#
############################################################################################################

all: $(UNCLUSTERED_SEQ) $(UNCLUSTERED_TAX) $(CLUSTER_SHARED) $(CLUSTER_TAX) $(REP_NG_FASTA) $(REP_TREE)

.PHONY: clean
clean:
	rm $(ROOT).trim.contigs.good.unique.fasta
	rm $(ROOT).trim.contigs.good.unique.good.filter.count_table
	rm $(ROOT).trim.contigs.good.unique.good.filter.unique.fasta
	rm $(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.count_table
	rm $(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.chimeras
	rm $(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table
	rm $(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.fasta
	rm $(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.V34.wang.taxonomy
	rm $(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.dist
	rm $(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta
	rm $(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.rep.fasta
	rm $(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.rep.phylip.dist
	rm $(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.list
	rm $(ROOT).trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.sensspec




#========================================================================================================
#
# Pipeline Steps
#
#========================================================================================================

######
#
# Get Reference Database
#
######

$(DB_ROOT).align:
	bash $(DB_GET) $(REF)

$(DB_AMP_ROOT).align: $(DB_ROOT).align
	bash $(CODE)/get_db_amp.sh $(DB_ROOT) $(DB_AMP_ROOT) $(DB_AMP_START) $(DB_AMP_END) $(MOTHUR)

.PHONY:  get.silva
get.silva:  $(DB_AMP_ROOT).align

######
#
# Get Trainset
#
######

$(TRAIN_ROOT).fasta $(TRAIN_ROOT).tax:
	bash $(TRAIN_GET) $(REF)


#get the V34 region of the RDP training set
$(TRAIN_SEQ): $(TRAIN_ROOT).fasta $(TRAIN_ROOT).tax $(DB_AMP_ROOT).align
	bash $(CODE)/get_train_amp.sh $(TRAIN_ROOT) $(DB_AMP_ROOT) $(TRAIN_AMP_ROOT)\
        $(TRAIN_AMP_START) $(TRAIN_AMP_END) $(MOTHUR)

.PHONY: get.trainset
get.trainset: $(TRAIN_ROOT).fasta $(TRAIN_SEQ) $(TRAIN_ROOT).tax

######
#
# Generate Contigs
#
######

# generate mgtst.files
$(FQ_FILES):
	R CMD BATCH $(MAKE_FQ_FILE)

# Merging read pairs and removing low quality contigs
$(CONTIG_FASTA) $(CONTIG_COUNT): $(FQ_FILES)
	$(MOTHUR) "#make.contigs(file=$(FQ_FILES), inputdir=$(FQ), outputdir=$(PROCESS), $(PROC));\
              screen.seqs(fasta=current, group=current, maxambig=0, maxlength=500, maxhomop=8);\
				      unique.seqs(fasta=current);\
				      count.seqs(name=current, group=current)"
	rm mothur/mgtst.scrap.contigs.fasta mothur/mgtst.scrap.contigs.qual mothur/mgtst.trim.contigs.fasta mothur/mgtst.trim.contigs.qual mothur/mgtst.trim.contigs.good.fasta

######
#
# Align Contigs to Ref Alignment
#
######

$(ALIGN): $(CONTIG_FASTA) $(DB_AMP_ROOT).align
	$(MOTHUR) "#align.seqs(fasta=$(CONTIG_FASTA), reference=$(DB_AMP_ROOT).align, flip=t, $(PROC))"

$(ALIGN_CLEAN) $(ALIGN_COUNT): $(ALIGN)
	$(MOTHUR) "#screen.seqs(fasta=$(ALIGN), count=$(CONTIG_COUNT), minlength=400, maxlength=500);\
              filter.seqs(fasta=current, vertical=T, trump=.);\
              unique.seqs(fasta=current, count=current)"
	rm mothur/mgtst.trim.contigs.good.unique.good.align mothur/mgtst.trim.contigs.good.unique.good.filter.fasta

######
#
# Unclustered
#
######

$(UNCLUSTERED_SEQ) $(UNCLUSTERED_TAX):	$(ALIGN) $(ALIGN_COUNT) $(TRAIN_SEQ) $(TRAIN_TAX)
	$(MOTHUR) "#degap.seqs(fasta=$(ALIGN));\
							screen.seqs(fasta=current, count=$(CONTIG_COUNT), minlength=430);\
							classify.seqs(fasta=current, count=current, reference=$(TRAIN_SEQ), taxonomy=$(TRAIN_TAX), cutoff=80, $(PROC))"
	rm mothur/mgtst.trim.contigs.good.unique.align mothur/mgtst.trim.contigs.good.unique.ng.fasta

######
#
# Cluster Contigs
#
######

$(PRECLUSTER_FASTA) $(PRECLUSTER_COUNT): $(ALIGN_CLEAN) $(ALIGN_CLEAN_COUNT)
	$(MOTHUR) "#pre.cluster(fasta=$(ALIGN_CLEAN), count=$(ALIGN_CLEAN_COUNT), diffs=5, $(PROC))"

$(CHIMERA_FASTA) $(CHIMERA_COUNT): $(PRECLUSTER_FAST) $(PRECLUSTER_COUNT)
	$(MOTHUR) "#chimera.uchime(fasta=$(PRECLUSTER_FASTA), count=$(PRECLUSTER_COUNT), dereplicate=T, $(CLUSTER_PROC));\
              remove.seqs(fasta=current, accnos=current, count=current)"

$(CLASSIFY): $(CHIMERA_FASTA) $(CHIMERA_COUNT) $(TRAIN_SEQ) $(TRAIN_TAX)
	$(MOTHUR) "#classify.seqs(fasta=$(CHIMERA_FASTA), count=$(CHIMERA_COUNT), reference=$(TRAIN_SEQ), taxonomy=$(TRAIN_TAX), cutoff=80, $(PROC))"

$(CLUSTER_FASTA) $(CLUSTER_TAX) $(CLUSTER_SHARED): $(CHIMERA_FASTA) $(CHIMERA_COUNT) $(CLASSIFY)
	$(MOTHUR) "#cluster.split(fasta=$(CHIMERA_FASTA), count=$(CHIMERA_COUNT), taxonomy=$(CLASSIFY), splitmethod=classify, taxlevel=4, $(CLUSTER_PROC));\
              make.shared(list=current, count=current, label=0.03);\
              classify.otu(list=current, count=current, taxonomy=$(CLASSIFY), label=0.03)"

######
#
# Representative Sequences
#
######

$(REP_FASTA) $(REP_NG_FASTA): $(CLUSTER)
	$(MOTHUR) "#get.oturep(column=$(CLUSTER_DIST), list=$(CLUSTER_LIST), count=$(CHIMERA_COUNT), fasta=$(CHIMERA_FASTA));\
				      degap.seqs(fasta=current)"

$(REP_TREE): $(REP_FASTA)
	$(MOTHUR) "#dist.seqs(fasta=$(REP_FASTA), output=lt, $(PROC));\
				      clearcut(phylip=current)"

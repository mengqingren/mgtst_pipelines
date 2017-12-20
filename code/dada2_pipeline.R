library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(parallel); packageVersion("parallel")
library(phyloseq); packageVersion("phyloseq")
library(msa); packageVersion("msa")
library(phangorn); packageVersion("phangorn")
cores <- 7

fastq_dir <- "data/"
fns <- list.files(fastq_dir,recursive = TRUE)
fastqs <- fns[grepl("*001.fastq.gz$", fns)]
fastqs <- sort(fastqs) # Sort should keep them paired in order
fnFs <- fastqs[grepl("_R1", fastqs)]
fnRs <- fastqs[grepl("_R2", fastqs)]


## Quality filter --------------------------------------------------------------
filt_dir <- "dada2/"
filtFs <- paste0(filt_dir, sapply(strsplit(fnFs, "\\."), 
                                  `[`, 1), "_filt.fastq.gz")
filtRs <- paste0(filt_dir, sapply(strsplit(fnRs, "\\."), 
                                  `[`, 1), "_filt.fastq.gz")
for (i in seq_along(fnFs)) {
      if (!(file.exists(filtFs[i]) && file.exists(filtRs[i]))) {
            fastqPairedFilter(paste0(fastq_dir, c(fnFs[i], fnRs[i])), 
                              c(filtFs[i], filtRs[i]),
                              maxN = 0, maxEE = 4, truncQ = 2, 
                              truncLen = c(280,220), trimLeft = c(10, 10), 
                              compress = TRUE, verbose = TRUE)
      }
}


## Dereplicate ----------------------------------------------------------------
# Name the derep-class objects by the sample names
sam_names <- sapply(strsplit(fnFs, "/"), tail, n = 1)
sam_names <- sapply(strsplit(sam_names, "_"), `[`, 1)

if (!file.exists("dada2/derepFs.rds")) {
      derepFs <- lapply(filtFs, derepFastq, verbose = TRUE)
      names(derepFs) <- sam_names
      saveRDS(derepFs, "dada2/derepFs.rds")
}else{
      derepFs <- readRDS("dada2/derepFs.rds")
}

if (!file.exists("dada2/derepRs.rds")) {
      derepRs <- lapply(filtRs, derepFastq, verbose = TRUE)
      names(derepRs) <- sam_names
      saveRDS(derepRs, "dada2/derepRs.rds")
}else{
      derepRs <- readRDS("dada2/derepRs.rds")
}

## Sequence Inference ----------------------------------------------------------
##### Forwared Reads
dadaF_file <- "dada2/dadaFs.rds"
if (!file.exists(dadaF_file)) {
  dadaFs <- dada(derepFs, 
                 err = inflateErr(tperr1,3), 
                 errorEstimationFunction = loessErrfun,
                 selfConsist = TRUE,
                 pool = FALSE, 
                 multithread = TRUE) 
  saveRDS(dadaFs, dadaF_file)
}else{
  dadaFs <- readRDS(dadaF_file)
}

##### Reverse Reads
dadaR_file <- "dada2/dadaRs.rds"
if (!file.exists(dadaR_file)) {
  dadaRs <- dada(derepRs, 
                 err = inflateErr(tperr1,3), 
                 errorEstimationFunction = loessErrfun,
                 selfConsist = TRUE,
                 pool = FALSE, 
                 multithread = TRUE) 
  saveRDS(dadaRs, dadaR_file)
}else{
  dadaRs <- readRDS(dadaR_file)
}

## Merge Paired Reads ---------------------------------------------------------
mergers_file <- "dada2/mergers.rds"
if (!file.exists(mergers_file)) {
      mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
      saveRDS(mergers, mergers_file)
}else{
      mergers <- readRDS(mergers_file)
}

## Construct Sequence Table ---------------------------------------------------
seqtab_file <- "dada2/seqtab.rds"
if (!file.exists(seqtab_file)) {
      seqtab <- makeSequenceTable(mergers)
      saveRDS(seqtab, seqtab_file)
}else{
      seqtab <- readRDS(seqtab_file)
}

## Remove Chimeras ------------------------------------------------------------
seqtab_nochim_file <- "dada2/seqtab_nochim.rds"
if (!file.exists(seqtab_nochim_file)) {
      seqtab_nochim <- removeBimeraDenovo(seqtab, verbose = TRUE)
      saveRDS(seqtab_nochim, seqtab_nochim_file)
}else{
      seqtab_nochim <- readRDS(seqtab_nochim_file)
}

## Assign Taxonomy -------------------------------------------------------------
taxa_file <- "dada2/taxa.rds"
if (!file.exists(taxa_file)) {
  taxa <- assignTaxonomy(seqtab_nochim, 
                         "reference/silva_nr_v123_train_set.fa.gz")
  colnames(taxa) <- paste0("Rank",1:6)
  saveRDS(taxa, taxa_file)
}else{
  seqtab_nochim <- readRDS(taxa_file)
}

## Representative Sequences ----------------------------------------------------
seq_file <- "dada2/sv_seqs.fasta"
if (!file.exists(seq_file)) {
    ## Extract sequences and name features
    sv_seqs <- getSequences(seqtab_nochim)
    names(sv_seqs) <- paste0("SV",1:length(sv_seqs))
    
    ## Write sequences to file
    DNAStringSet(sv_seqs) %>% writeXStringSet()
} else {
    sv_seqs <- readRDS(seq_file)
}


## Phylogenetic Tree -----------------------------------------------------------
##
## Tree construction methods based on Callahan BJ, Sankaran K, Fukuyama JA et
## al. Bioconductor Workflow for Microbiome Data Analysis: from raw reads to
## community analyses [version 2; referees: 3 approved]. F1000Research 2016,
## 5:1492 (doi: 10.12688/f1000research.8986.2)
##
## The phangorn package is then used to construct a phylogenetic tree. Here we
## first construct a neighbor-joining tree, and then fit a GTR+G+I maximum
## likelihood tree using the neighbor-joining tree as a starting point.

tree_file <- "dada2/dada_tree_GTR.rds" 
if (!file.exists(tree_file)) {
    ## Multiple sequence alignment 
    mult <- msa(sv_seqs, method = "ClustalW", type = "dna", order = "input")
   
    
    ## Generate phylogenetic tree
    phang.align <- as.phyDat(mult, type = "DNA", names = names(sv_seqs))
    dm <- dist.ml(phang.align)
    treeNJ <- NJ(dm) # Note, tip order != sequence order
    fit = pml(treeNJ, data = phang.align)
    
    ## negative edges length changed to 0!
    fitGTR <- update(fit, k = 4, inv = 0.2)
    fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                        rearrangement = "stochastic", 
                        control = pml.control(trace = 0))
    
    saveRDS(fitGTR, tree_file)
}
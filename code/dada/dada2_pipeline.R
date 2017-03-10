library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(parallel); packageVersion("parallel")
cores <- 7

setwd("~/Projects/16S_etec_mix_study/analysis/pipelines/dada2/")


project_dir <- "~/Projects/16S_etec_mix_study"
fastq_dir <- file.path(project_dir, "data/fastq/jhu_run2/")
fns <- list.files(fastq_dir,recursive = TRUE)
fastqs <- fns[grepl("*001.fastq.gz$", fns)]
fastqs <- sort(fastqs) # Sort should keep them paired in order
fnFs <- fastqs[grepl("_R1", fastqs)]
fnRs <- fastqs[grepl("_R2", fastqs)]


## Quality filter --------------------------------------------------------------
filt_dir <- "processed_data/"
filtFs <- paste0(filt_dir, sapply(strsplit(fnFs, "\\."), 
                                  `[`, 1), "_filt.fastq.gz")
filtRs <- paste0(filt_dir, sapply(strsplit(fnRs, "\\."), 
                                  `[`, 1), "_filt.fastq.gz")
for(i in seq_along(fnFs)) {
      if(!(file.exists(filtFs[i]) && file.exists(filtRs[i]))){
            fastqPairedFilter(paste0(fastq_dir, c(fnFs[i], fnRs[i])), 
                              c(filtFs[i], filtRs[i]),
                              maxN=0, maxEE=4, truncQ=2,truncLen = c(290,220),
                              trimLeft=c(10, 10), compress=TRUE, verbose=TRUE)
      }
}


## Dereplicate ----------------------------------------------------------------
# Name the derep-class objects by the sample names
sam_names <- sapply(strsplit(fnFs, "/"), tail, n=1)
sam_names <- sapply(strsplit(sam_names, "_"), `[`, 1)

if(!(file.exists("processed_data/derepFs-2016-11-07.rds"))){
      derepFs <- lapply(filtFs, derepFastq, verbose=TRUE)
      names(derepFs) <- sam_names
      saveRDS(derepFs, "processed_data/derepFs-2016-11-07.rds")
}else{
      derepFs <- readRDS("processed_data/derepFs-2016-11-07.rds")
}

if(!(file.exists("processed_data/derepRs-2016-11-07.rds"))){
      derepRs <- lapply(filtRs, derepFastq, verbose=TRUE)
      names(derepRs) <- sam_names
      saveRDS(derepRs, "processed_data/derepRs-2016-11-07.rds")
}else{
      derepRs <- readRDS("processed_data/derepRs-2016-11-07.rds")
}

### Single Inference ----------------------------------------------------------
#### Infer Errors -------------------------------------------------------------
##### Forwared Reads
dadaF_file <- "processed_data/dadaFs-single-inference-2016-11-07.rds"
if(!(file.exists(dadaF_file))){
  dadaFs <- dada(derepFs, 
                 err=inflateErr(tperr1,3), 
                 errorEstimationFunction=loessErrfun,
                 selfConsist = TRUE,
                 pool = FALSE, 
                 multithread = TRUE) 
  saveRDS(dadaFs, dadaF_file)
}else{
  dadaFs <- readRDS(dadaF_file)
}

##### Reverse Reads
dadaR_file <-"processed_data/dadaRs-single-inference-2016-11-07.rds"
if(!(file.exists(dadaR_file))){
  dadaRs <- dada(derepRs, 
                 err=inflateErr(tperr1,3), 
                 errorEstimationFunction=loessErrfun,
                 selfConsist = TRUE,
                 pool = FALSE, 
                 multithread = TRUE) 
  saveRDS(dadaRs, dadaR_file)
}else{
  dadaRs <- readRDS(dadaR_file)
}

## Merge Paired Reads ---------------------------------------------------------
mergers_file <-"processed_data/mergers-single-inference-2016-11-07.rds"
if(!(file.exists(mergers_file))){
      mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
      saveRDS(mergers, mergers_file)
}else{
      mergers <- readRDS(mergers_file)
}

## Construct Sequence Table ---------------------------------------------------
seqtab_file <-"processed_data/seqtab-single-inference-2016-11-07.rds"
if(!(file.exists(seqtab_file))){
      seqtab <- makeSequenceTable(mergers)
      saveRDS(seqtab, seqtab_file)
}else{
      seqtab <- readRDS(seqtab_file)
}

## Remove Chimeras ------------------------------------------------------------
seqtab_nochim_file <-"processed_data/seqtab_nochim-single-inference-2016-11-07.rds"
if(!(file.exists(seqtab_nochim_file))){
      seqtab_nochim <- removeBimeraDenovo(seqtab, verbose = TRUE)
      saveRDS(seqtab_nochim, seqtab_nochim_file)
}else{
      seqtab_nochim <- readRDS(seqtab_nochim_file)
}

## Assign Taxonomy
taxa <- assignTaxonomy(seqtab_nochim, "silva_nr_v123_train_set.fa.gz")
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

## Generate Phyloseq Object
library(phyloseq); packageVersion("phyloseq")

# Make sample_data
library(mgtst)
library(dplyr)
library(tidyr)
library(stringr)

data(sample_sheet)
sample_sheet_ns <- sample_sheet %>% mutate(pos_ns = str_replace(pos, "_",""))
sample_id <- data_frame(seq_ids = rownames(seqtab)) %>%
      separate(col = seq_ids, into = c("pcr_16S_plate","pos_ns"), sep = "-") %>%
      mutate(pcr_16S_plate = as.numeric(pcr_16S_plate))

sample_df <- left_join(sample_id, sample_sheet_ns) %>%
      filter(seq_lab == "JHU", barcode_lab == "JHU") %>%
      mutate(pcr_16S_plate = as.character(pcr_16S_plate)) %>% as.data.frame()

rownames(sample_df) <- rownames(seqtab_nochim)

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE), 
               sample_data(sample_df), 
               tax_table(taxa))
phyloseq_file <-"processed_data/phyloseq-single-inference-2016-11-07.rds"
saveRDS(ps, phyloseq_file)
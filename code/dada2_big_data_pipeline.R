## quality filter and trim
library(dada2, lib.loc = "~/R/x86_64-pc-linux-gnu-library/"); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(stringr); packageVersion("stringr")
library(parallel); packageVersion("parallel")
library(phyloseq, lib.loc = "~/R/x86_64-pc-linux-gnu-library/"); packageVersion("phyloseq")
library(msa, lib.loc = "~/R/x86_64-pc-linux-gnu-library/"); packageVersion("msa")
library(phangorn, lib.loc = "~/R/x86_64-pc-linux-gnu-library/"); packageVersion("phangorn")

## Remove primers 
fastq_dir <- "data"
fastqFs <- list.files(fastq_dir, pattern = ".*R1.*fastq.gz", full.names = TRUE)
fastqRs <- list.files(fastq_dir, pattern = ".*R2.*fastq.gz", full.names = TRUE)

cutadapt_primers <- "cutadapt -g CCTACGGGNGGCWGCAG -G GGACTACHVGGGTWTCTAAT"

for (infile_R1 in fastqFs) {
    
    outfile_R1 <- basename(infile_R1) %>%
        str_replace(".fastq.gz", "_primer_cut.fastq.gz") %>%
        file.path("dada2", .)
    
    infile_R2 <- str_replace(infile_R1, "R1", "R2")
    outfile_R2 <- basename(infile_R2) %>%
        str_replace(".fastq.gz", "_primer_cut.fastq.gz") %>%
        file.path("dada2", .)
    
    if (!(file.exists(outfile_R1) & file.exists(outfile_R2))) {
        cutadapt_cmd <-
            paste(cutadapt_primers,
                  "-o",
                  outfile_R1,
                  "-p",
                  outfile_R2,
                  infile_R1,
                  infile_R2)
        system(command = cutadapt_cmd, wait = FALSE)
    }
}

### Filtered forward and reverse files go into the primer_cut/filtered_F and
### primer_cut/filtered_R subdirectory respectively 

### Input file directory
fastq_dir <- "dada2"
filtpathF <- file.path("dada2", "filtered_F")
filtpathR <- file.path("dada2", "filtered_R")

fastqFs <- sort(list.files(fastq_dir, pattern = ".*R1.*_primer_cut.fastq.gz"))
fastqRs <- sort(list.files(fastq_dir, pattern = ".*R2.*_primer_cut.fastq.gz"))

## Already run - need to figure how to code to not run if files are present
if (length(fastqFs) != length(fastqRs)) {
    stop("Forward and reverse files do not match.")
}

# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(
    fwd = file.path(fastq_dir, fastqFs),
    filt = file.path(filtpathF, fastqFs),
    rev = file.path(fastq_dir, fastqRs),
    filt.rev = file.path(filtpathR, fastqRs),
    truncLen = c(240, 200),
    maxEE = 2,
    truncQ = 11,
    maxN = 0,
    rm.phix = TRUE,
    matchIDs = TRUE,
    compress = TRUE,
    verbose = TRUE,
    multithread = TRUE
)

## Run inference on individual seq runs ----------------------------------------

## Inferr seq variants
infer_seq_vars <- function(seq_run){
    seqtab_file <- paste0("dada2/", seq_run,"_seqtab.rds")
    ## only run if seqtab file does not exist
    if (file.exists(seqtab_file)) return(readRDS(seqtab_file))
        
    filtFs <- list.files("dada2/filtered_F",
                         pattern = paste0(seq_run, ".*fastq.gz"),
                         full.names = TRUE)
    
    filtRs <- list.files("dada2/filtered_R",
                         pattern = paste0(seq_run, ".*fastq.gz"),
                         full.names = TRUE)
    
    sample.names <- sapply(strsplit(basename(filtFs), "_S"), `[`, 1)
    
    sample.namesR <- sapply(strsplit(basename(filtRs), "_S"), `[`, 1)
    
    if (!identical(sample.names, sample.namesR))
        stop("Forward and reverse files do not match.")
    
    names(filtFs) <- sample.names
    names(filtRs) <- sample.names
    set.seed(100)
    
    # Learn forward error rates
    errF <- learnErrors(filtFs, nread = 1e6, multithread = TRUE)
    
    # Learn reverse error rates
    errR <- learnErrors(filtRs, nread = 1e6, multithread = TRUE)
    
    # Sample inference and merger of paired-end reads
    mergers <- vector("list", length(sample.names))
    names(mergers) <- sample.names
    for (sam in sample.names) {
        cat("Processing:", sam, "\n")
        derepF <- derepFastq(filtFs[[sam]])
        ddF <- dada(derepF, err = errF, multithread = TRUE)
        derepR <- derepFastq(filtRs[[sam]])
        ddR <- dada(derepR, err = errR, multithread = TRUE)
        merger <- mergePairs(ddF, derepF, ddR, derepR)
        mergers[[sam]] <- merger
    }
    
    seqtab <- makeSequenceTable(mergers)
    
    ## Write to disk
    saveRDS(seqtab, seqtab_file)
    
    seqtab
}






## Sequence tab files ----------------------------------------------------------

combined_seqtab_file <- "dada2/combined_seqtab.rds"
if (!file.exists(combined_seqtab_file)) {
    run_seqtabs <-
        c("nist_run1", "nist_run2", "jhu_run1", "jhu_run2") %>%
        set_names(.) %>%
        map(infer_seq_vars)
    
    combined_seqtab <- mergeSequenceTables(
        run_seqtabs$nist_run1,
        run_seqtabs$nist_run2,
        run_seqtabs$jhu_run1,
        run_seqtabs$jhu_run2
    )
    
    saveRDS(combined_seqtab, combined_seqtab_file)
} else{
    combined_seqtab <- readRDS(combined_seqtab_file)
}


## Remove Chimeras -------------------------------------------------------------
seqtab_nochim_file <- "dada2/seqtab_nochim.rds"
if (!file.exists(seqtab_nochim_file)) {
    seqtab_nochim <- removeBimeraDenovo(combined_seqtab, verbose = TRUE)
    saveRDS(seqtab_nochim, seqtab_nochim_file)
}else{
    seqtab_nochim <- readRDS(seqtab_nochim_file)
}

## Assign Taxonomy -------------------------------------------------------------
taxa_file <- "dada2/taxa.rds"
ref_url <- "http://dx.doi.org/10.5281/zenodo.158958"
ref_file <- "reference/silva_nr_v123_train_set.fa.gz"
if (!file.exists(taxa_file)) {
    ## Download reference db if not present
    if (!file.exists(ref_file)) download.file(ref_url, ref_file)
    
    taxa <- assignTaxonomy(seqtab_nochim, ref_file)
    colnames(taxa) <- paste0("Rank",1:6)
    saveRDS(taxa, taxa_file)
}else{
    taxa <- readRDS(taxa_file)
}

## Representative Sequences ----------------------------------------------------
seq_file <- "dada2/sv_seqs.fasta"
if (!file.exists(seq_file)) {
    ## Extract sequences and name features
    sv_seqs <- getSequences(seqtab_nochim)
    names(sv_seqs) <- paste0("SV",1:length(sv_seqs))
    
    ## Write sequences to file
    sv_seqs <- DNAStringSet(sv_seqs)
    writeXStringSet(sv_seqs, seq_file)
    
} else {
    sv_seqs <- readDNAStringSet(seq_file)
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
    fitGTR <- update(fit, k = 4, inv = 0.2)
    fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                        rearrangement = "stochastic", 
                        control = pml.control(trace = 0))
    
    saveRDS(fitGTR, tree_file)
}


## generate mgtst.files for mothur pipeline
library(mgtst)
library(dplyr)
library(stringr)
library(readr)

project_dir <- "~/Projects/16S_etec_mix_study/"
# fastq_dir <- paste0(project_dir, "analysis/pipelines/mothur/data/raw")
# fastq_dir <- paste0(project_dir, "data/raw_seq/nist_run1/170209_16s_MiSeq_Run_Folder")
fastq_dir <- paste0(project_dir, "data/raw_seq/nist_run2")
fastqs <- list.files(fastq_dir,recursive = TRUE, pattern = "NIST.*001.fastq", full.names = TRUE)
# fastqs <- fns[grepl("*001.fastq$", fns)]
fastqs <- sort(fastqs) # Sort should keep them paired in order
fnFs <- fastqs[grepl("_R1", fastqs)]
fnRs <- fastqs[grepl("_R2", fastqs)]

## Run
grp_df <- data_frame(read1 = fnFs, read2 = fnRs) %>%
    mutate(ill_id = basename(fnFs), ill_id = str_replace(ill_id, "_.*",""))


## study metadata
data(sample_sheet)
meta_df <- sample_sheet %>%
    mutate(pos_ns = str_replace(pos, "_",""),
           ill_id = paste("NIST", pcr_16S_plate, pos_ns, sep = "-")) %>%
    filter(seq_lab == "NIST", barcode_lab == "NIST") %>%
    mutate(pcr_16S_plate = as.character(pcr_16S_plate)) %>%
    select(id, ill_id) %>% left_join(grp_df) %>% select(-id)
# write_tsv(meta_df,"mgtst.files",col_names = FALSE)
# write_tsv(meta_df,"mgtstNIST1.files",col_names = FALSE)
write_tsv(meta_df,"mgtstNIST2.files",col_names = FALSE)

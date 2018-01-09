


library(dplyr, lib="~/R/x86_64-pc-linux-gnu-library/")
library(stringr)
library(readr)

## Get fastq files
fastqs <- list.files("data", recursive = TRUE, pattern = "*001.fastq.gz", 
                     full.names = TRUE) %>% sort()
fnFs <- fastqs[grepl("_R1", fastqs)]
fnRs <- fastqs[grepl("_R2", fastqs)]

## Generate DF and extract ids
grp_df <- data_frame(read1 = fnFs, read2 = fnRs) %>%
    mutate(id = basename(fnFs), id = str_replace(id, "_S.*","")) %>% 
    mutate(id = str_replace(id, "-","_"))

grp_df %>% select(id, read1, read2) %>% 
      write_tsv("mothur/mgtst.files",col_names = FALSE)

## generate mgtst.files for mothur pipeline, run from the mgtst_pipelines directory
library(tidyverse)
library(stringr)


## Get fastq files
get_fastqs <- function(fq_dir){
    fastqs <- list.files(fq_dir, recursive = TRUE, #pattern = "*001.fastq.gz", 
                         full.names = TRUE) %>% sort()
    fnFs <- fastqs[grepl("_R1", fastqs)]
    fnRs <- fastqs[grepl("_R2", fastqs)]
    
    ## Generate DF and extract ids
    grp_df <- data_frame(read1 = fnFs, read2 = fnRs) %>%
        mutate(id = basename(fnFs), id = str_replace(id, "_.*","")) %>% 
        mutate(id = str_replace(id, "-","_"))  
    
    ## Output data frame
    grp_df %>% select(id, read1, read2)
}

fq_df <- list(jhu_1  = "../Projects/16S_etec_mix_study/data/fastq/jhu_run1",
              jhu_2  = "../Projects/16S_etec_mix_study/data/fastq/jhu_run2",
              nist_1 = "../Projects/16S_etec_mix_study/data/fastq/nist_run1",
              nist_2 = "../Projects/16S_etec_mix_study/data/fastq/nist_run2") %>% 
    map_df(get_fastqs, .id = "run") %>% 
    mutate(id = paste(run, id, sep = "_")) %>% 
    select(-run)

fq_df %>% write_tsv("mothur/mgtst.files", col_names = FALSE)

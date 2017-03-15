#!/usr/bin/sh
# QIIME Open Reference Pipeline 
## Based on the Illumina Tutorial http://qiime.org/1.6.0/tutorials/illumina_overview_tutorial.html

## Merge forward and reverse reads
multiple_join_paired_ends.py -i ../data/ -o joined_pairs
rm joined_pairs/*/*un*
multiple_split_libraries_fastq.py -i joined_pairs -o split_libs --read_indicator fastqjoin.join \
	--demultiplexing_method sampleid_by_file --include_input_dir_path --remove_filepath_in_name

# Clustering
## Open Reference clustering
pick_open_reference_otus.py -o otus_uc_fast/ -i split_libs/seqs.fna -p ../code/qiime_params.txt


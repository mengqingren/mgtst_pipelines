#!/usr/bin/sh
# QIIME Open Reference Pipeline 
## Based on the Illumina Tutorial http://qiime.org/1.6.0/tutorials/illumina_overview_tutorial.html

#multiple_join_paired_ends.py -i ../data/ -o joined_pairs
#rm joined_pairs/*/*un*

#removed unpaired reads - 
#multiple_split_libraries_fastq.py -i joined_pairs -o split_libs --read_indicator fastqjoin.join --demultiplexing_method sampleid_by_file --include_input_dir_path --remove_filepath_in_name

## Chimera checking
#identify_chimeric_seqs.py -i split_libs/seqs.fna -o split_libs_chimera -m usearch61 --threads 6 --suppress_usearch61_ref

# Clustering
## Open Reference clustering
pick_open_reference_otus.py -o otus_uc_fast_no_chimera/ -i split_libs_chimera/seqs.fna_consensus_fixed.fasta -p ../code/qiime_params.txt


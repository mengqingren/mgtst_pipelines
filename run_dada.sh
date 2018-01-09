#!/usr/bin/sh
#SBATCH -J dada_full
#SBATCH --qos=large \
#SBATCH --mem=128gb \
#SBATCH --time=72:00:00 \
#SBATCH --cpus-per-task=32

module load Python2/common/2.7.12
module load R/3.4.0
R CMD BATCH code/dada2_big_data_pipeline.R

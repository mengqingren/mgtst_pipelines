#!/usr/bin/sh
#SBATCH -J qiime2_deblur
#SBATCH --qos=large \
#SBATCH --mem=128gb \
#SBATCH --time=36:00:00 \
#SBATCH --cpus-per-task=32

source ~/add_conda3
source activate qiime2-2017.12

cd qiime/

deblur workflow --seqs-fp split_libs/seqs.fna --output-dir deblur_output -t 430

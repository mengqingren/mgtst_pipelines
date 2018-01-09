#!/usr/bin/sh
#SBATCH -J qiime_open
#SBATCH --qos=large \
#SBATCH --mem=128gb \
#SBATCH --time=36:00:00 \
#SBATCH --cpus-per-task=32

module load qiime
module load ea_utils
cd qiime
bash ../code/qiime_pipeline.sh

#!/usr/bin/sh
#SBATCH -J mothur_pipe
#SBATCH --qos=large \
#SBATCH --mem=128g \
#SBATCH --time=72:00:00 \
#SBATCH --cpus-per-task=32

module load R/3.4.0
make -f code/mothur/Makefile
#make -f code/mothur/Makefile clean

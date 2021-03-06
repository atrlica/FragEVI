#!/bin/sh -l

#$ -l h_rt=24:00:00
#$ -N frag.proc
#$ -V
#$ -j y
#$ -l mem_total=125G
#$ -pe omp 8 -q "geo*"
#$ -v OMP_NUM_THREADS=8
#$ -m e
#$ -M atrlica@bu.edu
#$ -o /projectnb/buultra/atrlica/FragEVI/logs

# cd /projectnb/buultra/atrlica/BosAlbedo/working/logs

module load R_earth/3.1.0
Rscript /projectnb/buultra/atrlica/FragEVI/Rscripts/frag.processing.R

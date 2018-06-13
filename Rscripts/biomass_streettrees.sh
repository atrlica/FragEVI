#!/bin/sh -l

#$ -l h_rt=18:00:00
#$ -N biom.proc
#$ -V
#$ -j y
#$ -l mem_total=125G
#$ -pe omp 4 -q "geo*"
#$ -v OMP_NUM_THREADS=4
#$ -m e
#$ -M atrlica@bu.edu
#$ -o /projectnb/buultra/atrlica/FragEVI/logs

module load R_earth/3.1.0
Rscript /projectnb/buultra/atrlica/FragEVI/Rscripts/biomass_streettrees.R

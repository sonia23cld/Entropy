#!/bin/bash
#PBS -N Plot_estimates
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=8:mem=50gb


#trap 'clean_scratch' TERM EXIT

#if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi


#cd $SCRATCHDIR

module load conda-modules-py37 || exit 1

conda activate /storage/brno3-cerit/home/sonia_celestini/.conda/envs/r_entropy

#module add R-4.0.0-gcc
#cp /storage/brno3-cerit/home/sonia_celestini/Entropy/Plot_estimates.R . 

Rscript --vanilla /storage/brno3-cerit/home/sonia_celestini/CalSil_2024/R_scripts/plotadmix.R
echo "rplots finished at `date`"

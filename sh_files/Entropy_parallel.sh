#!/bin/bash -e
#PBS -N Entropy
#PBS -l walltime=500:00:00
#PBS -l select=1:ncpus=4:mem=100gb:scratch_local=40gb
#PBS -j oe

module add entropy-2.0

trap 'clean_scratch' TERM EXIT
if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi

DATADIR="/storage/brno12-cerit/home/sonia_celestini/CalSil_2024/Entropy"
qkfile="qk${k}inds.txt"
output="mcmcoutk${k}chain2_60k20k.hdf5"
 
cp $DATADIR/Arenosa_snp_raw.snps.fourfold.dp8nc.m0.5_Alldata76Pops.rep1.VCF_Pruned.mpgl $SCRATCHDIR || exit 1
cp $DATADIR/ploidy_inds.txt $SCRATCHDIR || exit 1
cp $DATADIR/$qkfile $SCRATCHDIR || exit 1

cd $SCRATCHDIR || exit 2
echo data loaded at `date`

entropy -i Arenosa_snp_raw.snps.fourfold.dp8nc.m0.5_Alldata76Pops.rep1.VCF_Pruned.mpgl -m 1 -n ploidy_inds.txt -k $k -q $qkfile -l 60000 -b 20000 -t 20 -D 1 -o $output -w 1 

ls


rm *mpgl
rm ploidy_inds.txt
rm qk*inds.txt
cp $SCRATCHDIR/* $DATADIR/ || export CLEAN_SCRATCH=false
printf "\nFinished\n\n"

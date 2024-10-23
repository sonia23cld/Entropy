Pipeline I used to run Entropy (from a cluster) with mixed-ploidy populations data - Shastry et al. 2021

Shastry V, Adams PE, Lindtke D, Mandeville EG, Parchman TL, Gompert Z, Alex Buerkle C. 2021. Model-based genotype and ancestry estimation for potential hybrids with mixed-ploidy. Molecular Ecology Resources [Internet] 21:1434–1451. Available from: http://dx.doi.org/10.1101/2020.07.31.231514!(https://github.com/user-attachments/assets/33327d01-3bc8-42c6-9f27-1926ca6a4759)

1. Create input data (mpgl file) from pruned vcf file using the R script inputdataformat_Entropy.R
2. Check that the data is correct:

'wc -l Arenosa_snp_raw.snps.fourfold.dp8nc.m0.5_Alldata76Pops.rep1.VCF_Pruned.mpgl' 
This number should be one greater than the number of loci in your data set

'awk -F' ' 'NR==3{print NF}' Arenosa_snp_raw.snps.fourfold.dp8nc.m0.5_Alldata76Pops.rep1.VCF_Pruned.mpgl'
This number should be one greater than the number of genotype classes times the number of individuals in your data set. This should be true for all lines in your input data file after line number 3.
I have: 1958 --> 260x5(tetra) + 219x3(dip)= 1957+1=1958

3. Run Entropy. I have written a script (Entropy_parallel.sh) to run entropy in parallel for several k. You just have to run it several times to have replicates (change chain name).
4. Plot admixture proportions. First, create the input files with plot_admix.sh, which runs plot_admxture.R. Then, download locally the results and analyse them with R script Entropy_plotadmix.R. It gives map with piecharts of admixture proportions and barplot admixture.

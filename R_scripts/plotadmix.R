## R script to plot admixture proportions from inds
## and output text files with Q-matrix coefficients

library(rhdf5)
library(coda)
set.seed(10041996)

setwd("/storage/brno12-cerit/home/sonia_celestini/CalSil_2024/Entropy/")
#colors
mycols<-c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#3288BD","blue","#5E4FA2", "violet")

# readings the ind names from previously stored text file
inds.z<-read.table("cinds_Arenosa_snp_raw.snps.fourfold.dp8nc.m0.5_Alldata76Pops.rep1.VCF_Pruned.txt",sep="\n")
# creating a population label based on ind label
locs.z<-substr(inds.z[,1],1,3)

# writing the population labels for each ind into a text file
write.table(as.integer(as.factor(locs.z)), file="indlabel.txt", row.names=F, 
            col.names=F, quote=F)

# writing the population labels for plotting
write.table(cbind(1:length(unique(locs.z)), unique(locs.z)), file="poplabel.txt",
            row.names=F, col.names=F, quote=F)

# creating a pdf of admxiture plots for all K - change the k based on what you used
#pdf("admix_arenosa.pdf", width=15, height=10)
par(mfrow=c(6,1))
for(k in 3:6){
  if(k!=10) {
    par(mar=c(1.2,2,0.7,1))}
  else {
    par(mar=c(2.4,2,0.7,1))}
  
  # averaging over estimates from each chain
  n=0
  files<- list()
  for (chain in 1:3) {
  if (file.exists(paste0('mcmcoutk', k, 'chain', chain, '_60k20k.hdf5'))) {
  if (file.info(paste0('mcmcoutk', k, 'chain', chain, '_60k20k.hdf5'))$size > 0) {
    qch<-h5read(paste0('mcmcoutk', k, 'chain', chain, '_60k20k.hdf5'),"q")
    files[[n+1]]<-apply(qch, c(2,3), mean)
    n=n+1
     }
    }
  }

  if (n == 1) {
    qest<- files[[1]]
  }
  if (n == 2) {
    qest<- (files[[1]]+files[[2]])/2
  }
  if (n == 3) {
    qest<- (files[[1]]+files[[2]]+files[[3]])/3
  }
  
  # creating a simple Q-matrix file for use with CLUMPAK, pong, etc
  write.table(format(t(qest), digits=5), file=paste0("simpleQ",k,".txt"), 
              row.names=F, col.names=F, quote=F)
  
  #barplot(qest, ylim=0:1, col=mycols[1:k], main=paste0("K=",k), border=NA, 
         # space=0, xaxt="n", yaxt="n")
}
dev.off()


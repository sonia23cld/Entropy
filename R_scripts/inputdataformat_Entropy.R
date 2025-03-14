## Script to extract relevant information from vcf files (PL field)
## and prepare input data for entropy in the right format 
# March 2019 -- Vivaswat Shastry

## USAGE: Rscript inputdataformat.R file1.vcf

library(vcfR)
library(MASS)

setwd('/storage/brno3-cerit/home/sonia_celestini/CalSil_2024/Entropy/VCF_Pruned/')
## flag to check if we have mixed ploidy 
mult.files<-F

#args<-commandArgs(TRUE)

vcffile<-as.character('Arenosa_snp_raw.snps.fourfold.dp8nc.m0.5_Alldata76Pops.rep1.VCF_Pruned.vcf.gz')	

#if(!file.exists("ploidy_inds.txt")){
#  cat("Error: ploidy_inds.txt file not found! \nInclude a file titled ploidy_inds.txt with ploidy values for each ind in a new line. \n")
#}

#if(length(args)>2){
#  mult.files<-T
#  ploidy1<-as.character(args[2])
#  vcffile2<-as.character(args[3])
#  ploidy2<-as.character(args[4])
#}

## vcf to mpgl portion

print("Creating input genotype likelihood data for entropy...")

fname<-strsplit(vcffile, ".vcf")[[1]][1]
a.vcf<-read.vcfR(vcffile)

a.ids<-colnames(a.vcf@gt)[-1]
#set new working directory for outputs
setwd('/storage/brno3-cerit/home/sonia_celestini/CalSil_2024/Entropy/')
write.table(a.ids,file=paste0("cinds_",fname,".txt"),row.names=F,col.names=F)

a.loci<-paste0(a.vcf@fix[,'CHROM'],';',a.vcf@fix[,'POS'])

a.pl<-extract.gt(a.vcf, element='PL') #I have NAs

## VS: maybe output some statistics from the vcf file...?
print(paste0(round(sum(is.na(a.pl))*100/(nrow(a.pl)*ncol(a.pl)),1),"% missing data in file"))
nonNAidx<-which(apply(a.pl,1,function(x){sum(is.na(x))})==0)[1] #this is checking which rows don't contain NAs and taking their index. Then it take only one row of them to calculate ploidy
#
if(is.na(nonNAidx)) {
  #	print("There are no loci with data for all individuals in the vcf file.")
  #	print("This is OK but the ploidy file will contain a few zeros which will need to be changed.")
  ploidy<-unlist(lapply(a.pl[10,],function(x){length(gregexpr(",",x,fixed=T)[[1]])}))
} else {
  ploidy<-unlist(lapply(a.pl[nonNAidx,],function(x){length(gregexpr(",",x,fixed=T)[[1]])}))
}

if(mult.files){
  fname2<-strsplit(vcffile2, ".vcf")
  a.vcf2<-read.vcfR(vcffile2)
  
  a.ids2<-colnames(a.vcf2@gt)[-1]
  write.table(a.ids2,file=paste0("inds_",fname2,".txt"),row.names=F,col.names=F)
  
  a.pl2<-extract.gt(a.vcf2, element='PL')
  a.pl2[is.na(a.pl2)]<-paste(rep("0",ploidy2+1),collapse=",")
  
  a.loci2<-paste0(a.vcf2@fix[,'CHROM'],';',a.vcf2@fix[,'POS'])
  
  nind<-dim(a.pl)[2]+dim(a.pl2)[2]
  commonloci<-a.loci %in% a.loci2
  nloci<-dim(a.pl)[1]+dim(a.pl2)[1]-length(commonloci.idx)
  
  total.pl<-matrix("0", nrow=nloci, ncol=nind)
  total.pl[1:length(commonloci),1:dim(a.pl)[2]]<-a.pl[commonloci,]
  total.pl[(length(commonloci)+1):dim(a.pl)[1],1:dim(a.pl)[2]]<-a.pl[!commonloci,]
  
  commonloci<-a.loci2 %in% a.loci
  total.pl[1:length(commonloci),(dim(a.pl)[2]+1):nind]<-a.pl2[commonloci,]
  total.pl[(dim(a.pl)[1]+1):nloci,(dim(a.pl)[2]+1):nind]<-a.pl2[!commonloci,]
  
  total.pl[(length(commonloci)+1):dim(a.pl)[1],(dim(a.pl)[2]+1):nind]<-paste(rep("0",ploidy+1),collapse=",")
  total.pl[(dim(a.pl)[1]+1):nloci,1:dim(a.pl)[2]]<-paste(rep("0",ploidy2+1),collapse=",")
  
  total.pl<-apply(total.pl, c(1,2), function(x){as.numeric(unlist(strsplit(x,split=",")))})
  
  print(paste0("Number of inds: ", nind, "    Number of loci: ", nloci)) 
  
  loci.names<-a.loci[a.loci %in% a.loci2]
  loci.names<-append(a.loci[!(a.loci %in% a.loci2)],loci.names)
  loci.names<-append(a.loci2[!(a.loci2 %in% a.loci)],loci.names)
  
  unlink("combined.mpgl")
  cat(paste(ncol(total.pl), nrow(total.pl), "\n"), file="combined.mpgl")
  cat(paste(a.ids), file="combined.mpgl", append=T)
  for(locus in 1:nloci){
    cat(paste0("\n",loci.names[locus]," "), file="combined.mpgl", append=T)
    write.table(unlist(total.pl[locus,]), file="combined.mpgl", append=T,
                col.names=F, row.names=F, eol=" ")
  }
} else {
  a.gl<-apply(a.pl, c(1,2), function(x){as.numeric(unlist(strsplit(x,split=",")))})
  
  #ploidy<-apply(a.gl,2,function(x){length(x[[1]])-1})
  #ploidy<-read.table("ploidy_inds.txt",header=F,stringsAsFactors=F)[,1]
  
  #if(length(unique(ploidy))==1){
  #	ploidy<-ploidy[1]
  #} else {
  #	print("Writing ploidy value for each individual into file...")
  #	print("[Check to make sure there are no zeros in this file]")
  write.table(ploidy, file="ploidy_inds.txt", row.names=F, col.names=F, sep="\n")
  #}
  
  nind<-dim(a.pl)[2]
  nloci<-dim(a.pl)[1]
  
  print(paste0("Number of inds: ", nind, "    Number of loci: ", nloci)) 
  
  for(i in 1:nrow(a.pl)){
    for(j in 1:ncol(a.pl)){
      if(is.na(a.gl[i,j])){
        a.gl[i,j][[1]]<-rep(0, ploidy[[j]]+1) #there was a bug here. It was ploidy[[1]]
      }
    }
  }
  unlink(paste0(fname, ".mpgl"))
  cat(paste(nind, nloci, "\n"), file=paste0(fname, ".mpgl"))
  cat(paste(a.ids), file=paste0(fname, ".mpgl"), append=T)
  for(locus in 1:nloci){
    cat(paste0("\n",a.vcf@fix[locus,'CHROM']," "), file=paste0(fname, ".mpgl"), append=T)
    write.table(unlist(t(a.gl[locus,])), file=paste0(fname, ".mpgl"), append=T,
                col.names=F, row.names=F, eol=" ")
  }
}


#--------------------------------------------------------------------

## mpgl to point estimates 

print("Converting genotype likelihoods into point estimates...")

if(mult.files){
  total.mean.gl<-matrix(NA, nrow=nrow(total.pl), ncol=ncol(total.pl))
  for(i in 1:nrow(total.pl)){
    for(j in 1:ncol(total.pl)){
      if(length(unique(unlist(total.pl[i,j])))!=1){
        temp<-10^(total.pl[i,j][[1]]/-10)/(sum(10^(total.pl[i,j][[1]]/-10)))
        if(length(temp)==2){
          total.mean.gl[i,j]<-sum(temp*(0:ploidy[j]))
        }
        else{
          total.mean.gl[i,j]<-sum(temp*(0:ploidy[j]))
        }
      }
    }
  }
  a.mean.gl<-total.mean.gl
}else{
  #a.mean.gl<-matrix(0, nrow=nloci, ncol=nind)
  #a.mean.gl<-apply(a.gl, c(1,2), function(x){temp<-sum(10^(-0.1*unlist(x)));sum((10^(-0.1*unlist(x)))*(0:max(ploidy)))/temp}) #here there is bug, you can not apply max(ploidy) also for 2x
  a.mean.gl2<-matrix(0, nrow=nloci, ncol=nind)
  for (i in 1:ncol(a.gl)) {
    a.mean.gl2[,i]<- sapply(a.gl[,i], function(x){temp<-sum(10^(-0.1*unlist(x)));sum((10^(-0.1*unlist(x)))*(0:ploidy[[i]]))/temp}) 
  }
  colnames(a.mean.gl2)<- colnames(a.gl)
  rownames(a.mean.gl2)<- rownames(a.gl)
}
#write.table(a.mean.gl, file="pntest_mean_gl_ProbablyWrong.txt", row.names=F, col.names=F, quote=F, sep=" ")
write.table(a.mean.gl2, file="pntest_mean_gl.txt", row.names=F, col.names=F, quote=F, sep=" ")

#--------------------------------------------------------------------

# point estimates to ancestry ldak estimates

print("Obtaining admixture proportions for chain initialization...")

do.pca<-function(gmat){
  print("1. centering genotype matrix")
  gmn<-apply(gmat,1,mean, na.rm=T)
  gmnmat<-matrix(gmn,nrow=nrow(gmat),ncol=ncol(gmat))
  gprime<-gmat-gmnmat ## remove mean
  
  print("2. finding pairwise covariance")
  gcovarmat<-matrix(NA,nrow=ncol(gmat),ncol=ncol(gmat))
  for(i in 1:ncol(gmat)){
    for(j in i:ncol(gmat)){
      if (i==j){
        gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
      }
      else{
        gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
        gcovarmat[j,i]<-gcovarmat[i,j]
      }
    }
  }
  missing.inds<-unique(which(is.na(gcovarmat),arr.ind=T)[,2])
  print("3. final principal comp analysis")
  return(list(prcomp(x=na.omit(gcovarmat),center=TRUE,scale=FALSE), missing.inds))
}
return.val<-do.pca(a.mean.gl2)
pcout<-return.val[1][[1]]
miss.inds<-return.val[2][[1]]
pcSummary<-summary(pcout)

print("4. k-means, lda and writing initial admixture proportions to file")

for(k in 3:6){
  init.admix<-matrix(0, nrow=nind, ncol=k)
  init.admix[miss.inds,]<-rep(1/k, k)
  
  kn<-kmeans(pcout$x[,1:5], k, iter.max=10, nstart=10, algorithm="Hartigan-Wong")
  ldakn<-lda(x=pcout$x[,1:5], grouping=kn$cluster, CV=T)
  init.admix<-ldakn$posterior
  write.table(round(init.admix, 5), quote=F, row.names=F, col.names=F, file=paste0("qk",k,"inds.txt"))
}

#--------------------------------------------------------------------

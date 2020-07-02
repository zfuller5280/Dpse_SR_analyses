#install.packages("poolfstat")
library("poolfstat")
library(boot)


poolfstat::computeFST

computeFST_windowed=function(pooldata,method="Anova",snp.index=1:nsnp,nperm=1000){
  #print(snp.index)
  npop=pooldata@npools 
  YY=pooldata@refallele.readcount[snp.index,]
  NN=pooldata@readcoverage[snp.index,]
  poolsize=pooldata@poolsizes
  nsnp = length(snp.index)
  #print(YY)
  #print(NN)
  YY_real = YY
  NN_real = NN
  print(poolsize)
  if(method=="Identity"){
    Q1=as.matrix((YY*(YY-1) + (NN-YY)*(NN-YY-1) )/(NN*(NN-1)))
    Q1 = (1/(matrix(1,nsnp,npop) %*% diag(poolsize-1)))*(Q1 %*% diag(poolsize) - 1)
    lambdaj=poolsize*(poolsize-1)
    lambdaj=lambdaj/sum(lambdaj)
    snp.Q1=rowSums(Q1%*%diag(lambdaj))
    hat.Q1=mean(snp.Q1,na.rm=T) 
    
    Q2=matrix(0,nrow(YY),npop*(npop-1)/2)
    omegajj=rep(0,npop*(npop-1)/2)
    cnt=0
    for(i in 1:(npop-1)){
      for(j in (i+1):npop){
        cnt=cnt+1
        omegajj[cnt]=poolsize[i]*poolsize[j]
        Q2[,cnt]=(YY[,i]*YY[,j] + (NN[,i]-YY[,i])*(NN[,j]-YY[,j]))/(NN[,i]*NN[,j])
      }
    }
    snp.Q2=rowSums(Q2%*%diag(omegajj/sum(omegajj)))
    hat.Q2=mean(snp.Q2,na.rm=TRUE)
    rslt_FST=(hat.Q1-hat.Q2)/(1-hat.Q2)
    snp.Q1=snp.Q1
    snp.Q2=snp.Q2
    snp.FST=(snp.Q1-snp.Q2)/(1-snp.Q2)
  }
  
  if (method=="Anova"){
    mtrx.n_i <- matrix(poolsize,nrow = nsnp,ncol = npop,byrow = TRUE)
    #print(mtrx.n_i)
    C_1 <- rowSums(NN)
    C_2 <- rowSums(NN^2)
    D_2 <- rowSums(NN / mtrx.n_i + (mtrx.n_i - 1) / mtrx.n_i,na.rm = TRUE)
    D_2.star <- rowSums(NN * (NN / mtrx.n_i + (mtrx.n_i - 1) / mtrx.n_i),na.rm = TRUE) / C_1
    n_c <- (C_1 - C_2 / C_1) / (D_2 - D_2.star)
    YY_alt <- NN - YY
    SSI <- rowSums(YY - YY^2 / NN + YY_alt - YY_alt^2 / NN,na.rm = TRUE)
    SSP <- rowSums(NN * ((YY / NN) - (rowSums(YY) / C_1))^2 + NN * ((YY_alt / NN) - (rowSums(YY_alt) / C_1))^2,na.rm = TRUE)
    MSI <- SSI / (C_1 - D_2)
    MSP <- SSP / (D_2 - D_2.star)
    F_ST <- (MSP - MSI)  / (MSP + (n_c - 1) * MSI)
    keep <- !is.na(F_ST)
    F_ST_multi <- sum(MSP[keep] - MSI[keep])  / sum(MSP[keep] + (n_c[keep] - 1) * MSI[keep])
    Q_1 <- 1 - MSI
    Q_2 <- 1 - MSI - (MSP - MSI) / n_c
    snp.FST = F_ST
    snp.Q1 = Q_1
    snp.Q2 = Q_2
    rslt_FST = F_ST_multi
  }
  
  bootstrapped_fst<-vector()
  for(x in seq(1:nperm)){
    perm<-sample(nrow(YY), nrow(YY),replace=T)
    YY = YY_real[perm,]
    NN = NN_real[perm,]
    if (method=="Anova"){
      mtrx.n_i <- matrix(poolsize,nrow = nsnp,ncol = npop,byrow = TRUE)
      C_1 <- rowSums(NN)
      C_2 <- rowSums(NN^2)
      D_2 <- rowSums(NN / mtrx.n_i + (mtrx.n_i - 1) / mtrx.n_i,na.rm = TRUE)
      D_2.star <- rowSums(NN * (NN / mtrx.n_i + (mtrx.n_i - 1) / mtrx.n_i),na.rm = TRUE) / C_1
      n_c <- (C_1 - C_2 / C_1) / (D_2 - D_2.star)
      YY_alt <- NN - YY
      SSI <- rowSums(YY - YY^2 / NN + YY_alt - YY_alt^2 / NN,na.rm = TRUE)
      SSP <- rowSums(NN * ((YY / NN) - (rowSums(YY) / C_1))^2 + NN * ((YY_alt / NN) - (rowSums(YY_alt) / C_1))^2,na.rm = TRUE)
      MSI <- SSI / (C_1 - D_2)
      MSP <- SSP / (D_2 - D_2.star)
      F_ST <- (MSP - MSI)  / (MSP + (n_c - 1) * MSI)
      keep <- !is.na(F_ST)
      F_ST_multi <- sum(MSP[keep] - MSI[keep])  / sum(MSP[keep] + (n_c[keep] - 1) * MSI[keep])
      Q_1 <- 1 - MSI
      Q_2 <- 1 - MSI - (MSP - MSI) / n_c
      bootstrapped_fst[x]<-F_ST_multi
    }
    if(method=="Identity"){
      Q1=as.matrix((YY*(YY-1) + (NN-YY)*(NN-YY-1) )/(NN*(NN-1)))
      Q1 = (1/(matrix(1,nsnp,npop) %*% diag(poolsize-1)))*(Q1 %*% diag(poolsize) - 1)
      lambdaj=poolsize*(poolsize-1)
      lambdaj=lambdaj/sum(lambdaj)
      snp.Q1=rowSums(Q1%*%diag(lambdaj))
      hat.Q1=mean(snp.Q1,na.rm=T) 
      
      Q2=matrix(0,nrow(YY),npop*(npop-1)/2)
      omegajj=rep(0,npop*(npop-1)/2)
      cnt=0
      for(i in 1:(npop-1)){
        for(j in (i+1):npop){
          cnt=cnt+1
          omegajj[cnt]=poolsize[i]*poolsize[j]
          Q2[,cnt]=(YY[,i]*YY[,j] + (NN[,i]-YY[,i])*(NN[,j]-YY[,j]))/(NN[,i]*NN[,j])
        }
      }
      snp.Q2=rowSums(Q2%*%diag(omegajj/sum(omegajj)))
      hat.Q2=mean(snp.Q2,na.rm=TRUE)
      bootstrapped_fst[x]<-(hat.Q1-hat.Q2)/(1-hat.Q2)
    }
  }
  
  rslt <- list(snp.FST,snp.Q1,snp.Q2, FST=rslt_FST, bootstrapped_samples=bootstrapped_fst)

  return(rslt)
}
environment(computeFST_windowed) <- asNamespace('poolfstat')
assignInNamespace("computeFST", computeFST_windowed, ns = "poolfstat")

fst_for_scaffold<-function(scaffold, index_1=1, index_2=2){
  infile<-paste(scaffold, ".merged.vcf.gz.indel_filtered.vcf",sep="")
  dat<-vcf2pooldata(vcf.file=infile,poolsize=c(8,8,2,2,1),poolnames=c("DpseST",'DpseSR',"DperST","DperSR","Dmir"))
  subset_dat<-pooldata.subset(dat, pool.index=c(index_1,index_2),min.cov.per.pool=8)
  subset_dat<-pooldata.subset(subset_dat, min.maf=(1/8))
  subset_dat@refallele.readcount<-subset_dat@refallele.readcount[rowSums(is.na(subset_dat@refallele.readcount)) != ncol(subset_dat@refallele.readcount), ]
  subset_dat@readcoverage<-subset_dat@readcoverage[rowSums(is.na(subset_dat@readcoverage)) != ncol(subset_dat@readcoverage), ]
  subset_dat@snp.info<-subset_dat@snp.info[rowSums(is.na(subset_dat@snp.info)) != ncol(subset_dat@snp.info), ]
  total_snps<-nrow(subset_dat@refallele.readcount)
  x<-1
  window_size<-100
  end_x<-window_size
  fst_windows<-matrix(nrow=total_snps/window_size,ncol=9)
  i<-1
  while(end_x<total_snps){
    fst<-computeFST_windowed(subset_dat, snp.index = x:end_x, method="Anova")
    low_bnd<-as.numeric(quantile(fst$bootstrapped_samples,c(0.025, 0.975),names = FALSE, na.rm=T)[1])
    up_bnd<-as.numeric(quantile(fst$bootstrapped_samples,c(0.025, 0.975),names = FALSE, na.rm=T)[2])
    window_start<-as.numeric(subset_dat@snp.info[x:end_x,2][1])
    window_end<-as.numeric(subset_dat@snp.info[x:end_x,2][window_size])
    window_size_bp<-window_end-window_start
    chrom<-subset_dat@snp.info[x:end_x,1][1]
    print(c(chrom, window_start, window_end, window_size_bp, x, end_x, fst$FST, low_bnd, up_bnd))
    fst_windows[i,]<-c(chrom, window_start, window_end, window_size_bp, x, end_x, fst$FST, low_bnd, up_bnd)
    i = i + 1
    end_x = end_x + window_size
    x = x + window_size
  }
  fst_windows<-data.frame(fst_windows)
  return(fst_windows)
}

make_AD_fst<-function(index_1=1, index_2=2){
  scaffolds<-c("XL_group1a","XL_group1e","XL_group3a","XL_group3b","XR_group3a","XR_group5","XR_group6","XR_group8")
  scaffold_fsts<-list()
  count<-1
  for(scaffold in scaffolds){
    scaffold_fsts[[count]]<-fst_for_scaffold(scaffold, index_1, index_2)
    count = count + 1
  }
  fst_raw_df = do.call(rbind, scaffold_fsts)
  fst_raw_df$X2<-as.numeric(as.character(fst_raw_df$X2))
  fst_raw_df$X3<-as.numeric(as.character(fst_raw_df$X3))
  fst_raw_df$X4<-as.numeric(as.character(fst_raw_df$X4))
  fst_raw_df$X7<-as.numeric(as.character(fst_raw_df$X7))
  fst_raw_df$X8<-as.numeric(as.character(fst_raw_df$X8))
  fst_raw_df$X9<-as.numeric(as.character(fst_raw_df$X9))
  head(fst_raw_df)
  
  #Chromosome XR
  sites<-fst_raw_df
  AD1<-subset(sites, sites[,1]=="XL_group3b" & sites[,2] > 1 & sites[,2] < 386183)
  AD2<-subset(sites, sites[,1]=="XL_group3a" & sites[,2] > 1362506 & sites[,2] < 2690836)
  AD3<-subset(sites, sites[,1]=="XL_group1a" & sites[,2] > 2660039 & sites[,2] < 5354993)
  AD4<-subset(sites, sites[,1]=="XR_group6" & sites[,2] > 1 & sites[,2] < 13314419)
  AD5<-subset(sites, sites[,1]=="XR_group8" & sites[,2] > 1 & sites[,2] < 3598432)
  AD6<-subset(sites, sites[,1]=="XR_group8" & sites[,2] > 3598433 & sites[,2] < 9212921)
  AD7<-subset(sites, sites[,1]=="XR_group3a" & sites[,2] > 1 & sites[,2] < 1468910)
  AD8<-subset(sites, sites[,1]=="XR_group5" & sites[,2] > 1 & sites[,2] < 740492)
  
  AD3<-(AD3[rev(rownames(AD3)),])
  AD6<-(AD6[rev(rownames(AD6)),])
  AD7<-(AD7[rev(rownames(AD7)),])
  AD8<-(AD8[rev(rownames(AD8)),])
  head(AD8)
  AD1$"Position"<-AD1[,2]
  AD2$"Position"<-(AD2[,2] - 1362506) + 386183
  AD3$"Position"<-(2694954 - (AD3[,2] - 2660039)) + 386183 + 1328330
  AD4$"Position"<-AD4[,2] + 386183 + 1328330 + 2694954
  AD5$"Position"<-AD5[,2] + 386183 + 1328330 + 2694954 + 13314419
  AD6$"Position"<-(5614488 - (AD6[,2] - 3598433)) + 386183 + 1328330 + 2694954 + 13314419 + 3598432
  AD7$"Position"<-(1468910 - AD7[,2]) + 386183 + 1328330 + 2694954 + 13314419 + 3598432 + 5614488
  AD8$"Position"<-(740492 - AD8[,2]) + 386183 + 1328330 + 2694954 + 13314419 + 3598432 + 5614488 + 1468910
  AD<-rbind(AD1,AD2,AD3,AD4,AD5,AD6,AD7,AD8)
  colnames(AD)<-c("Chrom","Chrom_start","Chrom_end","windowsize","snp_start","snp_end","fst","ci_u","ci_l","w_start")
  AD<-AD[complete.cases(AD), ]
  return(AD)
}
#vcf2pooldata(vcf.file="XR_group5.merged.vcf.gz.indel_filtered.vcf",poolsize=c(8,8,2,2,1))
fst_for_window<-function(scaffold, bp, window_size, index_1=1, index_2=2, min.maf=0){
  infile<-paste(scaffold, ".merged.vcf.gz.indel_filtered.vcf",sep="")
  dat<-vcf2pooldata(vcf.file=infile,poolsize=c(8,8,2,2,1),poolnames=c("DpseST","DpseSR","DperST","DperSR","Dmir"))
  subset_dat<-pooldata.subset(dat, pool.index=c(index_1,index_2),min.cov.per.pool=8)
  subset_dat<-pooldata.subset(subset_dat, min.maf=min.maf)
  subset_dat@refallele.readcount<-subset_dat@refallele.readcount[rowSums(is.na(subset_dat@refallele.readcount)) != ncol(subset_dat@refallele.readcount), ]
  subset_dat@readcoverage<-subset_dat@readcoverage[rowSums(is.na(subset_dat@readcoverage)) != ncol(subset_dat@readcoverage), ]
  subset_dat@snp.info<-subset_dat@snp.info[rowSums(is.na(subset_dat@snp.info)) != ncol(subset_dat@snp.info), ]
  total_snps<-nrow(subset_dat@refallele.readcount)
  dat_bps<-as.numeric(as.character(subset_dat@snp.info[,2]))
  closest_start_idx<-which.min(abs(dat_bps - bp-window_size))
  closest_end_idx<-which.min(abs(dat_bps - bp+window_size))
  start_idx<-min(c(closest_start_idx,closest_end_idx))
  stop_idx<-max(c(closest_start_idx,closest_end_idx))
  #print(c(start_idx,stop_idx))
  #start_idx<-closest_start_idx-1000
  #stop_idx<-closest_start_idx+1000
  #print(subset_dat)
  fst<-computeFST_windowed(subset_dat, snp.index = start_idx:stop_idx, method="Anova")
  #print(fst)
  low_bnd<-as.numeric(quantile(fst$bootstrapped_samples,c(0.025, 0.975),names = FALSE, na.rm=T)[1])
  up_bnd<-as.numeric(quantile(fst$bootstrapped_samples,c(0.025, 0.975),names = FALSE, na.rm=T)[2])
  window_start<-as.numeric(subset_dat@snp.info[start_idx:stop_idx,2][1])
  window_end<-as.numeric(subset_dat@snp.info[start_idx:stop_idx,2][nrow(subset_dat@snp.info[start_idx:stop_idx,])])
  window_size_bp<-window_end-window_start
  chrom<-subset_dat@snp.info[start_idx:stop_idx,1][1]
  print(c(chrom, window_start, window_end, window_size_bp, start_idx, stop_idx, fst$FST, low_bnd, up_bnd))
}

pop_name_convert<-function(population){
  if(population=="DpseST"){
    name = "10700X1"
  }
  if(population=="DpseSR"){
    name = "10700X2"
  }
  if(population=="DperST"){
    name = "11040X1"
  }
  if(population=="DperSR"){
    name = "11040X2"
  }
  if(population=="Dmir"){
    name = "Dmir_SP138"
  }
  return(name)
}

make_AD_dxy<-function(pop1, pop2){
  scaffolds<-c("XL_group1a","XL_group1e","XL_group3a","XL_group3b","XR_group3a","XR_group5","XR_group6","XR_group8")
  
  pop1_id<-pop_name_convert(pop1)
  pop2_id<-pop_name_convert(pop2)
  
  pop_names<-paste(pop1_id, pop2_id, sep="_")
  infile<-paste(pop_names,".cat_dxy",sep='')
  
  dxy_raw_df<-read.table(infile)
  
  #Chromosome XR
  sites<-dxy_raw_df
  AD1<-subset(sites, sites[,1]=="XL_group3b" & sites[,2] > 1 & sites[,2] < 386183)
  AD2<-subset(sites, sites[,1]=="XL_group3a" & sites[,2] > 1362506 & sites[,2] < 2690836)
  AD3<-subset(sites, sites[,1]=="XL_group1a" & sites[,2] > 2660039 & sites[,2] < 5354993)
  AD4<-subset(sites, sites[,1]=="XR_group6" & sites[,2] > 1 & sites[,2] < 13314419)
  AD5<-subset(sites, sites[,1]=="XR_group8" & sites[,2] > 1 & sites[,2] < 3598432)
  AD6<-subset(sites, sites[,1]=="XR_group8" & sites[,2] > 3598433 & sites[,2] < 9212921)
  AD7<-subset(sites, sites[,1]=="XR_group3a" & sites[,2] > 1 & sites[,2] < 1468910)
  AD8<-subset(sites, sites[,1]=="XR_group5" & sites[,2] > 1 & sites[,2] < 740492)
  
  AD3<-(AD3[rev(rownames(AD3)),])
  AD6<-(AD6[rev(rownames(AD6)),])
  AD7<-(AD7[rev(rownames(AD7)),])
  AD8<-(AD8[rev(rownames(AD8)),])
  head(AD8)
  AD1$"Position"<-AD1[,2]
  AD2$"Position"<-(AD2[,2] - 1362506) + 386183
  AD3$"Position"<-(2694954 - (AD3[,2] - 2660039)) + 386183 + 1328330
  AD4$"Position"<-AD4[,2] + 386183 + 1328330 + 2694954
  AD5$"Position"<-AD5[,2] + 386183 + 1328330 + 2694954 + 13314419
  AD6$"Position"<-(5614488 - (AD6[,2] - 3598433)) + 386183 + 1328330 + 2694954 + 13314419 + 3598432
  AD7$"Position"<-(1468910 - AD7[,2]) + 386183 + 1328330 + 2694954 + 13314419 + 3598432 + 5614488
  AD8$"Position"<-(740492 - AD8[,2]) + 386183 + 1328330 + 2694954 + 13314419 + 3598432 + 5614488 + 1468910
  AD<-rbind(AD1,AD2,AD3,AD4,AD5,AD6,AD7,AD8)
  colnames(AD)<-c("Chrom","Chrom_start","Chrom_end","windowsize","missing_sites","dxy","w_start")
  AD<-AD[complete.cases(AD), ]
  return(AD)
}

make_AD_pi<-function(df){
  scaffolds<-c("XL_group1a","XL_group1e","XL_group3a","XL_group3b","XR_group3a","XR_group5","XR_group6","XR_group8")
  
  #Chromosome XR
  sites<-df
  AD1<-subset(sites, sites$chrom=="XL_group3b" & sites$cds_start > 1 & sites$cds_start < 386183)
  AD2<-subset(sites, sites$chrom=="XL_group3a" & sites$cds_start > 1362506 & sites$cds_start < 2690836)
  AD3<-subset(sites, sites$chrom=="XL_group1a" & sites$cds_start > 2660039 & sites$cds_start < 5354993)
  AD4<-subset(sites, sites$chrom=="XR_group6" & sites$cds_start > 1 & sites$cds_start < 13314419)
  AD5<-subset(sites, sites$chrom=="XR_group8" & sites$cds_start > 1 & sites$cds_start < 3598432)
  AD6<-subset(sites, sites$chrom=="XR_group8" & sites$cds_start > 3598433 & sites$cds_start < 9212921)
  AD7<-subset(sites, sites$chrom=="XR_group3a" & sites$cds_start > 1 & sites$cds_start < 1468910)
  AD8<-subset(sites, sites$chrom=="XR_group5" & sites$cds_start > 1 & sites$cds_start < 740492)
  
  AD3<-(AD3[rev(rownames(AD3)),])
  AD6<-(AD6[rev(rownames(AD6)),])
  AD7<-(AD7[rev(rownames(AD7)),])
  AD8<-(AD8[rev(rownames(AD8)),])
  head(AD8)
  AD1$"Position"<-AD1$cds_start
  AD2$"Position"<-(AD2$cds_start - 1362506) + 386183
  AD3$"Position"<-(2694954 - (AD3$cds_start - 2660039)) + 386183 + 1328330
  AD4$"Position"<-AD4$cds_start + 386183 + 1328330 + 2694954
  AD5$"Position"<-AD5$cds_start + 386183 + 1328330 + 2694954 + 13314419
  AD6$"Position"<-(5614488 - (AD6$cds_start - 3598433)) + 386183 + 1328330 + 2694954 + 13314419 + 3598432
  AD7$"Position"<-(1468910 - AD7$cds_start) + 386183 + 1328330 + 2694954 + 13314419 + 3598432 + 5614488
  AD8$"Position"<-(740492 - AD8$cds_start) + 386183 + 1328330 + 2694954 + 13314419 + 3598432 + 5614488 + 1468910
  AD<-rbind(AD1,AD2,AD3,AD4,AD5,AD6,AD7,AD8)
  #colnames(AD)<-c("Chrom","Chrom_start","Chrom_end","windowsize","missing_sites","dxy","w_start")
  AD<-AD[complete.cases(AD), ]
  return(AD)
}


dpseST.Dmir_AD_dxy<-make_AD_dxy("DpseST","Dmir")
head(dpseST.Dmir_AD_dxy)
dpseSR.Dmir_AD_dxy<-make_AD_dxy("DpseSR","Dmir")
head(dpseSR.Dmir_AD_dxy)
dpseST.dpseSR_AD_dxy<-make_AD_dxy("DpseST","DpseSR")
head(dpseST.dpseSR_AD_dxy)

window_size<-125000
scaffold<-"XR_group5"

dpseST.dpseSR_AD_fst<-make_AD_fst(1, 2)
head(dpseST.dpseSR_AD_fst)
dpseST.Dmir_AD_fst<-make_AD_fst(1, 5)
head(dpseST.Dmir_AD_fst)
dpseSR.Dmir_AD_fst<-make_AD_fst(2, 5)
head(dpseSR.Dmir_AD_fst)

head(dpseST.dpseSR_AD_fst)
plot(dpseST.dpseSR_AD_fst$w_start, dpseST.dpseSR_AD_fst$fst)
mean(dpseST.Dmir_AD_fst$fst)
dev.off()
#Make plots of Fst for Figure 2. Read in Fst table with confidence intervals
par(mfrow=c(2,1))
plot(dpseST.dpseSR_AD_fst$w_start, dpseST.dpseSR_AD_fst$fst, col="gray", cex=0.5,ylim=c(0,1.0),xlim=c(0,29200000))
u_ci_avg<-smooth.spline(dpseST.dpseSR_AD_fst$w_start, dpseST.dpseSR_AD_fst$ci_l,spar=0.005)
lines(u_ci_avg,lty=3,lwd=2)
l_ci_avg<-smooth.spline(dpseST.dpseSR_AD_fst$w_start, dpseST.dpseSR_AD_fst$ci_u,spar=0.005)
lines(l_ci_avg,lty=3,lwd=2)
polygon(c(l_ci_avg$x, rev(l_ci_avg$x)), c(u_ci_avg$y, rev(l_ci_avg$y)), col = "grey", border = NA) 
lines(smooth.spline(dpseST.dpseSR_AD_fst$w_start, dpseST.dpseSR_AD_fst$fst,spar=0.005),lwd=2)
abline(v=c(bas_inv_bps[1],bas_inv_bps[2],med_inv_bps[1],med_inv_bps[2],dis_inv_bps[1],dis_inv_bps[2]))

#Inversion breakpoints
bas_inv_bps<-c(6162858, 8796222)
med_inv_bps<-c(9471558,18455166)
dis_inv_bps<-c(25039842,28331788)
25039842-18455166
8790276-6157103
18430099-9465432
28903141-23175372

#Inversion breakpoints (relative to scaffolds)
basal_prox<-c("XR_group6", 1747657)
basal_dist<-c("XR_group6", 4381014)
medial_prox<-c("XR_group6", 5056358)
medial_dist<-c("XR_group8", 706191)
#terminal_prox<-c("XR_group8", 4816529)
terminal_prox<-c("XR_group8", 5499034)
terminal_dist<-c("XR_group5", 415196)

#Subset the regions with 250kb surrounding the breakpoints
dpseST.dpseSR_basal_prox_fst<-fst_for_window(basal_prox[1], as.numeric(basal_prox[2]), (250000/2), min.maf = 1/8)
dpseST.dpseSR_basal_dist_fst<-fst_for_window(basal_dist[1],  as.numeric(basal_dist[2]), (250000/2), min.maf = 1/8)
dpseST.dpseSR_medial_prox_fst<-fst_for_window(medial_prox[1],  as.numeric(medial_prox[2]), (250000/2), min.maf = 1/8)
dpseST.dpseSR_medial_dist_fst<-fst_for_window(medial_dist[1],  as.numeric(medial_dist[2]), (250000/2), min.maf = 1/8)
dpseST.dpseSR_terminal_prox_fst<-fst_for_window(terminal_prox[1],  as.numeric(terminal_prox[2]), (250000/2), min.maf = 1/8)
dpseST.dpseSR_terminal_dist_fst<-fst_for_window(terminal_dist[1],  as.numeric(terminal_dist[2]), (250000/2), min.maf = 1/8)

dpseST.Dmir_basal_prox_fst<-fst_for_window(basal_prox[1], as.numeric(basal_prox[2]), (250000/2), index_1=1, index_2=5, min.maf = 1/8)
dpseST.Dmir_basal_dist_fst<-fst_for_window(basal_dist[1],  as.numeric(basal_dist[2]), (250000/2), index_1=1, index_2=5, min.maf = 1/8)
dpseST.Dmir_medial_prox_fst<-fst_for_window(medial_prox[1],  as.numeric(medial_prox[2]), (250000/2), index_1=1, index_2=5, min.maf = 1/8)
dpseST.Dmir_medial_dist_fst<-fst_for_window(medial_dist[1],  as.numeric(medial_dist[2]), (250000/2), index_1=1, index_2=5, min.maf = 1/8)
dpseST.Dmir_terminal_prox_fst<-fst_for_window(terminal_prox[1],  as.numeric(terminal_prox[2]), (250000/2), index_1=1, index_2=5, min.maf = 1/8)
dpseST.Dmir_terminal_dist_fst<-fst_for_window(terminal_dist[1],  as.numeric(terminal_dist[2]), (250000/2), index_1=1, index_2=5, min.maf = 1/8)

dpseSR.Dmir_basal_prox_fst<-fst_for_window(basal_prox[1], as.numeric(basal_prox[2]), (250000/2), index_1=2, index_2=5)
dpseSR.Dmir_basal_dist_fst<-fst_for_window(basal_dist[1],  as.numeric(basal_dist[2]), (250000/2), index_1=2, index_2=5)
dpseSR.Dmir_medial_prox_fst<-fst_for_window(medial_prox[1],  as.numeric(medial_prox[2]), (250000/2), index_1=2, index_2=5)
dpseSR.Dmir_medial_dist_fst<-fst_for_window(medial_dist[1],  as.numeric(medial_dist[2]), (250000/2), index_1=2, index_2=5)
dpseSR.Dmir_terminal_prox_fst<-fst_for_window(terminal_prox[1],  as.numeric(terminal_prox[2]), (250000/2), index_1=2, index_2=5)
dpseSR.Dmir_terminal_dist_fst<-fst_for_window(terminal_dist[1],  as.numeric(terminal_dist[2]), (250000/2), index_1=2, index_2=5)

dpseST.dpseSR_basal_fst<-mean(c(as.numeric(dpseST.dpseSR_basal_prox_fst[7]), as.numeric(dpseST.dpseSR_basal_dist_fst[7])))
dpseST.dpseSR_medial_fst<-mean(c(as.numeric(dpseST.dpseSR_medial_prox_fst[7]), as.numeric(dpseST.dpseSR_medial_dist_fst[7])))
dpseST.dpseSR_terminal_fst<-mean(c(as.numeric(dpseST.dpseSR_terminal_prox_fst[7]), as.numeric(dpseST.dpseSR_terminal_dist_fst[7])))

dpseST.dpseSR_basal_fst_tbl<-rbind(as.numeric(dpseST.dpseSR_basal_prox_fst), as.numeric(dpseST.dpseSR_basal_dist_fst))
dpseST.dpseSR_medial_fst_tbl<-rbind(as.numeric(dpseST.dpseSR_medial_prox_fst), as.numeric(dpseST.dpseSR_medial_dist_fst))
dpseST.dpseSR_terminal_fst_tbl<-rbind(as.numeric(dpseST.dpseSR_terminal_prox_fst), as.numeric(dpseST.dpseSR_terminal_dist_fst))

dpseST.Dmir_basal_fst<-mean(c(as.numeric(dpseST.Dmir_basal_prox_fst[7]), as.numeric(dpseST.Dmir_basal_dist_fst[7])))
dpseST.Dmir_medial_fst<-mean(c(as.numeric(dpseST.Dmir_medial_prox_fst[7]), as.numeric(dpseST.Dmir_medial_dist_fst[7])))
dpseST.Dmir_terminal_fst<-mean(c(as.numeric(dpseST.Dmir_terminal_prox_fst[7]), as.numeric(dpseST.Dmir_terminal_dist_fst[7])))

dpseSR.Dmir_basal_fst<-mean(c(as.numeric(dpseSR.Dmir_basal_prox_fst[7]), as.numeric(dpseSR.Dmir_basal_dist_fst[7])))
dpseSR.Dmir_medial_fst<-mean(c(as.numeric(dpseSR.Dmir_medial_prox_fst[7]), as.numeric(dpseSR.Dmir_medial_dist_fst[7])))
dpseSR.Dmir_terminal_fst<-mean(c(as.numeric(dpseSR.Dmir_terminal_prox_fst[7]), as.numeric(dpseSR.Dmir_terminal_dist_fst[7])))

dpseST.dpseSR_all_inv_fst<-rbind(dpseST.dpseSR_basal_prox_fst,dpseST.dpseSR_basal_dist_fst,dpseST.dpseSR_medial_prox_fst,dpseST.dpseSR_medial_dist_fst,dpseST.dpseSR_terminal_prox_fst,dpseST.dpseSR_terminal_dist_fst)
dpseST.Dmir_all_inv_fst<-rbind(dpseST.Dmir_basal_prox_fst,dpseST.Dmir_basal_dist_fst,dpseST.Dmir_medial_prox_fst,dpseST.Dmir_medial_dist_fst,dpseST.Dmir_terminal_prox_fst,dpseST.Dmir_terminal_dist_fst)
dpseSR.Dmir_all_inv_fst<-rbind(dpseSR.Dmir_basal_prox_fst,dpseSR.Dmir_basal_dist_fst,dpseSR.Dmir_medial_prox_fst,dpseSR.Dmir_medial_dist_fst,dpseSR.Dmir_terminal_prox_fst,dpseSR.Dmir_terminal_dist_fst)


fact<-mean(c(dpseST.Dmir_basal_fst,dpseST.Dmir_medial_fst,dpseST.Dmir_terminal_fst))

#Get overall Fst for all inversion breakpoints
dpseST.dpseSR_all_inv_fst_mean<-mean(as.numeric(as.character(dpseST.dpseSR_all_inv_fst[,7])))
dpseST.dpseSR_all_inv_fst_cil<-mean(as.numeric(as.character(dpseST.dpseSR_all_inv_fst[,8])))
dpseST.dpseSR_all_inv_fst_ciu<-mean(as.numeric(as.character(dpseST.dpseSR_all_inv_fst[,9])))

dpseST.Dmir_all_inv_fst_mean<-mean(as.numeric(as.character(dpseST.Dmir_all_inv_fst[,7])))
dpseST.Dmir_all_inv_fst_cil<-mean(as.numeric(as.character(dpseST.Dmir_all_inv_fst[,8])))
dpseST.Dmir_all_inv_fst_ciu<-mean(as.numeric(as.character(dpseST.Dmir_all_inv_fst[,9])))

dpseSR.Dmir_all_inv_fst_mean<-mean(as.numeric(as.character(dpseSR.Dmir_all_inv_fst[,7])))
dpseSR.Dmir_all_inv_fst_cil<-mean(as.numeric(as.character(dpseSR.Dmir_all_inv_fst[,8])))
dpseSR.Dmir_all_inv_fst_ciu<-mean(as.numeric(as.character(dpseSR.Dmir_all_inv_fst[,9])))
c(dpseST.dpseSR_all_inv_fst_mean,dpseST.dpseSR_all_inv_fst_cil,dpseST.dpseSR_all_inv_fst_ciu)
c(dpseST.Dmir_all_inv_fst_mean,dpseST.Dmir_all_inv_fst_cil,dpseST.Dmir_all_inv_fst_ciu)
c(dpseST.dpseSR_all_inv_fst_mean,dpseST.dpseSR_all_inv_fst_cil,dpseST.dpseSR_all_inv_fst_ciu)

#Divergence time estimates for all inverted regions
-log(1-dpseST.dpseSR_all_inv_fst_mean) * (2000000/-log(1-dpseST.Dmir_all_inv_fst_mean))/2
-log(1-dpseST.dpseSR_all_inv_fst_cil) * (2000000/-log(1-dpseST.Dmir_all_inv_fst_mean))/2
-log(1-dpseST.dpseSR_all_inv_fst_ciu) * (2000000/-log(1-dpseST.Dmir_all_inv_fst_mean))/2

#Fst transformed estimates of divergence times
-log(1-dpseST.dpseSR_basal_fst) * (2000000/-log(1-dpseST.Dmir_basal_fst))/2
-log(1-dpseST.dpseSR_medial_fst) * (2000000/-log(1-dpseST.Dmir_medial_fst))/2
-log(1-dpseST.dpseSR_terminal_fst) * (2000000/-log(1-dpseST.Dmir_terminal_fst))/2

#Get Fst estimates for each breakpoint
mean(dpseST.dpseSR_basal_fst_tbl[,7])
mean(dpseST.dpseSR_basal_fst_tbl[,8])
mean(dpseST.dpseSR_basal_fst_tbl[,9])

mean(dpseST.dpseSR_medial_fst_tbl[,7])
mean(dpseST.dpseSR_medial_fst_tbl[,8])
mean(dpseST.dpseSR_medial_fst_tbl[,9])

mean(dpseST.dpseSR_terminal_fst_tbl[,7])
mean(dpseST.dpseSR_terminal_fst_tbl[,8])
mean(dpseST.dpseSR_terminal_fst_tbl[,9])

-log(1-mean(dpseST.dpseSR_terminal_fst_tbl[,8])) * (2000000/-log(1-dpseST.Dmir_terminal_fst))/2
-log(1-mean(dpseST.dpseSR_terminal_fst_tbl[,9])) * (2000000/-log(1-dpseST.Dmir_terminal_fst))/2


gtf<-read.table("~/dpse_3.2.cds.subset.gtf", sep="\t") 
gene_id<-sapply(strsplit(as.character(gtf$V9), " ", fixed=TRUE), tail, 1)
gene_id<-(as.vector(unlist(as.vector(strsplit(gene_id, " ")))))
gene_id<-gsub( ";", "", as.character(gene_id))
gtf$gene_id<-gene_id
mins<-aggregate(gtf$V4, by = list(gtf$gene_id), min)
maxs<-aggregate(gtf$V5, by = list(gtf$gene_id), max)
chroms<-gtf[!duplicated(gtf$gene_id),c("V1","gene_id")]
cds_coords<-merge(mins, maxs, by="Group.1")
names(cds_coords)<-c("gene_id","cds_start","cds_stop")
cds_coords<-merge(cds_coords, chroms, by='gene_id')
names(cds_coords)<-c("gene_id","cds_start","cds_stop","chrom")

pi<-read.table("~/x_chrom_pi_estimates.txt",header=T)
pi<-merge(pi, cds_coords, by.x="product", by.y="gene_id")
pi$file<-gsub(".sorted.vcf","",pi$file)
pi$sample<-gsub("^([^_]*_[^_]*)_(.*)$", "\\2", pi$file)
dpseST.pi<-subset(pi,pi$sample=="10700X1")
dpseSR.pi<-subset(pi,pi$sample=="10700X2")
Dmir.pi<-subset(pi,pi$sample=="Dmir_SP138")

dpseST.pi<-make_AD_pi(dpseST.pi)
dpseSR.pi<-make_AD_pi(dpseSR.pi)
Dmir.pi<-make_AD_pi(Dmir.pi)
head(dpseST.pi)

mean(dpseST.pi$piS)
mean(dpseSR.pi$piS)

#Define the regions to compare based on inversion breakpoints
bas_inv_bps<-c(6162858, 8796222)
med_inv_bps<-c(9471558,18455166)
dis_inv_bps<-c(25039842,28331788)

outside_reg<-c(1,6162858-1000000)
bas_reg<-c(6162858, 8796222)
colinear1_reg<-c(8796222,9471558)
med_reg<-c(9471558,18455166)
colinear2_reg<-c(18455166,25039842)
dis_reg<-c(25039842,28331788)

get_reg_pi<-function(df){
  outside_reg<-c(1,6162858-1000000)
  bas_reg<-c(6162858, 8796222)
  colinear1_reg<-c(8796222,9471558)
  med_reg<-c(9471558,18455166)
  colinear2_reg<-c(18455166,25039842)
  dis_reg<-c(25039842,28331788)
  
  raw_df<-df
  df<-subset(raw_df, raw_df$piS>0)
  outside_piS<-subset(df, df$Position %in% seq(outside_reg[1],outside_reg[2]))[,"piS"]
  bas_piS<-subset(df, df$Position %in% seq(bas_reg[1],bas_reg[2]))[,"piS"]
  med_piS<-subset(df, df$Position %in% seq(med_reg[1],med_reg[2]))[,"piS"]
  colinear1_piS<-subset(df, df$Position %in% seq(colinear1_reg[1],colinear1_reg[2]))
  colinear2_piS<-subset(df, df$Position %in% seq(colinear2_reg[1],colinear2_reg[2]))
  dis_piS<-subset(df, df$Position %in% seq(dis_reg[1],dis_reg[2]))[,"piS"]
  colinear_piS<-rbind(colinear1_piS, colinear2_piS)[,"piS"]
  
  df<-subset(raw_df, raw_df$piN>0)
  outside_piN<-subset(df, df$Position %in% seq(outside_reg[1],outside_reg[2]))[,"piN"]
  bas_piN<-subset(df, df$Position %in% seq(bas_reg[1],bas_reg[2]))[,"piN"]
  med_piN<-subset(df, df$Position %in% seq(med_reg[1],med_reg[2]))[,"piN"]
  colinear1_piN<-subset(df, df$Position %in% seq(colinear1_reg[1],colinear1_reg[2]))
  colinear2_piN<-subset(df, df$Position %in% seq(colinear2_reg[1],colinear2_reg[2]))
  dis_piN<-subset(df, df$Position %in% seq(dis_reg[1],dis_reg[2]))[,"piN"]
  colinear_piN<-rbind(colinear1_piN, colinear2_piN)[,"piN"]
  
  return(list(as.vector(outside_piS), as.vector(bas_piS), as.vector(med_piS), as.vector(dis_piS), as.vector(colinear_piS),
              as.vector(outside_piN), as.vector(bas_piN), as.vector(med_piN), as.vector(dis_piN), as.vector(colinear_piN)))
}

par(mfrow=c(2,2))
dpseST.pi_regions<-get_reg_pi(data.frame(dpseST.pi))
dpseSR.pi_regions<-get_reg_pi(data.frame(dpseSR.pi))
boxplot(dpseST.pi_regions[1:5],ylim=c(0,0.013))
boxplot(dpseSR.pi_regions[1:5],ylim=c(0,0.013))

boxplot(dpseST.pi_regions[6:10],ylim=c(0,0.0021))
boxplot(dpseSR.pi_regions[6:10],ylim=c(0,0.0021))

wilcox.test(dpseST.pi_regions[[1]],dpseSR.pi_regions[[1]])
wilcox.test(dpseST.pi_regions[[2]],dpseSR.pi_regions[[2]])
wilcox.test(dpseST.pi_regions[[3]],dpseSR.pi_regions[[3]])
wilcox.test(dpseST.pi_regions[[4]],dpseSR.pi_regions[[4]])
wilcox.test(dpseST.pi_regions[[5]],dpseSR.pi_regions[[5]])

wilcox.test(dpseST.pi_regions[[6]],dpseSR.pi_regions[[6]])
wilcox.test(dpseST.pi_regions[[7]],dpseSR.pi_regions[[7]])
wilcox.test(dpseST.pi_regions[[8]],dpseSR.pi_regions[[8]])
wilcox.test(dpseST.pi_regions[[9]],dpseSR.pi_regions[[9]])
wilcox.test(dpseST.pi_regions[[10]],dpseSR.pi_regions[[10]])


#Make windowed Fst plot for Figure 3
head(dpseST.dpseSR_AD_fst)
get_reg_fst<-function(df){
  outside_reg<-c(1,6162858-1000000)
  bas_reg<-c(6162858, 8796222)
  colinear1_reg<-c(8796222,9471558)
  med_reg<-c(9471558,18455166)
  colinear2_reg<-c(18455166,25039842)
  dis_reg<-c(25039842,28331788)
  
  outside_fst<-subset(df, df$w_start %in% seq(outside_reg[1],outside_reg[2]))[,"fst"]
  bas_fst<-subset(df, df$w_start %in% seq(bas_reg[1],bas_reg[2]))[,"fst"]
  med_fst<-subset(df, df$w_start %in% seq(med_reg[1],med_reg[2]))[,"fst"]
  colinear1_fst<-subset(df, df$w_start %in% seq(colinear1_reg[1],colinear1_reg[2]))
  colinear2_fst<-subset(df, df$w_start %in% seq(colinear2_reg[1],colinear2_reg[2]))
  dis_fst<-subset(df, df$w_start %in% seq(dis_reg[1],dis_reg[2]))[,"fst"]
  colinear_fst<-rbind(colinear1_fst, colinear2_fst)[,"fst"]
  
  chrom_avg<-data.frame(fst=as.vector(df[,"fst"]),region="chrom")
  outside_avg<-data.frame(fst=as.vector(outside_fst),region="outside")
  bas_avg<-data.frame(fst=as.vector(bas_fst),region="basal")
  med_avg<-data.frame(fst=as.vector(med_fst),region="medial")
  dis_avg<-data.frame(fst=as.vector(dis_fst),region="distal")
  colinear_avg<-data.frame(fst=as.vector(colinear_fst),region="colinear")
  
  fst_df<-rbind(chrom_avg, outside_avg, bas_avg, med_avg, colinear_avg, dis_avg)
  
  return(fst_df)
}
dev.off()
fst_regs<-get_reg_fst(dpseST.dpseSR_AD_fst)
fst_regs$comp<-"ST_SR"
#fst_regs<-append(list(dpseST.dpseSR_AD_fst$fst),fst_regs)
ggplot(data = fst_regs) + aes(x = region, y = fst, colour=comp) +  geom_boxplot() + theme_bw() + geom_hline(yintercept=c(dpseST.dpseSR_basal_fst,dpseST.dpseSR_medial_fst,dpseST.dpseSR_terminal_fst))
boxplot(fst_regs,notch=T)
head(fst_regs)
c(dpseST.dpseSR_basal_fst,dpseST.dpseSR_medial_fst,dpseST.dpseSR_terminal_fst)
abline(h=c(dpseST.dpseSR_basal_fst,dpseST.dpseSR_medial_fst,dpseST.dpseSR_terminal_fst))

fst_inversions<-unlist(c(fst_regs[[3]],fst_regs[[4]],fst_regs[[5]]))
wilcox.test(fst_inversions,fst_regs[[6]],alternative="greater")


make_AD_windowed_pi<-function(df, window){
  scaffolds<-c("XL_group1a","XL_group1e","XL_group3a","XL_group3b","XR_group3a","XR_group5","XR_group6","XR_group8")
  
  #Chromosome XR
  sites<-df
  AD1<-subset(sites, sites$V1=="XL_group3b" & sites$V2 > 1 & sites$V2 < 386183)
  AD2<-subset(sites, sites$V1=="XL_group3a" & sites$V2 > 1362506 & sites$V2 < 2690836)
  AD3<-subset(sites, sites$V1=="XL_group1a" & sites$V2 > 2660039 & sites$V2 < 5354993)
  AD4<-subset(sites, sites$V1=="XR_group6" & sites$V2 > 1 & sites$V2 < 13314419)
  AD5<-subset(sites, sites$V1=="XR_group8" & sites$V2 > 1 & sites$V2 < 3598432)
  AD6<-subset(sites, sites$V1=="XR_group8" & sites$V2 > 3598433 & sites$V2 < 9212921)
  AD7<-subset(sites, sites$V1=="XR_group3a" & sites$V2 > 1 & sites$V2 < 1468910)
  AD8<-subset(sites, sites$V1=="XR_group5" & sites$V2 > 1 & sites$V2 < 740492)
  
  AD3<-(AD3[rev(rownames(AD3)),])
  AD6<-(AD6[rev(rownames(AD6)),])
  AD7<-(AD7[rev(rownames(AD7)),])
  AD8<-(AD8[rev(rownames(AD8)),])
  head(AD8)
  AD1$"Position"<-AD1$V2
  AD2$"Position"<-(AD2$V2 - 1362506) + 386183
  AD3$"Position"<-(2694954 - (AD3$V2 - 2660039)) + 386183 + 1328330
  AD4$"Position"<-AD4$V2 + 386183 + 1328330 + 2694954
  AD5$"Position"<-AD5$V2 + 386183 + 1328330 + 2694954 + 13314419
  AD6$"Position"<-(5614488 - (AD6$V2 - 3598433)) + 386183 + 1328330 + 2694954 + 13314419 + 3598432
  AD7$"Position"<-(1468910 - AD7$V2) + 386183 + 1328330 + 2694954 + 13314419 + 3598432 + 5614488
  AD8$"Position"<-(740492 - AD8$V2) + 386183 + 1328330 + 2694954 + 13314419 + 3598432 + 5614488 + 1468910
  AD<-rbind(AD1,AD2,AD3,AD4,AD5,AD6,AD7,AD8)
  #colnames(AD)<-c("Chrom","Chrom_start","Chrom_end","windowsize","missing_sites","dxy","w_start")
  AD<-AD[complete.cases(AD), ]
  
  names(AD)<-c("chrom","midpoint","N","callable_sites","pi","Position")
  AD<-subset(AD, !(is.na(as.numeric(as.character(AD$pi)))))
  return(AD)
}

make_AD_windowed_D<-function(df, window){
  scaffolds<-c("XL_group1a","XL_group1e","XL_group3a","XL_group3b","XR_group3a","XR_group5","XR_group6","XR_group8")
  
  #Chromosome XR
  sites<-df
  AD1<-subset(sites, sites$V1=="XL_group3b" & sites$V2 > 1 & sites$V2 < 386183)
  AD2<-subset(sites, sites$V1=="XL_group3a" & sites$V2 > 1362506 & sites$V2 < 2690836)
  AD3<-subset(sites, sites$V1=="XL_group1a" & sites$V2 > 2660039 & sites$V2 < 5354993)
  AD4<-subset(sites, sites$V1=="XR_group6" & sites$V2 > 1 & sites$V2 < 13314419)
  AD5<-subset(sites, sites$V1=="XR_group8" & sites$V2 > 1 & sites$V2 < 3598432)
  AD6<-subset(sites, sites$V1=="XR_group8" & sites$V2 > 3598433 & sites$V2 < 9212921)
  AD7<-subset(sites, sites$V1=="XR_group3a" & sites$V2 > 1 & sites$V2 < 1468910)
  AD8<-subset(sites, sites$V1=="XR_group5" & sites$V2 > 1 & sites$V2 < 740492)
  
  AD3<-(AD3[rev(rownames(AD3)),])
  AD6<-(AD6[rev(rownames(AD6)),])
  AD7<-(AD7[rev(rownames(AD7)),])
  AD8<-(AD8[rev(rownames(AD8)),])
  head(AD8)
  AD1$"Position"<-AD1$V2
  AD2$"Position"<-(AD2$V2 - 1362506) + 386183
  AD3$"Position"<-(2694954 - (AD3$V2 - 2660039)) + 386183 + 1328330
  AD4$"Position"<-AD4$V2 + 386183 + 1328330 + 2694954
  AD5$"Position"<-AD5$V2 + 386183 + 1328330 + 2694954 + 13314419
  AD6$"Position"<-(5614488 - (AD6$V2 - 3598433)) + 386183 + 1328330 + 2694954 + 13314419 + 3598432
  AD7$"Position"<-(1468910 - AD7$V2) + 386183 + 1328330 + 2694954 + 13314419 + 3598432 + 5614488
  AD8$"Position"<-(740492 - AD8$V2) + 386183 + 1328330 + 2694954 + 13314419 + 3598432 + 5614488 + 1468910
  AD<-rbind(AD1,AD2,AD3,AD4,AD5,AD6,AD7,AD8)
  #colnames(AD)<-c("Chrom","Chrom_start","Chrom_end","windowsize","missing_sites","dxy","w_start")
  AD<-AD[complete.cases(AD), ]
  
  names(AD)<-c("chrom","midpoint","N","callable_sites","D","Position")
  AD<-subset(AD, !(is.na(as.numeric(as.character(AD$D)))))
  return(AD)
}

pi_10kb_dpseST<-make_AD_windowed_pi(read.table("~/10700X1.pooled.pi"),10000)
pi_10kb_dpseSR<-make_AD_windowed_pi(read.table("~/10700X2.pooled.pi"),10000)
pi_10kb_dpseSTSR<-make_AD_windowed_pi(read.table("~/SR.merged.pooled.pi"),10000)
head(pi_10kb_dpseSTSR)


D_10kb_dpseST<-make_AD_windowed_D(read.table("~/10700X1.pooled.D"),10000)
D_10kb_dpseSR<-make_AD_windowed_D(read.table("~/10700X2.pooled.D"),10000)
D_10kb_dpseSTSR<-make_AD_windowed_D(read.table("~/SR.merged.pooled.D"),10000)

intergenic_regs<-read.table("intergenic_10kb_windows.bed")
names(intergenic_regs)<-c("chrom","from","to")
head(intergenic_regs)

library(dplyr)
pi_10kb_dpseST<-left_join(pi_10kb_dpseST, intergenic_regs) %>%
  filter(midpoint > from, midpoint < to ) %>%
  select(-from, -to)
pi_10kb_dpseST$pi<-as.numeric(as.character(pi_10kb_dpseST$pi))
pi_10kb_dpseSR<-left_join(pi_10kb_dpseSR, intergenic_regs) %>%
  filter(midpoint > from, midpoint < to ) %>%
  select(-from, -to)
pi_10kb_dpseSR$pi<-as.numeric(as.character(pi_10kb_dpseSR$pi))
pi_10kb_dpseSTSR<-left_join(pi_10kb_dpseSTSR, intergenic_regs) %>%
  filter(midpoint > from, midpoint < to ) %>%
  select(-from, -to)
pi_10kb_dpseSTSR$pi<-as.numeric(as.character(pi_10kb_dpseSTSR$pi))

D_10kb_dpseST<-left_join(D_10kb_dpseST, intergenic_regs) %>%
  filter(midpoint > from, midpoint < to ) %>%
  select(-from, -to)
D_10kb_dpseST$D<-as.numeric(as.character(D_10kb_dpseST$D))
D_10kb_dpseSR<-left_join(D_10kb_dpseSR, intergenic_regs) %>%
  filter(midpoint > from, midpoint < to ) %>%
  select(-from, -to)
D_10kb_dpseSR$D<-as.numeric(as.character(D_10kb_dpseSR$D))
D_10kb_dpseSTSR<-left_join(D_10kb_dpseSTSR, intergenic_regs) %>%
  filter(midpoint > from, midpoint < to ) %>%
  select(-from, -to)
D_10kb_dpseSTSR$D<-as.numeric(as.character(D_10kb_dpseSTSR$D))


head(pi_10kb_dpseSTSR)
dev.off()
plot(pi_10kb_dpseST$Position, pi_10kb_dpseST$pi,col="blue",ylab="Pi",xlab="Position",main="Intergenic Pi",cex=0.5)
points(pi_10kb_dpseSR$Position, pi_10kb_dpseSR$pi,col="green",cex=0.5)
points(pi_10kb_dpseSTSR$Position, pi_10kb_dpseSTSR$pi,col="black",cex=0.5)
lines(smooth.spline(pi_10kb_dpseST$Position, pi_10kb_dpseST$pi,spar=0.1),col="blue",lwd=2)
lines(smooth.spline(pi_10kb_dpseSR$Position, pi_10kb_dpseSR$pi,spar=0.1),col="green",lwd=2)
lines(smooth.spline(pi_10kb_dpseSTSR$Position, pi_10kb_dpseSTSR$pi,spar=0.1),col="black",lwd=2)

plot(D_10kb_dpseST$Position, as.numeric(as.character(D_10kb_dpseST$D)),col="blue",ylab="D",xlab="Position",main="D",cex=0.5,ylim=c(-3.5, 1.5))
points(D_10kb_dpseSR$Position, as.numeric(as.character(D_10kb_dpseSR$D)),col="green",ylab="Pi",xlab="Position",main="D",cex=0.5)
points(D_10kb_dpseSTSR$Position, as.numeric(as.character(D_10kb_dpseSTSR$D)),col="black",ylab="Pi",xlab="Position",main="D",cex=0.5)
lines(smooth.spline(D_10kb_dpseST$Position, as.numeric(as.character(D_10kb_dpseST$D)),,spar=0.1),col="blue",lwd=2)
lines(smooth.spline(D_10kb_dpseSR$Position, as.numeric(as.character(D_10kb_dpseSR$D)),spar=0.1),col="green",lwd=2)
lines(smooth.spline(D_10kb_dpseSTSR$Position,  as.numeric(as.character(D_10kb_dpseSTSR$D)),spar=0.1),col="black",lwd=2)
mean(pi_10kb_dpseST$V10,breaks=20)
mean(pi_10kb_dpseSR$V10,breaks=20)
boxplot(D_10kb_dpseST$D,D_10kb_dpseSR$D)
hist(pi_10kb_dpseSR$V10)

ST_SR_merged_pi<-merge(pi_10kb_dpseST, pi_10kb_dpseSR, by=c("V1", "Position"))
head(ST_SR_merged_pi)
plot(ST_SR_merged_pi$Position, as.numeric(as.character(ST_SR_merged_pi$V4.x))/as.numeric(as.character(ST_SR_merged_pi$V4.y)), ylim=c(-5,5))
mean(as.numeric(as.character(ST_SR_merged_pi$V4.x))/as.numeric(as.character(ST_SR_merged_pi$V4.y)))

pi_1kb_dpseST<-make_AD_windowed_pi(read.table("~/dpseST_1kb.pool_pi.txt"),1000)
pi_1kb_dpseSR<-make_AD_windowed_pi(read.table("~/dpseSR_1kb.pool_pi.txt"),1000)
pi_1kb_dpseSTSR<-make_AD_windowed_pi(read.table("~/dpseST_dpseSR_1kb.pool_pi.txt"),1000)

plot(pi_1kb_dpseST$Position, pi_1kb_dpseST$V4,col="blue",type="n")
lines(smooth.spline(pi_1kb_dpseST$Position, pi_1kb_dpseST$V4,spar=0.0005),col="blue")
lines(smooth.spline(pi_1kb_dpseSR$Position, pi_1kb_dpseSR$V4,spar=0.0005),col="green")
lines(smooth.spline(pi_1kb_dpseSTSR$Position, pi_1kb_dpseSTSR$V4,spar=0.0005),col="black")

get_reg_pi_intergenic<-function(df){
  outside_reg<-c(1,6162858-1000000)
  bas_reg<-c(6162858, 8796222)
  colinear1_reg<-c(8796222,9471558)
  med_reg<-c(9471558,18455166)
  colinear2_reg<-c(18455166,25039842)
  dis_reg<-c(25039842,28331788)
  
  outside_piS<-subset(df, df$Position %in% seq(outside_reg[1],outside_reg[2]))[,"pi"]
  bas_piS<-subset(df, df$Position %in% seq(bas_reg[1],bas_reg[2]))[,"pi"]
  med_piS<-subset(df, df$Position %in% seq(med_reg[1],med_reg[2]))[,"pi"]
  colinear1_piS<-subset(df, df$Position %in% seq(colinear1_reg[1],colinear1_reg[2]))
  colinear2_piS<-subset(df, df$Position %in% seq(colinear2_reg[1],colinear2_reg[2]))
  dis_piS<-subset(df, df$Position %in% seq(dis_reg[1],dis_reg[2]))[,"pi"]
  colinear_piS<-rbind(colinear1_piS, colinear2_piS)[,"pi"]
  
  chrom_avg<-data.frame(pi=as.vector(df[,"pi"]),region="chrom")
  outside_avg<-data.frame(pi=outside_piS,region="outside")
  bas_avg<-data.frame(pi=bas_piS,region="basal")
  med_avg<-data.frame(pi=med_piS,region="medial")
  dis_avg<-data.frame(pi=dis_piS,region="distal")
  colinear_avg<-data.frame(pi=colinear_piS,region="collinear")
  
  pi_df<-rbind(chrom_avg, outside_avg, bas_avg, med_avg, colinear_avg, dis_avg)
  return(pi_df)
}

dpseST_reg_pi<-get_reg_pi_intergenic(pi_10kb_dpseST)
dpseST_reg_pi$group<-"dpseST"
dpseSR_reg_pi<-get_reg_pi_intergenic(pi_10kb_dpseSR)
dpseSR_reg_pi$group<-"dpseSR"
dpseSTSR_reg_pi<-get_reg_pi_intergenic(pi_10kb_dpseSTSR)
dpseSTSR_reg_pi$group<-"dpseSTSR"

reg_pi<-rbind(dpseST_reg_pi, dpseSR_reg_pi, dpseSTSR_reg_pi)

ggplot(data = reg_pi) + aes(x = region, y = pi, colour = group) +  geom_boxplot() + theme_bw() 

mean(reg_pi[(reg_pi$region=="outside") & (reg_pi$group=="dpseST"),"pi"])/mean(reg_pi[(reg_pi$region=="outside") & (reg_pi$group=="dpseSR"),"pi"])
mean(reg_pi[(reg_pi$region=="basal") & (reg_pi$group=="dpseST"),"pi"])/mean(reg_pi[(reg_pi$region=="basal") & (reg_pi$group=="dpseSR"),"pi"])
mean(reg_pi[(reg_pi$region=="medial") & (reg_pi$group=="dpseST"),"pi"])/mean(reg_pi[(reg_pi$region=="medial") & (reg_pi$group=="dpseSR"),"pi"])
mean(reg_pi[(reg_pi$region=="distal") & (reg_pi$group=="dpseST"),"pi"])/mean(reg_pi[(reg_pi$region=="distal") & (reg_pi$group=="dpseSR"),"pi"])
mean(reg_pi[(reg_pi$region=="collinear") & (reg_pi$group=="dpseST"),"pi"])/mean(reg_pi[(reg_pi$region=="collinear") & (reg_pi$group=="dpseSR"),"pi"])

t.test(reg_pi[(reg_pi$region=="outside") & (reg_pi$group=="dpseST"),"pi"],reg_pi[(reg_pi$region=="outside") & (reg_pi$group=="dpseSR"),"pi"])
wilcox.test(reg_pi[(reg_pi$region=="outside") & (reg_pi$group=="dpseST"),"pi"],reg_pi[(reg_pi$region=="outside") & (reg_pi$group=="dpseSR"),"pi"])

t.test(reg_pi[(reg_pi$region=="distal") & (reg_pi$group=="dpseST"),"pi"],reg_pi[(reg_pi$region=="distal") & (reg_pi$group=="dpseSR"),"pi"])

wilcox.test(reg_pi[reg_pi$group=="dpseST","pi"],reg_pi[reg_pi$group=="dpseSR","pi"])
#Mean and 95% CI for total pi
SR_XR_pi_mean <- boot(data.frame(reg_pi[reg_pi$group=="dpseSR",]), function(d, i){mean(d[i, 'pi'])}, R = 10000)
SR_XR_pi_ci <- boot.ci(SR_XR_pi_mean, type = 'bca')
mean(reg_pi[reg_pi$group=="dpseSR","pi"])
SR_XR_pi_ci 

ST_XR_pi_mean <- boot(data.frame(reg_pi[reg_pi$group=="dpseST",]), function(d, i){mean(d[i, 'pi'])}, R = 10000)
ST_XR_pi_ci <- boot.ci(ST_XR_pi_mean, type = 'bca')
mean(reg_pi[reg_pi$group=="dpseST","pi"])
ST_XR_pi_ci 

mean(reg_pi[reg_pi$group=="dpseST","pi"]/reg_pi[reg_pi$group=="dpseSR","pi"])

get_reg_td_intergenic<-function(df){
  outside_reg<-c(1,6162858-1000000)
  bas_reg<-c(6162858, 8796222)
  colinear1_reg<-c(8796222,9471558)
  med_reg<-c(9471558,18455166)
  colinear2_reg<-c(18455166,25039842)
  dis_reg<-c(25039842,28331788)
  
  outside_td<-subset(df, df$Position %in% seq(outside_reg[1],outside_reg[2]))[,"D"]
  bas_td<-subset(df, df$Position %in% seq(bas_reg[1],bas_reg[2]))[,"D"]
  med_td<-subset(df, df$Position %in% seq(med_reg[1],med_reg[2]))[,"D"]
  colinear1_td<-subset(df, df$Position %in% seq(colinear1_reg[1],colinear1_reg[2]))
  colinear2_td<-subset(df, df$Position %in% seq(colinear2_reg[1],colinear2_reg[2]))
  dis_td<-subset(df, df$Position %in% seq(dis_reg[1],dis_reg[2]))[,"D"]
  colinear_td<-rbind(colinear1_td, colinear2_td)[,"D"]
  
  chrom_avg<-data.frame(td=as.vector(df[,"D"]),region="chrom")
  outside_avg<-data.frame(td=outside_td,region="outside")
  bas_avg<-data.frame(td=bas_td,region="basal")
  med_avg<-data.frame(td=med_td,region="medial")
  dis_avg<-data.frame(td=dis_td,region="distal")
  colinear_avg<-data.frame(td=colinear_td,region="collinear")
  
  td_df<-rbind(chrom_avg, outside_avg, bas_avg, med_avg, colinear_avg, dis_avg)
  return(td_df)
}

dpseST_reg_td<-get_reg_td_intergenic(D_10kb_dpseST)
dpseST_reg_td$group<-"dpseST"
dpseST_reg_td$td<-as.numeric(as.character(dpseST_reg_td$td))

dpseSR_reg_td<-get_reg_td_intergenic(D_10kb_dpseSR)
dpseSR_reg_td$group<-"dpseSR"
dpseSR_reg_td$td<-as.numeric(as.character(dpseSR_reg_td$td))

dpseSTSR_reg_td<-get_reg_td_intergenic(D_10kb_dpseSTSR)
dpseSTSR_reg_td$group<-"dpseSTSR"
dpseSTSR_reg_td$td<-as.numeric(as.character(dpseSTSR_reg_td$td))

reg_td<-rbind(dpseST_reg_td, dpseSR_reg_td, dpseSTSR_reg_td)
tail(reg_td)
reg_td$td<-as.numeric(as.character(reg_td$td))
ggplot(data = reg_td) + aes(x = region, y = td, colour = group) +  geom_boxplot() + theme_bw()

SR_XR_pi_mean <- boot(data.frame(reg_td[reg_td$group=="dpseSR",]), function(d, i){mean(d[i, 'td'])}, R = 10000)
SR_XR_pi_ci <- boot.ci(SR_XR_pi_mean, type = 'bca')
mean(reg_td[reg_td$group=="dpseSR","td"])
SR_XR_pi_ci 

ST_XR_pi_mean <- boot(data.frame(reg_td[reg_td$group=="dpseST",]), function(d, i){mean(d[i, 'td'])}, R = 10000)
ST_XR_pi_ci <- boot.ci(ST_XR_pi_mean, type = 'bca')
mean(reg_td[reg_td$group=="dpseST","td"])
ST_XR_pi_ci 

wilcox.test(reg_td[reg_td$group=="dpseST","td"], reg_td[reg_td$group=="dpseSR","td"])
t.test(reg_td[reg_td$group=="dpseSR","td"], reg_td[(reg_td$group=="dpseSR") & (reg_td$region=="distal"),"td"],alternative="less")
min(reg_td[reg_td$group=="dpseSR","td"])


raw_ann<-read.table("~/SR.pooled.ann.table")
head(raw_ann)
make_AD_ann<-function(df){
  scaffolds<-c("XL_group1a","XL_group1e","XL_group3a","XL_group3b","XR_group3a","XR_group5","XR_group6","XR_group8")
  
  #Chromosome XR
  sites<-df
  AD1<-subset(sites, sites$V1=="XL_group3b" & sites$V2 > 1 & sites$V2 < 386183)
  AD2<-subset(sites, sites$V1=="XL_group3a" & sites$V2 > 1362506 & sites$V2 < 2690836)
  AD3<-subset(sites, sites$V1=="XL_group1a" & sites$V2 > 2660039 & sites$V2 < 5354993)
  AD4<-subset(sites, sites$V1=="XR_group6" & sites$V2 > 1 & sites$V2 < 13314419)
  AD5<-subset(sites, sites$V1=="XR_group8" & sites$V2 > 1 & sites$V2 < 3598432)
  AD6<-subset(sites, sites$V1=="XR_group8" & sites$V2 > 3598433 & sites$V2 < 9212921)
  AD7<-subset(sites, sites$V1=="XR_group3a" & sites$V2 > 1 & sites$V2 < 1468910)
  AD8<-subset(sites, sites$V1=="XR_group5" & sites$V2 > 1 & sites$V2 < 740492)
  
  AD3<-(AD3[rev(rownames(AD3)),])
  AD6<-(AD6[rev(rownames(AD6)),])
  AD7<-(AD7[rev(rownames(AD7)),])
  AD8<-(AD8[rev(rownames(AD8)),])
  head(AD8)
  AD1$"Position"<-AD1$V2
  AD2$"Position"<-(AD2$V2 - 1362506) + 386183
  AD3$"Position"<-(2694954 - (AD3$V2 - 2660039)) + 386183 + 1328330
  AD4$"Position"<-AD4$V2 + 386183 + 1328330 + 2694954
  AD5$"Position"<-AD5$V2 + 386183 + 1328330 + 2694954 + 13314419
  AD6$"Position"<-(5614488 - (AD6$V2 - 3598433)) + 386183 + 1328330 + 2694954 + 13314419 + 3598432
  AD7$"Position"<-(1468910 - AD7$V2) + 386183 + 1328330 + 2694954 + 13314419 + 3598432 + 5614488
  AD8$"Position"<-(740492 - AD8$V2) + 386183 + 1328330 + 2694954 + 13314419 + 3598432 + 5614488 + 1468910
  AD<-rbind(AD1,AD2,AD3,AD4,AD5,AD6,AD7,AD8)
  #colnames(AD)<-c("Chrom","Chrom_start","Chrom_end","windowsize","missing_sites","dxy","w_start")
  AD<-AD[complete.cases(AD), ]
  
  #names(AD)<-c("chrom","midpoint","N","callable_sites","D","Position")
  #AD<-subset(AD, !(is.na(as.numeric(as.character(AD$D)))))
  return(AD)
}
AD_ann<-make_AD_ann(raw_ann)
tail(AD_ann)
names(AD_ann)<-c("scaff", "scaff_pos","ref","alt","ST_p","SR_p","Dmir_p","ST_q","SR_q","Dmir_q","ann","gene","csq","genes","Position")
length(AD_ann[(AD_ann$ST_q)|(AD_ann$SR_q),1])

#Get the proportion of sites that fixed on XR in 10kb windows

get_window_fx<-function(start, stop){
  ST_fx<-0
  SR_fx<-0
  window_start<-start
  window_end<-stop
  window<-AD_ann[AD_ann$Position>start & AD_ann$Position<stop,]
  if ( !is.null(window)){
    SR_fx<-length(window[window$SR_q==8 & window$ST_q==0,1])/10000
    ST_fx<-length(window[window$ST_q==8 & window$SR_q==0,1])/10000
    window_start<-window$Position[1]
    window_end<-window[length(window$Position),"Position"]   
  }
  
  return(c(as.integer(window_start),as.integer(window_end),ST_fx,SR_fx))
}

get_window_shared<-function(start, stop){
  ST_fx<-0
  SR_fx<-0
  window_start<-start
  window_end<-stop
  window<-AD_ann[AD_ann$Position>start & AD_ann$Position<stop,]
  if ( !is.null(window)){
    shr<-length(window[window$SR_q>0 & window$ST_q>0,1])/10000
    window_start<-window$Position[1]
    window_end<-window[length(window$Position),"Position"]   
  }
  
  return(c(as.integer(window_start),as.integer(window_end),shr))
}

head(AD_ann)

starts=NULL
stops=NULL
ST_fxs=NULL
SR_fxs=NULL
shared=NULL
count=1
for(i in seq(1,30000000,10000)){
  start<-i
  stop<-i+10000
  window_res<-get_window_fx(start,stop)
  window_shared<-get_window_shared(start, stop)
  print(window_shared)
  starts[count]<-window_res[1]
  stops[count]<-window_res[2]
  ST_fxs[count]<-window_res[3]
  SR_fxs[count]<-window_res[4]
  shared[count]<-window_shared[3]
  count<-count+1
}
shared
daf_fx<-data.frame(starts,stops,ST_fxs,SR_fxs,shared)
mean(daf_fx$SR_fxs,na.rm=T)
plot(daf_fx$starts, daf_fx$SR_fxs)
plot(daf_fx$starts, daf_fx$shared)
plot(daf_fx$starts, log10(daf_fx$SR_fxs/daf_fx$shared))

SR_XR_fixed <- boot(data.frame(daf_fx), function(d, i){mean(d[i, 'SR_fxs'],na.rm=T)}, R = 10000)
SR_XR_fixed_ci <- boot.ci(SR_XR_fixed, type = 'bca')
mean(daf_fx$SR_fxs,na.rm=T)
SR_XR_fixed_ci

ST_XR_fixed <- boot(data.frame(daf_fx), function(d, i){mean(d[i, 'ST_fxs'],na.rm=T)}, R = 10000)
ST_XR_fixed_ci <- boot.ci(ST_XR_fixed, type = 'bca')
mean(daf_fx$ST_fxs,na.rm=T)
ST_XR_fixed_ci

get_reg_fixed<-function(df){
  outside_reg<-c(1,6162858-1000000)
  bas_reg<-c(6162858, 8796222)
  colinear1_reg<-c(8796222,9471558)
  med_reg<-c(9471558,18455166)
  colinear2_reg<-c(18455166,25039842)
  dis_reg<-c(25039842,28331788)
  
  outside_fx<-subset(df, df$starts %in% seq(outside_reg[1],outside_reg[2]))
  bas_fx<-subset(df, df$starts %in% seq(bas_reg[1],bas_reg[2]))
  med_fx<-subset(df, df$starts %in% seq(med_reg[1],med_reg[2]))
  colinear1_fx<-subset(df, df$starts %in% seq(colinear1_reg[1],colinear1_reg[2]))
  colinear2_fx<-subset(df, df$starts %in% seq(colinear2_reg[1],colinear2_reg[2]))
  dis_fx<-subset(df, df$starts %in% seq(dis_reg[1],dis_reg[2]))
  colinear_fx<-rbind(colinear1_fx, colinear2_fx)
  
  chrom_avg<-data.frame(fx=df$SR_fxs/df$shared,region="chrom")
  outside_avg<-data.frame(fx=outside_fx$SR_fxs/outside_fx$shared,region="outside")
  bas_avg<-data.frame(fx=bas_fx$SR_fxs/bas_fx$shared,region="basal")
  med_avg<-data.frame(fx=med_fx$SR_fxs/med_fx$shared,region="medial")
  dis_avg<-data.frame(fx=dis_fx$SR_fxs/dis_fx$shared,region="distal")
  colinear_avg<-data.frame(fx=colinear_fx$SR_fxs/colinear_fx$shared,region="collinear")
  
  td_df<-rbind(chrom_avg, outside_avg, bas_avg, med_avg, dis_avg, colinear_avg)
  return(td_df)
}
fx_reg<-get_reg_fixed(daf_fx)
fx_reg<-fx_reg[(is.finite(fx_reg$fx)),]
fx_reg$fx<-as.numeric(as.character(fx_reg$fx))

aggregate(fx_reg, by=list(fx_reg$region), FUN=mean)

#Find the number of fixed, derived sites that are found within protein coding genes
fixed_sites<-subset(AD_ann,(AD_ann$SR_q==8 & AD_ann$ST_q==0 & AD_ann$Dmir_q==0)|(AD_ann$SR_q==0 & AD_ann$ST_q==8 & AD_ann$Dmir_q==1))

shared_sites<-subset(AD_ann,(AD_ann$SR_q>0 & AD_ann$ST_q>0))
plot(shared_sites$Position)
unique(fixed_sites$AllTypes)
nonint<-subset(fixed_sites,fixed_sites$ann!="I")
length(nonint[,1])
nonsyn<-subset(fixed_sites,fixed_sites$ann=="N")

cis_utrs<-fixed_sites[grep("5_prime_UTR_variant",fixed_sites$csq),]
length(unique(cis_utrs$gene))
table(fixed_sites$ann)
lof<-subset(fixed_sites,fixed_sites$ann=="L")
length(unique(nonsyn$gene))
head(nonsyn)






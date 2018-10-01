#External libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("RUVSeq")
#install.packages("statmod")
library(RUVSeq)
library(RColorBrewer)
library(gplots)
library(statmod)

#Distance Function
dist.pear<-function(x) as.dist(1-cor(t(x)))

#Color Function
colfunc<-colorRampPalette(c("blue","gray","red"))

##Read in files and reformat names
sr_rc<-read.delim("data/SR_fc_3.2.txt",row.names=1,comment.char="#")
#sr_rc<-read.delim("sr_fc.txt",row.names=1,comment.char="#")
colnames(sr_rc)<-gsub("\\_sequence.txt.fastq.subjunc.BAM","",names(sr_rc))
colnames(sr_rc)<-gsub("*reads.","",names(sr_rc))
colnames(sr_rc)<-gsub("*D_pseudoobscura_Re.Sequencing_PI_SW_Schaeffer_Sex_Ratio_Strains_Sample_","",names(sr_rc))
colnames(sr_rc)<-gsub("_[A-Z][1-9].*","",names(sr_rc))
colnames(sr_rc)<-gsub("PSU.*","",names(sr_rc))
colnames(sr_rc)<-gsub("[.]","",names(sr_rc))
colnames(sr_rc)<-gsub("alignments","",names(sr_rc))
names(sr_rc)

#Filter genes by requiring more than 10 reads in at least three samples
gene_counts<-sr_rc[,-(1:5)]
filter<-apply(gene_counts,1,function(x) length(x[x>10])>=3)
filtered<-gene_counts[filter,]
filter_merge<-merge(filtered,sr_rc,by=0)
genes<-rownames(filtered)
length(genes)

#Set factors of groups for DE analysis
colnames(filtered)<-gsub("s_","s",names(filtered))
colnames(filtered)
x<-as.factor(gsub("*[0-9]","",names(filtered)))
x
#filtered<-subset(filtered, rownames(filtered) %in% all_XR$Row.names)
#Store data in object of S4 class from EDASeq package
set<-newSeqExpressionSet(as.matrix(filtered),phenoData=data.frame(x,row.names=colnames(filtered)))
set

#Make expression level and PCA plots
colors<-brewer.pal(12,"Set3")
plotRLE(set,outline=FALSE,ylim=c(-4,4),col=colors[x])
plotPCA(set,col=colors[x],cex=.8)

#Normalize between lanes
set<-betweenLaneNormalization(set, which="upper")
plotRLE(set,outline=FALSE,ylim=c(-4,4),col=colors[x])
plotPCA(set,col=colors[x],cex=.8)
colnames(set)

differences<-makeGroups(x)
set2<-RUVs(set, rownames(filtered), k=1, differences)
plotRLE(set1,outline=FALSE,ylim=c(-4,4),col=colors[x])
plotPCA(set1,col=colors[x],cex=.8)

#Create model design matrix
design<-model.matrix(~x+W_1,data=pData(set2))
design

#Normalize factors and estimate dispersion
y<-DGEList(counts=counts(set2), group=x)
y<-calcNormFactors(y, method="RLE")
y<-estimateDisp(y, design, robust=TRUE)
y<-estimateGLMTagwiseDisp(y, design)
plotBCV(y)
fit<-glmQLFit(y, design,robust=TRUE)
lrt<-glmQLFTest(fit, coef=2)
sig_table<-lrt$table
sig_table$FDR<-p.adjust(lrt$table$PValue, method="BH")
sig<-subset(sig_table,sig_table$FDR<0.05)
length(sig[,1])
topTags(lrt)
AD_sig_genes<-(subset(AD_exp,AD_exp$FDR<0.05)$Row.names)
#cat(rownames(sig),sep="\n")

#Merge with all genes and create XR and NonXR tables
all_merged<-merge(sr_rc,sig_table,by=0)
all_XR<-(subset(all_merged,gsub("*_.*","",all_merged$Chr)=="XR"))
all_XL<-(subset(all_merged,gsub("*_.*","",all_merged$Chr)=="XL"))
all_nonXR<-(subset(all_merged,gsub("*_.*","",all_merged$Chr)!="XR"))
head(all_merged)

length(genes)
length(all_XR[,1])
length(all_XL[,1])

#Merge significant genes for Chr3 and NonMullerC
merged<-merge(sr_rc,sig,by=0)
sig_XR<-(subset(merged,gsub("*_.*","",merged$Chr)=="XR"))
sig_nonXR<-(subset(merged,gsub("*_.*","",merged$Chr)!="XR"))
sig_XL<-(subset(merged,gsub("*_.*","",merged$Chr)=="XL"))
rpkm<-rpkm(y,gene.length=filter_merge$Length)
sig_mat<-subset(rpkm,row.names(rpkm) %in% row.names(sig))
sig_XRmat<-subset(rpkm,row.names(rpkm) %in% (subset(AD_exp,AD_exp$FDR<0.05)$Row.names))

#Plot Heatmap
dev.off()
heatmap.2(log2(sig_XRmat[,1:12]+1),col=colfunc(300),trace='none',
          ,scale='row',distfun=dist.pear,cexCol=0.5,labRow=F)
lrt
#subset(lrt,rownames(lrt$table) %in% all_XR$Row.names)
#See if there are more upregulated or down regulared on the chromosome
de1 <- decideTestsDGE(lrt[subset(AD_exp,AD_exp$FDR<0.05)$Row.names,], adjust.method="BH", p.value=0.05)
summary(de1)
cat(subset(AD_exp,AD_exp$FDR<0.05)[de1=="1",]$Row.names,sep="\n")
de1tags12 <- rownames(lrt[all_XR$Row.names,])[as.logical(de1)] 
plotSmear(lrt[all_XR$Row.names,],de.tags=de1tags12,ylim=c(-10,10),cex=.75)
p.adjust(lrt$table$PValue, method="BH")
volcanoData <- cbind(lrt[all_XR$Row.names,]$table$logFC, -log(p.adjust(lrt[all_XR$Row.names,]$table$PValue,method="BH")))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19, col = ifelse(volcanoData[,2] > 3,'red','black'),xlim=c(-10,10))
sig_table[rownames(sig_table)=="FBgn0072656",]
sig_table[rownames(sig_table)=="FBgn0077431",]
sig_table[rownames(sig_table)=="FBgn0077730",]
sig_table[rownames(sig_table)=="FBgn0077782",]
sig_table[rownames(sig_table)=="FBgn0082126",]
sig_table[rownames(sig_table)=="FBgn0249618",]
sig_FC<-(subset(volcanoData,volcanoData[,2] > 3))

mean(subset(sig_FC,sig_FC[,1] < 1)[,1])

std <- function(x) sd(x)/sqrt(length(x))
gene_exp<-function(gene_name){
  avg_1<-mean(log2(sig_XRmat[rownames(sig_XRmat)==gene_name,][1:6]+1))
  avg_2<-mean(log2(sig_XRmat[rownames(sig_XRmat)==gene_name,][7:12]+1))
  se_1<-std(log2(sig_XRmat[rownames(sig_XRmat)==gene_name,][1:6]+1))
  se_2<-std(log2(sig_XRmat[rownames(sig_XRmat)==gene_name,][7:12]+1))
  return(c(avg_1,avg_2,se_1,se_2))
}
ST_means<-c(gene_exp("FBgn0072656")[1],gene_exp("FBgn0077431")[1],gene_exp("FBgn0077730")[1],gene_exp("FBgn0077782")[1],gene_exp("FBgn0082126")[1],gene_exp("FBgn0249618")[1])
SR_means<-c(gene_exp("FBgn0072656")[2],gene_exp("FBgn0077431")[2],gene_exp("FBgn0077730")[2],gene_exp("FBgn0077782")[2],gene_exp("FBgn0082126")[2],gene_exp("FBgn0249618")[2])
ST_se<-c(gene_exp("FBgn0072656")[3],gene_exp("FBgn0077431")[3],gene_exp("FBgn0077730")[3],gene_exp("FBgn0077782")[3],gene_exp("FBgn0082126")[3],gene_exp("FBgn0249618")[3])
SR_se<-c(gene_exp("FBgn0072656")[4],gene_exp("FBgn0077431")[4],gene_exp("FBgn0077730")[4],gene_exp("FBgn0077782")[4],gene_exp("FBgn0082126")[4],gene_exp("FBgn0249618")[4])
barcenters<-barplot(height=rbind(ST_means,SR_means),beside=TRUE,ylim=c(0,5))
segments(barcenters,rbind(ST_means,SR_means)-rbind(ST_se,SR_se)*2,barcenters,rbind(ST_means,SR_means)+rbind(ST_se,SR_se)*2)
arrows(barcenters,rbind(ST_means,SR_means)-rbind(ST_se,SR_se)*2,barcenters,rbind(ST_means,SR_means)+rbind(ST_se,SR_se)*2,lwd=1.5,angle=90,code=3,length=0.05)


sig_XR$Chr<-gsub("*;.*","",sig_XR$Chr)
sig_XR$Start<-gsub("*;.*","",sig_XR$Start)
sig_XR$End<-gsub(".*;","",sig_XR$End)
sig_XL$Chr<-gsub("*;.*","",sig_XL$Chr)
sig_XL$Start<-gsub("*;.*","",sig_XL$Start)
sig_XL$End<-gsub(".*;","",sig_XL$End)

all_genes_lrt<-sig_table
all_genes_lrt<-merge(sr_rc,all_genes_lrt,by=0)
all_genes_lrt$Chr<-gsub("*;.*","",all_genes_lrt$Chr)
all_genes_lrt$Start<-gsub("*;.*","",all_genes_lrt$Start)
all_genes_lrt$End<-gsub(".*;","",all_genes_lrt$End)
head(all_genes_lrt)
XR_genes<-subset(all_genes_lrt,gsub("*_.*","",all_genes_lrt$Chr)=="XR")

#Muller AD (Chr XR) Genes
scafs<-c(0,386183,1714514,4409469,17723888,26936807,28405715,29141086)
AD_1<-subset(all_genes_lrt,all_genes_lrt$Chr=="XL_group3b" & as.numeric(all_genes_lrt$Start) > 1 & as.numeric(all_genes_lrt$End) < 386183)
AD_2<-subset(all_genes_lrt,all_genes_lrt$Chr=="XL_group3a" & as.numeric(all_genes_lrt$Start) > 1362506 & as.numeric(all_genes_lrt$End) < 2690836)
AD_3<-subset(all_genes_lrt,all_genes_lrt$Chr=="XL_group1a" & as.numeric(all_genes_lrt$Start) > 2660039 & as.numeric(all_genes_lrt$End) < 5354993)
AD_4<-subset(all_genes_lrt,all_genes_lrt$Chr=="XR_group6" & as.numeric(all_genes_lrt$Start) > 1 & as.numeric(all_genes_lrt$End) < 13314419)
AD_5<-subset(all_genes_lrt,all_genes_lrt$Chr=="XR_group8" & as.numeric(all_genes_lrt$Start) > 1 & as.numeric(all_genes_lrt$End) < 3598432)
AD_6<-subset(all_genes_lrt,all_genes_lrt$Chr=="XR_group8" & as.numeric(all_genes_lrt$Start) > 3598433 & as.numeric(all_genes_lrt$End) < 9212921)
AD_7<-subset(all_genes_lrt,all_genes_lrt$Chr=="XR_group3a" & as.numeric(all_genes_lrt$Start) > 1 & as.numeric(all_genes_lrt$End) < 1468910)
AD_8<-subset(all_genes_lrt,all_genes_lrt$Chr=="XR_group5" & as.numeric(all_genes_lrt$Start) > 1 & as.numeric(all_genes_lrt$End) < 740492)

AD_1$Position<-as.numeric(AD_1$Start)
AD_2$Position<-as.numeric(AD_2$Start)-1362506+386183
AD_3$Position<-(5354993-(as.numeric(AD_3$Start))+1328330+386183)
AD_4$Position<-as.numeric(AD_4$Start)+2694954+1328330+386183
AD_5$Position<-as.numeric(AD_5$Start)+13314419+2694954+1328330+386183
AD_6$Position<-(9212921-as.numeric(AD_6$Start))+3598432+13314419+2694954+1328330+386183
AD_7$Position<-(1468910-as.numeric(AD_7$Start))+5614488+3598432+13314419+2694954+1328330+386183
AD_8$Position<-(740492-as.numeric(AD_8$Start))+1468910+5614488+3598432+13314419+2694954+1328330+386183
AD_exp<-rbind(AD_1,AD_2,AD_3,AD_4,AD_5,AD_6,AD_7,AD_8)
plot(AD_exp$Position,-log(AD_exp$FDR),cex=.75,col=ifelse(-log(AD_exp$FDR) > 3,'red','gray'),pch=ifelse(-log(AD_exp$FDR) > 3,19,1))
abline(v=scafs,col="gray",lty=5)
#plot(AD_exp$Position,abs(AD_exp$logFC),cex=.75)
SR_inv<-c(6157103,8790276,9465432,18430099,23175372,28903141)
inverted_size<-(8790276-6157103) + (18430099-9465432) + (28903141 - 23175372) 
collienar_size <- 29146149 - inverted_size - 6157103

abline(v=SR_inv,col="black")
bas_inv_genes<-subset(AD_exp,AD_exp$Position>=6157103 & AD_exp$Position<=8790276)
med_inv_genes<-subset(AD_exp,AD_exp$Position>=9465432 & AD_exp$Position<=18430099)
dis_inv_genes<-subset(AD_exp,AD_exp$Position>=23175372 & AD_exp$Position<=28903141)
inv_genes<-rbind(bas_inv_genes,med_inv_genes,dis_inv_genes)
noninv_genes<-subset(AD_exp,!(AD_exp$Row.names %in% inv_genes$Row.names) & AD_exp$Position>6157103)
cat(AD_sig_genes,sep="\n")
length(unique(subset(inv_genes, inv_genes$FDR<0.05)[,1]))
length(unique(subset(noninv_genes, noninv_genes$FDR<0.05)[,1]))

#Muller AD (Chr XR) Genes
#scafs<-c(0,386183,1714514,4409469,17723888,26936807,28405715,29141086)
A_1<-subset(all_genes_lrt,all_genes_lrt$Chr=="XL_group1a" & as.numeric(all_genes_lrt$Start) > 714715 & as.numeric(all_genes_lrt$End) < 2660038)
A_2<-subset(all_genes_lrt,all_genes_lrt$Chr=="XL_group1a" & as.numeric(all_genes_lrt$Start) > 1 & as.numeric(all_genes_lrt$End) < 711448)
A_3<-subset(all_genes_lrt,all_genes_lrt$Chr=="XL_group1e" & as.numeric(all_genes_lrt$Start) > 1 & as.numeric(all_genes_lrt$End) < 12253060)
A_4<-subset(all_genes_lrt,all_genes_lrt$Chr=="XL_group1a" & as.numeric(all_genes_lrt$Start) > 5354994 & as.numeric(all_genes_lrt$End) < 9151669)
A_5<-subset(all_genes_lrt,all_genes_lrt$Chr=="XL_group3a" & as.numeric(all_genes_lrt$Start) > 1 & as.numeric(all_genes_lrt$End) < 1362505)

A_1$Position<-2660038-as.numeric(A_1$Start)
A_2$Position<-as.numeric(A_2$Start)+1945323
A_3$Position<-as.numeric(A_3$Start)+1945323+711448
A_4$Position<-(as.numeric(A_4$Start)-5354994)+1945323+711448+12253060
A_5$Position<-(1362505-as.numeric(A_5$Start))+1945323+711448+12253060+3796675
A_exp<-rbind(A_1,A_2,A_3,A_4,A_5)
plot(A_exp$Position,-log(A_exp$FDR),cex=.75,col=ifelse(-log(A_exp$FDR) > 3,'red','gray'),pch=ifelse(-log(A_exp$FDR) > 3,19,1))
length(subset(AD_exp,AD_exp$FDR<0.05)[,1])

inv_scafs<-rbind(AD_4,AD_5,AD_6,AD_7,AD_8)
noninv_scafs<-rbind(AD_1,AD_2,AD_3)
length(subset(inv_scafs,inv_scafs$FDR<0.05)[,1])
length(subset(noninv_scafs,noninv_scafs$FDR<0.05)[,1])

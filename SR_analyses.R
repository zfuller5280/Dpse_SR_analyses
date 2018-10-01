#Library
library(boot)

#Read in the data
AD_dxy<-read.table("data/XR_dpseSTdpseSR_stats.txt",sep="\t",header=T)
A_dxy<-read.table("data/XL_dpseSTdpseSR_stats.txt",sep="\t",header=T)
head(A_dxy)
AD_dxy[AD_dxy==0] <- .000001
AD_dxy[AD_dxy==1] <- (1-.000001)

#Inversion breakpoints
bas_inv_bps<-c(6157103, 8790276)
med_inv_bps<-c(9465432,18430099)
dis_inv_bps<-c(23175372,28903141)
23175372-18430099

#Subset the regions with 250kb surrounding the breakpoints
bas_inv_reg<-subset(AD_dxy,((AD_dxy$w_start > bas_inv_bps[1] - 125000) & (AD_dxy$w_start < bas_inv_bps[1] + 125000))|
                      ((AD_dxy$w_start > bas_inv_bps[2] - 125000) & (AD_dxy$w_start < bas_inv_bps[2] + 125000)))
med_inv_reg<-subset(AD_dxy,((AD_dxy$w_start > med_inv_bps[1] - 125000) & (AD_dxy$w_start < med_inv_bps[1] + 125000))|
                      ((AD_dxy$w_start > med_inv_bps[2] - 125000) & (AD_dxy$w_start < med_inv_bps[2] + 125000)))
dis_inv_reg<-subset(AD_dxy,((AD_dxy$w_start > dis_inv_bps[1] - 125000) & (AD_dxy$w_start < dis_inv_bps[1] + 125000))|
                      ((AD_dxy$w_start > dis_inv_bps[2] - 125000) & (AD_dxy$w_start < dis_inv_bps[2] + 125000)))
#Put all the inverted regions together
all_inv<-rbind(bas_inv_reg,med_inv_reg,dis_inv_reg)

bas_inv_var<-subset(AD_var_table,((AD_var_table$Position > bas_inv_bps[1] - 125000) & (AD_var_table$Position < bas_inv_bps[1] + 125000))|
                      ((AD_var_table$Position > bas_inv_bps[2] - 125000) & (AD_var_table$Position < bas_inv_bps[2] + 125000)))
med_inv_var<-subset(AD_var_table,((AD_var_table$Position > med_inv_bps[1] - 125000) & (AD_var_table$Position < med_inv_bps[1] + 125000))|
                      ((AD_var_table$Position > med_inv_bps[2] - 125000) & (AD_var_table$Position < med_inv_bps[2] + 125000)))
dis_inv_var<-subset(AD_var_table,((AD_var_table$Position > dis_inv_bps[1] - 125000) & (AD_var_table$Position < dis_inv_bps[1] + 125000))|
                      ((AD_var_table$Position > dis_inv_bps[2] - 125000) & (AD_var_table$Position < dis_inv_bps[2] + 125000)))
med_inv_genes<-(unique(med_inv_var$Gene))
dis_inv_genes<-(unique(dis_inv_var$Gene))

#Fst transformed estimates of divergence times
mean(-log(1-bas_inv_reg$DpseST.DpseSR_fst) * 2000000/(mean(-log(1-bas_inv_reg$DpseST.Dmir_fst))))
mean(-log(1-med_inv_reg$DpseST.DpseSR_fst) * 2000000/(mean(-log(1-med_inv_reg$DpseST.Dmir_fst))))
mean(-log(1-dis_inv_reg$DpseST.DpseSR_fst) * 2000000/(mean(-log(1-dis_inv_reg$DpseST.Dmir_fst))))
t.test(-log(1-all_inv$DpseST.DpseSR_fst) * 2000000/(mean(-log(1-all_inv$DpseST.Dmir_fst))))

#Mean and 95% CI for Dpse ST and SR Fst in all inverted breakpoint regions
boot_mean <- boot(data.frame(all_inv), function(d, i){mean(d[i, 'DpseST.DpseSR_fst'])}, R = 10000)
boot_ci <- boot.ci(boot_mean, type = 'bca')
mean(all_inv$DpseST.DpseSR_fst)
boot_ci

#Mean and 95% CI for Dpse ST and Dmir Fst in all inverted breakpoint regions
boot_mean <- boot(data.frame(all_inv), function(d, i){mean(d[i, 'DpseST.Dmir_fst'])}, R = 10000)
boot_ci <- boot.ci(boot_mean, type = 'bca')
mean(all_inv$DpseST.Dmir_fst)
boot_ci

#Mean and 95% CI for Dpse ST and SR Fst transformed estimates of divergence time
fst_age<-(-log(1-all_inv$DpseST.DpseSR_fst) * 2000000/(mean(-log(1-all_inv$DpseST.Dmir_fst))))
boot_mean <- boot(data.frame(fst_age), function(d, i){mean(d[i, 'fst_age'])}, R = 10000)
boot_ci <- boot.ci(boot_mean, type = 'bca')
mean(fst_age)
boot_ci

#Mean and 95% CI for Dpse ST and SR dxy in all inverted breakpoint regions
boot_mean <- boot(data.frame(all_inv), function(d, i){mean(d[i, 'DpseST.DpseSR_dxy'])}, R = 10000)
boot_ci <- boot.ci(boot_mean, type = 'bca')
mean(all_inv$DpseST.DpseSR_dxy)
boot_ci

#Mean and 95% CI for Dpse ST and Dmir dxy in all inverted breakpoint regions
boot_mean <- boot(data.frame(all_inv), function(d, i){mean(d[i, 'DpseST.Dmir_dxy'])}, R = 10000)
boot_ci <- boot.ci(boot_mean, type = 'bca')
mean(all_inv$DpseST.Dmir_dxy)
boot_ci

#Test if divergence is greater between species than between arrangements
wilcox.test(all_inv$DpseST.Dmir_dxy,all_inv$DpseST.DpseSR_dxy,alternative="greater")

#Mean and 95% CI for Dpse ST and SR Fst in basal inverted region
basal_boot_mean <- boot(data.frame(bas_inv_reg), function(d, i){mean(d[i, 'DpseST.DpseSR_fst'])}, R = 10000)
basal_boot_ci <- boot.ci(basal_boot_mean, type = 'bca')
mean(bas_inv_reg$DpseST.DpseSR_fst)
basal_boot_ci

#Mean and 95% CI for Dpse ST and SR Fst in medial inverted region
medial_boot_mean <- boot(data.frame(med_inv_reg), function(d, i){mean(d[i, 'DpseST.DpseSR_fst'])}, R = 10000)
medial_boot_ci <- boot.ci(medial_boot_mean, type = 'bca')
mean(med_inv_reg$DpseST.DpseSR_fst)
medial_boot_ci

#Mean and 95% CI for Dpse ST and SR Fst in terminal inverted region
terminal_boot_mean <- boot(data.frame(dis_inv_reg), function(d, i){mean(d[i, 'DpseST.DpseSR_fst'])}, R = 10000)
terminal_boot_ci <- boot.ci(terminal_boot_mean, type = 'bca')
mean(dis_inv_reg$DpseST.DpseSR_fst)
terminal_boot_ci

#Test for difference in Fst between medial and basal inversion
wilcox.test(med_inv_reg$DpseST.DpseSR_fst,bas_inv_reg$DpseST.DpseSR_fst)
#Test for difference in Fst between medial and terminal inversion
wilcox.test(med_inv_reg$DpseST.DpseSR_fst,dis_inv_reg$DpseST.DpseSR_fst)
#Test for difference in Fst between basal and terminal inversion
wilcox.test(bas_inv_reg$DpseST.DpseSR_fst,dis_inv_reg$DpseST.DpseSR_fst)

#Mean and 95% CI for Dpse ST and SR dxy in basal inverted region
basal_boot_mean <- boot(data.frame(bas_inv_reg), function(d, i){mean(d[i, 'DpseST.DpseSR_dxy'])}, R = 10000)
basal_boot_ci <- boot.ci(basal_boot_mean, type = 'bca')
mean(bas_inv_reg$DpseST.DpseSR_dxy)
basal_boot_ci

#Mean and 95% CI for Dpse ST and SR dxy in medial inverted region
medial_boot_mean <- boot(data.frame(med_inv_reg), function(d, i){mean(d[i, 'DpseST.DpseSR_dxy'])}, R = 10000)
medial_boot_ci <- boot.ci(medial_boot_mean, type = 'bca')
mean(med_inv_reg$DpseST.DpseSR_dxy)
medial_boot_ci

#Mean and 95% CI for Dpse ST and SR dxy in terminal inverted region
terminal_boot_mean <- boot(data.frame(dis_inv_reg), function(d, i){mean(d[i, 'DpseST.DpseSR_dxy'])}, R = 10000)
terminal_boot_ci <- boot.ci(terminal_boot_mean, type = 'bca')
mean(dis_inv_reg$DpseST.DpseSR_dxy)
terminal_boot_ci

#Test for difference in dxy between medial and basal inversion
wilcox.test(med_inv_reg$DpseST.DpseSR_dxy,bas_inv_reg$DpseST.DpseSR_dxy)
#Test for difference in dxy between medial and terminal inversion
wilcox.test(med_inv_reg$DpseST.DpseSR_dxy,dis_inv_reg$DpseST.DpseSR_dxy)
#Test for difference in dxy between basal and terminal inversion
wilcox.test(bas_inv_reg$DpseST.DpseSR_dxy,dis_inv_reg$DpseST.DpseSR_dxy)

#Mean and 95% CI for Dpse ST and SR Fst transformed estimate of divergence time in terminal inversion
terminal_fst_age<-(-log(1-dis_inv_reg$DpseST.DpseSR_fst) * 2000000/(mean(-log(1-dis_inv_reg$DpseST.Dmir_fst))))
terminal_boot_mean <- boot(data.frame(terminal_fst_age), function(d, i){mean(d[i, 'terminal_fst_age'])}, R = 10000)
terminal_boot_ci <- boot.ci(terminal_boot_mean, type = 'bca')
mean(terminal_fst_age)
terminal_boot_ci

#Get regions of all inversions and collinear regions of XR
inverted_regions<-subset(AD_dxy,((AD_dxy$w_start > bas_inv_bps[1]) & (AD_dxy$w_start < bas_inv_bps[2])) | ((AD_dxy$w_start > med_inv_bps[1]) & (AD_dxy$w_start < med_inv_bps[2])) |
                           ((AD_dxy$w_start > dis_inv_bps[1]) & (AD_dxy$w_start < dis_inv_bps[2])))
colinear_regions<-subset(AD_dxy,!(AD_dxy$w_start %in% inverted_regions$w_start) & AD_dxy$w_start>bas_inv_bps[1])

basal_regions<-subset(AD_dxy,((AD_dxy$w_start > bas_inv_bps[1]) & (AD_dxy$w_start < bas_inv_bps[2])))
medial_regions<-subset(AD_dxy,((AD_dxy$w_start > med_inv_bps[1]) & (AD_dxy$w_start < med_inv_bps[2])))
terminal_regions<-subset(AD_dxy,((AD_dxy$w_start > dis_inv_bps[1]) & (AD_dxy$w_start < dis_inv_bps[2])))
#Make new column for distance to closest inversion breakpoint within each inverted region
basal_regions$dist_to_bp<-pmin(abs(basal_regions$w_start - bas_inv_bps[1]),abs(basal_regions$w_stop - bas_inv_bps[2]))
medial_regions$dist_to_bp<-pmin(abs(medial_regions$w_start - med_inv_bps[1]),abs(medial_regions$w_stop - med_inv_bps[2]))
terminal_regions$dist_to_bp<-pmin(abs(terminal_regions$w_start - dis_inv_bps[1]),abs(terminal_regions$w_stop - dis_inv_bps[2]))

#Get correlation coefficients for distance to nearest bp and fst (also output in the summary of the linear model below)
cor(basal_regions$dist_to_bp,basal_regions$DpseST.DpseSR_fst)^2
cor(medial_regions$dist_to_bp,medial_regions$DpseST.DpseSR_fst)^2
cor(terminal_regions$dist_to_bp,terminal_regions$DpseST.DpseSR_fst)^2
#Get lm stats for distance to nearest bp and fst and test significance of relationship
summary(lm(basal_regions$DpseST.DpseSR_fst~basal_regions$dist_to_bp))
summary(lm(medial_regions$DpseST.DpseSR_fst~medial_regions$dist_to_bp))
summary(lm(terminal_regions$DpseST.DpseSR_fst~terminal_regions$dist_to_bp))

#Get Fst in center 500kb of basal and medial inverted regions
central_basal<-(bas_inv_bps[2] - bas_inv_bps[1])/2 + bas_inv_bps[1]
central_basal_region<-subset(basal_regions,basal_regions$w_start > central_basal - 250000 & basal_regions$w_stop < central_basal + 250000)
mean(central_basal_region$DpseST.DpseSR_fst)
central_medial<-(med_inv_bps[2] - med_inv_bps[1])/2 + med_inv_bps[1]
central_medial_region<-subset(medial_regions,medial_regions$w_start > central_medial - 250000 & medial_regions$w_stop < central_medial + 250000)
mean(central_medial_region$DpseST.DpseSR_fst)

#Test for differences in Fst between collinear regions and inverted regions and compare means
t.test(inverted_regions$DpseST.DpseSR_fst,colinear_regions$DpseST.DpseSR_fst) #Just report mean values
wilcox.test(inverted_regions$DpseST.DpseSR_fst,colinear_regions$DpseST.DpseSR_fst,alternative = "greater") #Test for difference

#Make plots of Fst for Figure 2. Read in Fst table with confidence intervals
par(mfrow=c(2,1))
plot(AD_dxy$w_start, AD_dxy$DpseST.DpseSR_fst, col="gray", cex=0.5,ylim=c(0,0.5),xlim=c(0,29200000))
u_ci_avg<-smooth.spline(AD_dxy$w_start, AD_dxy$DpseST.DpseSR_fst_ci_u,spar=0.005)
lines(u_ci_avg,lty=3,lwd=2)
l_ci_avg<-smooth.spline(AD_dxy$w_start, AD_dxy$DpseST.DpseSR_fst_ci_l,spar=0.005)
lines(l_ci_avg,lty=3,lwd=2)
polygon(c(l_ci_avg$x, rev(l_ci_avg$x)), c(u_ci_avg$y, rev(l_ci_avg$y)), col = "grey", border = NA) 
lines(smooth.spline(AD_dxy$w_start, AD_dxy$DpseST.DpseSR_fst,spar=0.005),lwd=2)
abline(v=c(bas_inv_bps[1],bas_inv_bps[2],med_inv_bps[1],med_inv_bps[2],dis_inv_bps[1],dis_inv_bps[2]))

plot(A_dxy$w_start, A_dxy$DpseST.DpseSR_fst, col="gray", cex=0.5,ylim=c(0,0.5),xlim=c(0,29200000))
u_ci_avg<-smooth.spline(A_dxy$w_start, A_dxy$DpseST.DpseSR_fst_ci_u,spar=0.005)
lines(u_ci_avg,lty=3,lwd=2)
l_ci_avg<-smooth.spline(A_dxy$w_start, A_dxy$DpseST.DpseSR_fst_ci_l,spar=0.005)
lines(l_ci_avg,lty=3,lwd=2)
polygon(c(l_ci_avg$x, rev(l_ci_avg$x)), c(u_ci_avg$y, rev(l_ci_avg$y)), col = "grey", border = NA) 
lines(smooth.spline(A_dxy$w_start, A_dxy$DpseST.DpseSR_fst,spar=0.005),lwd=2)

#Mean and 95% CI for Dpse ST and SR Fst on XL
XL_boot_mean <- boot(data.frame(A_dxy), function(d, i){mean(d[i, 'DpseST.DpseSR_fst'])}, R = 10000)
XL_boot_ci <- boot.ci(XL_boot_mean, type = 'bca')
mean(A_dxy$DpseST.DpseSR_fst)
XL_boot_ci

#Get Fst in proximal region approaching basal inversion on XR
prox_boot_mean<-boot(data.frame(AD_dxy[AD_dxy$w_start<bas_inv_bps[1],]),function(d, i){mean(d[i, 'DpseST.DpseSR_fst'])}, R = 10000)
prox_boot_ci<-boot.ci(prox_boot_mean, type='bca')
mean((AD_dxy[AD_dxy$w_start<bas_inv_bps[1],'DpseST.DpseSR_fst']))
prox_boot_ci
wilcox.test((AD_dxy[AD_dxy$w_start<bas_inv_bps[1],'DpseST.DpseSR_fst']),colinear_regions$DpseST.DpseSR_fst)

#Make boxplots for figure
dev.off()
par(mfrow=c(2,1))
boxplot(A_dxy$DpseST.DpseSR_fst,AD_dxy$DpseST.DpseSR_fst,colinear_regions$DpseST.DpseSR_fst,
        basal_regions$DpseST.DpseSR_fst,medial_regions$DpseST.DpseSR_fst,terminal_regions$DpseST.DpseSR_fst,notch=T)
abline(h=c(mean(bas_inv_reg$DpseST.DpseSR_fst),mean(med_inv_reg$DpseST.DpseSR_fst),mean(dis_inv_reg$DpseST.DpseSR_fst)))
c(mean(bas_inv_reg$DpseST.DpseSR_fst),mean(med_inv_reg$DpseST.DpseSR_fst),mean(dis_inv_reg$DpseST.DpseSR_fst))

boxplot(A_dxy$DpseST.DpseSR_dxy,AD_dxy$DpseST.DpseSR_dxy,colinear_regions$DpseST.DpseSR_dxy,
        basal_regions$DpseST.DpseSR_dxy,medial_regions$DpseST.DpseSR_dxy,terminal_regions$DpseST.DpseSR_dxy,notch=T,ylim=c(0,0.015))
abline(h=c(mean(bas_inv_reg$DpseST.DpseSR_dxy),mean(med_inv_reg$DpseST.DpseSR_dxy),mean(dis_inv_reg$DpseST.DpseSR_dxy)))
c(mean(bas_inv_reg$DpseST.DpseSR_dxy),mean(med_inv_reg$DpseST.DpseSR_dxy),mean(dis_inv_reg$DpseST.DpseSR_dxy))


#Get the proportion of sites that fixed on XR in 10kb windows
AD_var_table<-read.table("data/XR.ann.var.table.csv", sep="\t",header=T)

get_window_fx<-function(start, stop){
  ST_fx<-0
  SR_fx<-0
  window_start<-start
  window_end<-stop
  window<-AD_var_table[AD_var_table$Position>start & AD_var_table$Position<stop,]
  if ( !is.null(window)){
    SR_fx<-length(window[window$SR_q==8 & window$ST_q==0,1])/10000
    ST_fx<-length(window[window$ST_q==8 & window$SR_q==0,1])/10000
    window_start<-window$Position[1]
    window_end<-window[length(window$Position),"Position"]   
  }
  
  return(c(as.integer(window_start),as.integer(window_end),ST_fx,SR_fx))
}

starts=NULL
stops=NULL
ST_fxs=NULL
SR_fxs=NULL
count=1
for(i in seq(1,30000000,10000)){
  start<-i
  stop<-i+10000
  window_res<-get_window_fx(start,stop)
  print(window_res)
  starts[count]<-window_res[1]
  stops[count]<-window_res[2]
  ST_fxs[count]<-window_res[3]
  SR_fxs[count]<-window_res[4]
  count<-count+1
}

daf_fx<-data.frame(starts,stops,ST_fxs,SR_fxs)

#Get chromosome-wide mean and CI of proportion of fixed sites
fix_boot_mean<-boot(daf_fx,function(d, i){mean(d[i, 'SR_fxs'],na.rm=T)}, R = 10000)
fix_boot_ci<-boot.ci(fix_boot_mean, type='bca')
mean(daf_fx$SR_fxs,na.rm=T)
fix_boot_ci$bca

#Get inverted region mean and CI of proportion of fixed sites
inverted_regions_fx<-subset(daf_fx,((daf_fx$starts > bas_inv_bps[1]) & (daf_fx$starts < bas_inv_bps[2])) | ((daf_fx$starts > med_inv_bps[1]) & (daf_fx$starts < med_inv_bps[2])) |
                               ((daf_fx$starts > dis_inv_bps[1]) & (daf_fx$starts < dis_inv_bps[2])))
colinear_regions_fx<-subset(daf_fx,!(daf_fx$starts %in% inverted_regions$starts) & daf_fx$starts>bas_inv_bps[1])
fixInv_boot_mean<-boot(inverted_regions_fx,function(d, i){mean(d[i, 'SR_fxs'],na.rm=T)}, R = 10000)
fixInv_boot_ci<-boot.ci(fixInv_boot_mean, type='bca')
mean(inverted_regions_fx$SR_fxs,na.rm=T)
fixInv_boot_ci$bca

#Test for difference in proportion of fixed sites between inverted regions and collinear regions
wilcox.test(inverted_regions_fx$SR_fxs,colinear_regions_fx$SR_fxs)

#Find the 10kb window with the greatest number of fixed derived sites unique to SR
head(daf_fx[order(daf_fx$SR_fxs,decreasing=T),])

#Find the number of fixed, derived sites that are found within protein coding genes
fixed_sites<-subset(AD_var_table,AD_var_table$SR_q==8 & AD_var_table$ST_q==0)
head(fixed_sites)
nonint<-subset(fixed_sites,fixed_sites$Type!="I")
length(nonint[,1])
nonsyn<-subset(fixed_sites,fixed_sites$Type=="N")
lof<-subset(fixed_sites,fixed_sites$Type=="L")
length(unique(nonsyn$Gene))

length(lof[,1])
#Get the number of genes that contain multiple fixed amino acid differences
nonsyn$PL<-1
nsyn_sums<-aggregate(nonsyn$PL, by=list(Gene=nonsyn$Gene),FUN=sum)
head(nsyn_sums[order(nsyn_sums$x,decreasing=TRUE),])
length(nsyn_sums[nsyn_sums$x>1,1])

#Get the number of fixed derived sites that are in intergenic (non-coding) regions
nonint<-subset(fixed_sites,fixed_sites$Type!="I")

agg_counts<-aggregate(Type~Gene,data=nonint,FUN=length)
head(agg_counts)
med_agg<-subset(agg_counts,agg_counts$Gene %in% med_inv_genes)
dis_agg<-subset(agg_counts,agg_counts$Gene %in% dis_inv_genes)
tail(med_agg[order(as.numeric(med_agg$Type)),])
tail(dis_agg[order(as.numeric(dis_agg$Type)),])

#Get proportion of fixed sites for XL
A_var_table<-read.table("data/XL.ann.var.table.csv", sep="\t",header=T)
get_A_window_fx<-function(start, stop){
  ST_fx<-0
  SR_fx<-0
  window_start<-start
  window_end<-stop
  window<-A_var_table[A_var_table$Position>start & A_var_table$Position<stop,]
  
  if ( !is.null(window)){
    SR_fx<-length(window[window$SR_q==8 & window$ST_q==0,1])/10000
    ST_fx<-length(window[window$ST_q==8 & window$SR_q==0,1])/10000
    window_start<-window$Position[1]
    window_end<-window[length(window$Position),"Position"]    
  }
  return(c(as.integer(window_start),as.integer(window_end),ST_fx,SR_fx))
}

A_starts=NULL
A_stops=NULL
ST_A_fxs=NULL
SR_A_fxs=NULL
count=1
for(i in seq(1,20068944,10000)){
  start<-i
  stop<-i+10000
  window_res<-get_A_window_fx(start,stop)
  print(window_res)
  A_starts[count]<-window_res[1]
  A_stops[count]<-window_res[2]
  ST_A_fxs[count]<-window_res[3]
  SR_A_fxs[count]<-window_res[4]
  count<-count+1
}
A_daf_fx<-data.frame(A_starts,A_stops,ST_A_fxs,SR_A_fxs)
head(A_daf_fx)

fixed_sites_A<-subset(A_var_table,A_var_table$SR_q==8 & A_var_table$ST_q==0)
head(fixed_sites_A)
nonint_A<-subset(fixed_sites_A,fixed_sites_A$Type!="I")
length(nonint_A[,1])
nonsyn<-subset(fixed_sites_A,fixed_sites_A$Type=="N")
lof<-subset(fixed_sites,fixed_sites$Type=="L")
length(unique(nonsyn$Gene))





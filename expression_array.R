
require(Biobase)
library(Mfuzz)
library(goseq)
library(GO.db)
library(gplots)
library(ggplot2)
library(siggenes)
options(bitmapType="cairo")
options(scipen=999)
library(limma)

# READ IN data
data<-read.table("Illu-Quant.txt",row.names=1,header=T,sep="\t")

# Parsing names
colnames(data)<-gsub("AVG_Signal.HCT116.","",colnames(data))
colnames(data)<-gsub("p007","WT",colnames(data))
colnames(data)<-gsub("\\.\\.",".",colnames(data),perl=T)

# Design test
#lev <- c("wt.0hr","wt.6hr","wt.24hr","mu.0hr","mu.6hr","mu.24hr")
#f <- factor(c("wt.0hr","wt.0hr","wt.6hr","wt.24hr","mu.0hr","mu.0hr","mu.6hr","mu.24hr"), levels=lev)
#design <- model.matrix(~0+f)
#colnames(design) <- lev

# Design matrix
samples <-gsub(".R.+","",colnames(data),perl=T)
f <- factor (samples, levels=unique(samples))
design <- model.matrix(~0+f)
colnames(design) <- unique(samples)

# ESET transformation
eset<-new("ExpressionSet", exprs=data.matrix(data))

# linear model FIT
fit <- lmFit(eset, design)

# genes respond at either the 24 hour or 48 hour or 72 hour times in the wild-type? 
cont.wt <- makeContrasts(
      "WT.24h-WT.DMSO",
      "WT.48h-WT.DMSO",
      "WT.72h-WT.DMSO",
  levels=design)
 fit_wt <- contrasts.fit(fit, cont.wt)
 fit_wt <- eBayes(fit_wt)
 wt_table=topTableF(fit_wt,number=35000, adjust="BH")
 table(wt_table$adj.P.Val<0.05)
wtnames=rownames(wt_table[wt_table$adj.P.Val<0.05,])

# Which genes respond differently over time in the mutants relative to the wild-type?
cont.mt <- makeContrasts(
     p53_24h=(Tp53.24h-Tp53.DMSO)-(WT.24h-WT.DMSO),
     p53_48h=(Tp53.48h-Tp53.DMSO)-(WT.48h-WT.DMSO),
     p53_72h=(Tp53.72h-Tp53.DMSO)-(WT.72h-WT.DMSO),
     dnmt1_24h=(DNMT1.24h-DNMT1.DMSO)-(WT.24h-WT.DMSO),
     dnmt1_48h=(DNMT1.48h-DNMT1.DMSO)-(WT.48h-WT.DMSO),
     dnmt1_72h=(DNMT1.72h-DNMT1.DMSO)-(WT.72h-WT.DMSO),
     levels=design)
 fit_mt <- contrasts.fit(fit, cont.mt)
 fit_mt <- eBayes(fit_mt)
 mt_table=topTableF(fit_mt,number=35000, adjust="BH")
 table(mt_table$adj.P.Val<0.05)
mtnames=rownames(mt_table[mt_table$adj.P.Val<0.05,])
saveRDS(mtnames,"mtnames.RDS")

# Fuzzy clustering on wt
data<-read.table("~/Downloads/Illu-Quant.txt",row.names=1,header=T,sep="\t")
colnames(data)<-gsub("AVG_Signal.HCT116.","",colnames(data))
colnames(data)<-gsub("p007","WT",colnames(data))
colnames(data)<-gsub("\\.\\.",".",colnames(data),perl=T)

mtnames<-readRDS("~/Downloads/mtnames.RDS")

wt=data[rownames(data) %in% mtnames,]
wt=(data.matrix(wt[,1:4])+data.matrix(wt[,13:16]))/2
wt<-new("ExpressionSet", exprs=wt)
wt.s<-standardise(wt)
cl_wt<-mfuzz(wt.s,c=7,m=mestimate(wt.s))
pdf('mfuzz.pdf')
mfuzz.plot(wt.s,cl=cl_wt,mfrow=c(3,3),new.window=T)
dev.off()

saveRDS(cl_wt,"cl_wt.rds")

# ploting trajectories
cl_wt=readRDS("cl_wt.rds")
data<-read.table("~/Downloads/Illu-Quant.txt",row.names=1,header=T,sep="\t")
colnames(data)<-gsub("AVG_Signal.HCT116.","",colnames(data))
colnames(data)<-gsub("p007","WT",colnames(data))
colnames(data)<-gsub("\\.\\.",".",colnames(data),perl=T)

mtnames<-readRDS("~/Downloads/mtnames.RDS")

plot.trajectories2<-function(n){
common=rownames(data) %in% (names(cl_wt$cluster)[cl_wt$cluster==n])
wt1<-data[common,1:4]
wt2<-data[common,13:16]
dm1<-data[common,5:8]
dm2<-data[common,17:20]
tp1<-data[common,9:12]
tp2<-data[common,21:24]

confi_wt1=apply(wt1,2,function(x) t.test(x)$conf.int)
confi_wt2=apply(wt2,2,function(x) t.test(x)$conf.int)
confi_dm1=apply(dm1,2,function(x) t.test(x)$conf.int)
confi_dm2=apply(dm2,2,function(x) t.test(x)$conf.int)
confi_tp1=apply(tp1,2,function(x) t.test(x)$conf.int)
confi_tp2=apply(tp2,2,function(x) t.test(x)$conf.int)
c.min=min(confi_wt1,confi_wt2,confi_dm1,confi_dm2,confi_tp1,confi_tp2)
c.max=max(confi_wt1,confi_wt2,confi_dm1,confi_dm2,confi_tp1,confi_tp2)
c.name=paste("Cluster",n)

plot(apply(wt1,2,mean),main=paste(c.name,"Trajectories"),type="l",ylim=c(c.min,c.max),xaxt='n',xlab="",ylab="")
polygon( c(1,2,3,4,rev(c(1,2,3,4)) ),c(confi_wt1[1,],rev(confi_wt1[2,])), col = alpha('salmon',.2), border = NA)
lines(apply(wt1,2,mean),col="salmon",lwd=1.6)
polygon( c(1,2,3,4,rev(c(1,2,3,4)) ),c(confi_wt2[1,],rev(confi_wt2[2,])), col = alpha('darkred',.2), border = NA)
lines(apply(wt2,2,mean),col="darkred",lwd=1.6)

polygon( c(1,2,3,4,rev(c(1,2,3,4)) ),c(confi_dm1[1,],rev(confi_dm1[2,])), col = alpha('lightcyan',.2), border = NA)
lines(apply(dm1,2,mean),col="lightcyan",lwd=1.6)
polygon( c(1,2,3,4,rev(c(1,2,3,4)) ),c(confi_dm2[1,],rev(confi_dm2[2,])), col = alpha('darkblue',.2), border = NA)
lines(apply(dm2,2,mean),col="darkblue",lwd=1.6)

polygon( c(1,2,3,4,rev(c(1,2,3,4)) ),c(confi_tp1[1,],rev(confi_tp1[2,])), col = alpha('olivedrab1',.2), border = NA)
lines(apply(tp1,2,mean),col="olivedrab1",lwd=1.6)
polygon( c(1,2,3,4,rev(c(1,2,3,4)) ),c(confi_tp2[1,],rev(confi_tp2[2,])), col = alpha('darkgreen',.2), border = NA)
lines(apply(tp2,2,mean),col="darkgreen",lwd=1.6)
Axis(side=1, labels=c(0,24,48,72),at=c(1,2,3,4))
legend("topright", inset=c(-0.2,0), legend=c("Wt","","DNMT1","","TP53",""), pch="_",col=c('salmon','darkred','lightcyan','darkblue','olivedrab1','darkgreen'), bty = "n",cex = 0.65)
}

pdf("wt_p53KD_expression_c7_sameScale.pdf")
	par(mfrow=c(3,3),c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
	for(j in 1:7){
	plot.trajectories2(j)}

dev.off()

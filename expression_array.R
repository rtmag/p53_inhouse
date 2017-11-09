
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

# Fuzzy clustering
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

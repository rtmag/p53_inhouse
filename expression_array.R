
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

#### Which genes respond differently at 48 and 72h, in the mutants relative to the wild-type?
cont.mt <- makeContrasts(
     p53_48h=(Tp53.48h-Tp53.DMSO)-(WT.48h-WT.DMSO),
     levels=design)
 fit_mt <- contrasts.fit(fit, cont.mt)
 fit_mt <- eBayes(fit_mt)
 mt_table=topTableF(fit_mt,number=35000, adjust="BH")
 table(mt_table$adj.P.Val<0.01)
mtnames= rownames(mt_table[mt_table$adj.P.Val<0.01,]) 

cont.mt <- makeContrasts(
     p53_72h=(Tp53.72h-Tp53.DMSO)-(WT.72h-WT.DMSO),
     levels=design)
 fit_mt <- contrasts.fit(fit, cont.mt)
 fit_mt <- eBayes(fit_mt)
 mt_table=topTableF(fit_mt,number=35000, adjust="BH")
 table(mt_table$adj.P.Val<0.01)
 mtnames=mtnames[mtnames %in% rownames(mt_table[mt_table$adj.P.Val<0.01,]) ]

 cont.mt <- makeContrasts(
     dnmt1_48h=(DNMT1.48h-DNMT1.DMSO)-(WT.48h-WT.DMSO),
     levels=design)
 fit_mt <- contrasts.fit(fit, cont.mt)
 fit_mt <- eBayes(fit_mt)
 mt_table=topTableF(fit_mt,number=35000, adjust="BH")
 table(mt_table$adj.P.Val<0.01)
 mtnames=mtnames[mtnames %in% rownames(mt_table[mt_table$adj.P.Val<0.01,]) ]
 
 cont.mt <- makeContrasts(
     dnmt1_72h=(DNMT1.72h-DNMT1.DMSO)-(WT.72h-WT.DMSO),
     levels=design)
 fit_mt <- contrasts.fit(fit, cont.mt)
 fit_mt <- eBayes(fit_mt)
 mt_table=topTableF(fit_mt,number=35000, adjust="BH")
 table(mt_table$adj.P.Val<0.01)
 mtnames=mtnames[mtnames %in% rownames(mt_table[mt_table$adj.P.Val<0.01,]) ]
 
##########

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
c.name=paste("Cluster",n,"Trajectories","\n(",sum(common),"genes )")

plot(apply(wt1,2,mean),main=paste(c.name),type="l",ylim=c(c.min,c.max),xaxt='n',xlab="",ylab="")
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
Axis(side=1, labels=c("0h","24h","48h","72h"),at=c(1,2,3,4))
#legend("topright", inset=c(-0.2,0), legend=c("Wt","","DNMT1","","TP53",""), pch="_",col=c('salmon','darkred','lightcyan','darkblue','olivedrab1','darkgreen'), bty = "n",cex = 0.65)
}

pdf("wt_p53KD_expression_c7_sameScale.pdf")
par(mfrow=c(3,3),c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
	for(j in 1:7){
	plot.trajectories2(j)}
	
plot.new()
legend("center",  legend=c("WT","DNMT1-KO","TP53-KO"), fill=c('darkred','darkblue','darkgreen'), bty = "n",cex = 2.3)
dev.off()

# GO	
	
common<-cl_wt$cluster==2
common<-as.numeric(common)
names(common)<-names(cl_wt$cluster)
pwf=nullp(common,"hg19","geneSymbol")
GO.BP=goseq(pwf,"hg19","geneSymbol",test.cats=c("GO:BP"))
q.val=p.adjust(GO.BP$over_represented_pvalue,method="BH")
indix=q.val<.05
GO.BP$term[indix]
write.table(x=GO.BP$term[indix],file="cluster2_GO.txt",row.names=F,col.names=F,quote=F,sep="\t")
		
common<-cl_wt$cluster==6
common<-as.numeric(common)
names(common)<-names(cl_wt$cluster)
pwf=nullp(common,"hg19","geneSymbol")
GO.BP=goseq(pwf,"hg19","geneSymbol",test.cats=c("GO:BP"))
q.val=p.adjust(GO.BP$over_represented_pvalue,method="BH")
indix=q.val<.05
GO.BP$term[indix]
write.table(x=GO.BP$term[indix],file="cluster6_GO.txt",row.names=F,col.names=F,quote=F,sep="\t")
		
common<-cl_wt$cluster==7
common<-as.numeric(common)
names(common)<-names(cl_wt$cluster)
pwf=nullp(common,"hg19","geneSymbol")
GO.BP=goseq(pwf,"hg19","geneSymbol",test.cats=c("GO:BP"))
q.val=p.adjust(GO.BP$over_represented_pvalue,method="BH")
indix=q.val<.05
GO.BP$term[indix]
write.table(x=GO.BP$term[indix],file="cluster7_GO.txt",row.names=F,col.names=F,quote=F,sep="\t")

# Write up
write.table(names(which(cl_wt$cluster==1)),file="cluster1_names.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(names(which(cl_wt$cluster==2)),file="cluster2_names.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(names(which(cl_wt$cluster==3)),file="cluster3_names.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(names(which(cl_wt$cluster==4)),file="cluster4_names.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(names(which(cl_wt$cluster==5)),file="cluster5_names.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(names(which(cl_wt$cluster==6)),file="cluster6_names.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(names(which(cl_wt$cluster==7)),file="cluster7_names.txt",row.names=F,col.names=F,quote=F,sep="\t")

### PLOT GENE TRAJECTORY

		
		
plot.gene.trajectories<-function(n,legend=FALSE){
common=rownames(data)==n
wt1<-data[common,1:4]
wt2<-data[common,13:16]
dm1<-data[common,5:8]
dm2<-data[common,17:20]
tp1<-data[common,9:12]
tp2<-data[common,21:24]

c.min=min(wt1,wt2,dm1,dm2,tp1,tp2)
c.max=max(wt1,wt2,dm1,dm2,tp1,tp2)
c.name=paste(n)

plot(apply(wt1,2,mean),main=paste(c.name,"after Doxorubicin tr"),type="l",ylim=c(c.min,c.max),xaxt='n',xlab="",ylab="")
lines(apply(wt1,2,mean),col="salmon",lwd=1.6)
lines(apply(wt2,2,mean),col="darkred",lwd=1.6)
	
lines(apply(dm1,2,mean),col="lightcyan",lwd=1.6)
lines(apply(dm2,2,mean),col="darkblue",lwd=1.6)

lines(apply(tp1,2,mean),col="olivedrab1",lwd=1.6)
lines(apply(tp2,2,mean),col="darkgreen",lwd=1.6)
Axis(side=1, labels=c('0h','24h','48h','72h'),at=c(1,2,3,4))
if(legend==TRUE){
	legend("topright", legend=c("Wt","TP53KO","DNMT1KO"), fill=c('salmon','olivedrab1','darkblue'),cex=.6, bty = "n")}
}
	
pdf("Genes_trajectories_diffScale.pdf")
par(mfrow=c(3,2))
plot.gene.trajectories("RAD51",legend=T)
plot.gene.trajectories("LMNB1",legend=T)
plot.gene.trajectories("LMNB2",legend=T)
plot.gene.trajectories("SOD1",legend=T)
plot.gene.trajectories("SOD2",legend=T)
plot.gene.trajectories("SLC9A1",legend=T)
dev.off()


pdf("Genes_trajectories_lowestQvalue.pdf")
par(mfrow=c(3,2))
plot.gene.trajectories("CTSL2",legend=T)
plot.gene.trajectories("CDCA5",legend=T)
plot.gene.trajectories("LOC399942",legend=T)
plot.gene.trajectories("CDKN3",legend=T)
plot.gene.trajectories("RRM2",legend=T)
plot.gene.trajectories("FEN1",legend=T)
plot.gene.trajectories("C11ORF82",legend=T)
plot.gene.trajectories("CDC2",legend=T)
plot.gene.trajectories("BIRC5",legend=T)
plot.gene.trajectories("OIP5",legend=T)
plot.gene.trajectories("RFC4",legend=T)
plot.gene.trajectories("LOC653874",legend=T)
plot.gene.trajectories("TK1",legend=T)
plot.gene.trajectories("CDC45L",legend=T)
plot.gene.trajectories("FAM64A",legend=T)
plot.gene.trajectories("GINS2",legend=T)
plot.gene.trajectories("KPNA2",legend=T)
plot.gene.trajectories("TMSB15A",legend=T)
plot.gene.trajectories("DSCC1",legend=T)
plot.gene.trajectories("TMEM194A",legend=T)
plot.gene.trajectories("LOC729816",legend=T)
plot.gene.trajectories("TUBA1A",legend=T)
plot.gene.trajectories("FANCI",legend=T)
plot.gene.trajectories("E2F2",legend=T)
plot.gene.trajectories("CDK2",legend=T)
plot.gene.trajectories("EXO1",legend=T)
plot.gene.trajectories("APOBEC3B",legend=T)
plot.gene.trajectories("EZH2",legend=T)
plot.gene.trajectories("FOXM1",legend=T)
plot.gene.trajectories("NUDT1",legend=T)
dev.off()

		
######################################################
plot.gene.trajectories.small<-function(n){
common=rownames(data)==n
wt1<-data[common,1:4]
wt2<-data[common,13:16]
dm1<-data[common,5:8]
dm2<-data[common,17:20]
tp1<-data[common,9:12]
tp2<-data[common,21:24]

c.min=min(wt1,wt2,dm1,dm2,tp1,tp2)
c.max=max(wt1,wt2,dm1,dm2,tp1,tp2)
c.name=paste(n)

plot(apply(wt1,2,mean),main=paste(c.name),type="l",ylim=c(c.min,c.max),xaxt='n',yaxt='n',xlab="",ylab="")
lines(apply(wt1,2,mean),col="salmon",lwd=1.6)
lines(apply(wt2,2,mean),col="darkred",lwd=1.6)
	
lines(apply(dm1,2,mean),col="lightcyan",lwd=1.6)
lines(apply(dm2,2,mean),col="darkblue",lwd=1.6)

lines(apply(tp1,2,mean),col="olivedrab1",lwd=1.6)
lines(apply(tp2,2,mean),col="darkgreen",lwd=1.6)
}		

# par(mar=c(5.1,4.1,4.1,2.1))
lsgenes = names(cl_wt$cluster)[cl_wt$cluster==2 | cl_wt$cluster==5]
lsdata = data[rownames(data) %in% lsgenes,]
wt_48<-rowMeans(lsdata[,c(3,15)])
wt_72<-rowMeans(lsdata[,c(4,16)])
dm_48<-rowMeans(lsdata[,c(7,19)])
dm_72<-rowMeans(lsdata[,c(8,20)])
tp_48<-rowMeans(lsdata[,c(11,23)])
tp_72<-rowMeans(lsdata[,c(12,24)])
		
lsgenes = lsgenes[which( (dm_48 > wt_48) & (tp_48 > wt_48) & (dm_72 > wt_72) & (tp_72 > wt_72) ) ]

		
pdf("all_genes_c2_c5.pdf")
par(mfrow=c(10,10))
par(mar=c(1.1,1.1,1.1,1.1))
for(j in 1:length(lsgenes)){
    plot.gene.trajectories.small(lsgenes[j])
}
dev.off()

write.table(lsgenes,file="cluster2_5_lsgenes_names.txt",row.names=F,col.names=F,quote=F,sep="\t")
		
###
c134 = read.csv("c1_3_4_reactome.csv")
c134_top = c134[1:15,c(7)]
names(c134_top) = c134[1:15,c(2)]
c134_top = sort(-log10(c134_top))
names(c134_top[2]) = "TP53 regulates transcription of additional cell cycle genes\nwhose exact role in the p53 pathway remain uncertain"
		
pdf("c134_pathway.pdf",width=16)		
par(mar=c(2.1,42.1,1.1,1.1))
barplot(c134_top,horiz=T,las=2,xlim=c(0,5))
dev.off()
		
		

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

# Design matrix
samples <-gsub(".R.+","",colnames(data),perl=T)
f <- factor (samples, levels=unique(samples))
design <- model.matrix(~0+f)
colnames(design) <- unique(samples)

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
 table(wt_table$adj.P.Val<0.01)
mtnames=rownames(wt_table[wt_table$adj.P.Val<0.01,])

# Fuzzy clustering on wt
wt=data[rownames(data) %in% mtnames,]
wt=(data.matrix(wt[,1:4])+data.matrix(wt[,13:16]))/2
wt<-new("ExpressionSet", exprs=wt)
wt.s<-standardise(wt)
cl_wt<-mfuzz(wt.s,c=16,m=mestimate(wt.s))
pdf('mfuzz.pdf')
mfuzz.plot(wt.s,cl=cl_wt,mfrow=c(4,4),new.window=F,time.labels=c("0h","24h","48h","72h"))
dev.off()

saveRDS(cl_wt,"cl_wt.rds")

###################################

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
par(mfrow=c(4,4),c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
for(j in 1:16){
	plot.trajectories2(j)}
dev.off()

# Write up
for(j in 1:16){
write.table(names(which(cl_wt$cluster==j)),file=paste("cluster",j,"_GeneNames.txt"sep=""),row.names=F,col.names=F,quote=F,sep="\t")
}

pdf("Color_Labels.pdf")
plot.new()
legend("center",  legend=c("WT","DNMT1-KO","TP53-KO"), fill=c('darkred','darkblue','darkgreen'), bty = "n",cex = 2.3)
dev.off()

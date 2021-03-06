require(Biobase)
library(Mfuzz)
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
#cl_wt<-mfuzz(wt.s,c=16,m=mestimate(wt.s))
#saveRDS(cl_wt,"cl_wt.rds")
cl_wt = readRDS("cl_wt.rds")
pdf('mfuzz.pdf')
mfuzz.plot(wt.s,cl=cl_wt,mfrow=c(4,4),new.window=F,time.labels=c("0h","24h","48h","72h"))
dev.off()



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
write.table(names(which(cl_wt$cluster==j)),file=paste("cluster",j,"_GeneNames.txt",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
}

pdf("Color_Labels.pdf")
plot.new()
legend("center",  legend=c("WT","DNMT1-KO","TP53-KO"), fill=c('darkred','darkblue','darkgreen'), bty = "n",cex = 2.3)
dev.off()

################################
#individual genes

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

if( abs(c.min - c.max)<3 ){ 
    c.min = mean(c(c.min,c.max))-(1.5)
    c.max = mean(c(c.min,c.max))+(1.5)
}

c.name=paste(n)

plot(apply(wt1,2,mean),main=paste(c.name,""),type="l",ylim=c(c.min,c.max),xaxt='n',xlab="",ylab="")
lines(apply(wt1,2,mean),col="red",lwd=1.6)
lines(apply(wt2,2,mean),col="red",lwd=1.6)
	
lines(apply(dm1,2,mean),col="blue",lwd=1.6)
lines(apply(dm2,2,mean),col="blue",lwd=1.6)

lines(apply(tp1,2,mean),col="green",lwd=1.6)
lines(apply(tp2,2,mean),col="green",lwd=1.6)
Axis(side=1, labels=c('0h','24h','48h','72h'),at=c(1,2,3,4))
if(legend==TRUE){
	legend("topright", legend=c("Wt","TP53KO","DNMT1KO"), fill=c('salmon','olivedrab1','darkblue'),cex=.6, bty = "n")}
}

options(bitmapType="cairo")
options(scipen=999)
data<-read.table("Illu-Quant.txt",row.names=1,header=T,sep="\t")

pdf("genes_ofInterest.pdf")
par(mfrow=c(3,3))
plot.gene.trajectories("CCNB1",legend=F)
plot.gene.trajectories("MAPK3",legend=F)
plot.gene.trajectories("MAPK1",legend=F)

plot.gene.trajectories("EIF2AK2",legend=F)
plot.gene.trajectories("ITGA1",legend=F)
plot.gene.trajectories("ITGB1",legend=F)
		
plot.gene.trajectories("DNMT1",legend=F)
plot.gene.trajectories("EZH2",legend=F)
plot.gene.trajectories("TP53",legend=F)
#
plot.gene.trajectories("IL6",legend=F)
plot.gene.trajectories("IL8",legend=F)
plot.gene.trajectories("IL10",legend=F)
		
plot.gene.trajectories("IL18",legend=F)
plot.gene.trajectories("E2F1",legend=F)
plot.gene.trajectories("E2F2",legend=F)

plot.gene.trajectories("E2F3",legend=F)
plot.gene.trajectories("E2F4",legend=F)
plot.gene.trajectories("E2F5",legend=F)
#
plot.gene.trajectories("E2F6",legend=F)
plot.gene.trajectories("E2F7",legend=F)
plot.gene.trajectories("E2F8",legend=F)

plot.gene.trajectories("BRCA1",legend=F)
plot.gene.trajectories("BRCA2",legend=F)
plot.gene.trajectories("IKBKB",legend=F)

plot.gene.trajectories("CDKN2A",legend=F)
plot.gene.trajectories("CDKN2B",legend=F)
plot.gene.trajectories("CDKN1A",legend=F)
dev.off()
		
pdf("Color_Labels_gene_ofInterest.pdf")
plot.new()
legend("center",  legend=c("WT","DNMT1-KO","TP53-KO"), fill=c('red','blue','green'), bty = "n",cex = 2.3)
dev.off()

pdf("cluster2_genes.pdf")
par(mfrow=c(3,3))
cl_2 = names(which(cl_wt$cluster==2))
for( i in 1:length(cl_2) ){
	plot.gene.trajectories(cl_2[i],legend=F)
	}
dev.off()

##
pdf("HPV_response_genes.pdf")
par(mfrow=c(3,4))
plot.gene.trajectories("KRT6C",legend=F)
plot.gene.trajectories("SERPINA3",legend=F)
plot.gene.trajectories("GJB6",legend=F)

plot.gene.trajectories("FBXO32",legend=F)
plot.gene.trajectories("FABP5",legend=F)
plot.gene.trajectories("FBXO32",legend=F)
		
plot.gene.trajectories("SPRR2F",legend=F)
plot.gene.trajectories("KRT6B",legend=F)
plot.gene.trajectories("TUBB6",legend=F)

plot.gene.trajectories("NT5DC2",legend=F)
plot.gene.trajectories("KRT6B",legend=F)
plot.gene.trajectories("AKR1B10",legend=F)
dev.off()

par(mfrow=c(3,3))
plot.gene.trajectories("MACROD1",legend=F)
plot.gene.trajectories("RARA",legend=F)
plot.gene.trajectories("SP1",legend=F)

plot.gene.trajectories("ESR1",legend=F)
plot.gene.trajectories("JUN",legend=F)
plot.gene.trajectories("EGR1",legend=F)

plot.gene.trajectories("SP1",legend=F)
plot.gene.trajectories("RELA",legend=F)
plot.gene.trajectories("NFKB1",legend=F)

pdf("interesting_genes_hg19_liftover.pdf")		
par(mfrow=c(3,3))
plot.gene.trajectories("BCAR3",legend=F)
plot.gene.trajectories("MRPL24",legend=F)
plot.gene.trajectories("MYOF",legend=F)

plot.gene.trajectories("MACROD1",legend=F)
plot.gene.trajectories("TMBIM4",legend=F)
plot.gene.trajectories("TMBIM4",legend=F)

plot.gene.trajectories("GOLGA8B",legend=F)
plot.gene.trajectories("CLYBL",legend=F)
plot.gene.trajectories("LAMA5",legend=F)
		
plot.gene.trajectories("KLHL24",legend=F)
dev.off()
		
pdf("Valerio_genes.pdf")		
par(mfrow=c(3,3))
plot.gene.trajectories("TP53I3",legend=F)
plot.gene.trajectories("TP73",legend=F)
plot.gene.trajectories("PARP14",legend=F)
		
plot.gene.trajectories("E2F4",legend=F)
plot.gene.trajectories("TRIM22",legend=F)
plot.gene.trajectories("RBL2",legend=F)
		
plot.gene.trajectories("TRIM4",legend=F)
plot.gene.trajectories("CSNK1E",legend=F)
plot.gene.trajectories("S100A14",legend=F)
dev.off()
		
pdf("Alu_neighbour_genes.pdf")	
par(mfrow=c(3,3))

plot.gene.trajectories("GNG12",legend=F)
plot.gene.trajectories("PI4K2A",legend=F)
plot.gene.trajectories("RHOD",legend=F)
		
plot.gene.trajectories("CAPN2",legend=F)
plot.gene.trajectories("RAP1A",legend=F)
plot.gene.trajectories("CHD9",legend=F)
		
plot.gene.trajectories("INO80C",legend=F)
plot.gene.trajectories("REV1",legend=F)
plot.gene.trajectories("GSN",legend=F)
#		
plot.gene.trajectories("PTK2",legend=F)
plot.gene.trajectories("ACTN4",legend=F)
plot.gene.trajectories("KATNAL2",legend=F)
		
plot.gene.trajectories("BCL2A1",legend=F)
plot.gene.trajectories("INCENP",legend=F)
plot.gene.trajectories("TMC6",legend=F)
		
plot.gene.trajectories("TMC8",legend=F)
plot.gene.trajectories("SUV39H1",legend=F)
plot.gene.trajectories("CBX5",legend=F)	
dev.off()

#		
#########################################
# cluster 2
wt0<-rowMeans(data[rownames(data) %in% names(which(cl_wt$cluster==2)),c(1,13)])
wt24<-rowMeans(data[rownames(data) %in% names(which(cl_wt$cluster==2)),c(2,14)])
wt72<-rowMeans(data[rownames(data) %in% names(which(cl_wt$cluster==2)),c(4,16)])

data_anova = data[rownames(data) %in% names(which(cl_wt$cluster==2)),]
		
table((wt0-wt24)>1 & (wt72-wt24)>1 )
		
interesting_genes = rownames(data_anova[(wt0-wt24)>1 & (wt72-wt24)>1,])
	
pdf("cluster2_2FC_filtered_genes.pdf")
par(mfrow=c(3,3))
for( i in 1:length(interesting_genes) ){
	plot.gene.trajectories(interesting_genes[i],legend=F)
	}
dev.off()
		
write.table(interesting_genes,file="cluster2_2FC_filtered_genes.txt",row.names=F,col.names=F,quote=F,sep="\t")

#########################################
# cluster 8
wt0<-rowMeans(data[rownames(data) %in% names(which(cl_wt$cluster==8)),c(1,13)])
wt24<-rowMeans(data[rownames(data) %in% names(which(cl_wt$cluster==8)),c(2,14)])
wt48<-rowMeans(data[rownames(data) %in% names(which(cl_wt$cluster==8)),c(3,15)])

data_anova = data[rownames(data) %in% names(which(cl_wt$cluster==8)),]
		
table((wt0-wt24)>1 & (wt48-wt24)>1 )
		
interesting_genes = rownames(data_anova[(wt0-wt24)>1 & (wt48-wt24)>1,])
	
pdf("cluster8_2FC_filtered_genes.pdf")
par(mfrow=c(3,3))
for( i in 1:length(interesting_genes) ){
	plot.gene.trajectories(interesting_genes[i],legend=F)
	}
dev.off()
		
write.table(interesting_genes,file="cluster8_2FC_filtered_genes.txt",row.names=F,col.names=F,quote=F,sep="\t")
##############################################################################################
# NO CLUSTERING
wt0<-rowMeans(data[rownames(data) %in% names(cl_wt$cluster),c(1,13)])
wt24<-rowMeans(data[rownames(data) %in% names(cl_wt$cluster),c(2,14)])
wt48<-rowMeans(data[rownames(data) %in% names(cl_wt$cluster),c(3,15)])
wt72<-rowMeans(data[rownames(data) %in% names(cl_wt$cluster),c(4,16)])

data_anova = data[rownames(data) %in% names(cl_wt$cluster),]
		
table((wt0-wt24)>1 & ((wt48-wt24)>1 | (wt72-wt24)>1 ) )
		
interesting_genes = rownames(data_anova[(wt0-wt24)>1 & ((wt48-wt24)>1 | (wt72-wt24)>1 ),])
	
pdf("NoClustering_2FC_filtered_genes.pdf")
par(mfrow=c(3,3))
for( i in 1:length(interesting_genes) ){
	plot.gene.trajectories(interesting_genes[i],legend=F)
	}
dev.off()
		
write.table(interesting_genes,file="NoClustering_2FC_filtered_genes.txt",row.names=F,col.names=F,quote=F,sep="\t")
		
pdf("Known_p53_targets.pdf")
par(mfrow=c(3,3))
plot.gene.trajectories("CDKN1A",legend=F)
plot.gene.trajectories("BBC3",legend=F)
plot.gene.trajectories("C12ORF5",legend=F)
		
plot.gene.trajectories("BTG2",legend=F)
plot.gene.trajectories("MDM2",legend=F)
plot.gene.trajectories("FAS",legend=F)
		
plot.gene.trajectories("PLK3",legend=F)
plot.gene.trajectories("PLK2",legend=F)
plot.gene.trajectories("SESN1",legend=F)
dev.off()

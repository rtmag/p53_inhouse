options(bitmapType="cairo")
options(scipen=999)
library(gplots)
library(ggplot2)

data<-read.table("Illu-Quant-expression.txt",row.names=1,header=T,sep="\t")

colnames(data)<-gsub("AVG_Signal.HCT116.","",colnames(data))
colnames(data)<-gsub("p007","WT",colnames(data))
colnames(data)<-gsub("\\.\\.",".",colnames(data),perl=T)

plot.gene.trajectories.points<-function(n,legend=FALSE){
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

plot(apply(wt1,2,mean),main=paste(c.name,""),ylim=c(c.min,c.max),xaxt='n',xlab="",ylab="")
points(apply(wt1,2,mean),col="red",cex=2.6,pch=19)
points(apply(wt2,2,mean),col="red",cex=2.6,pch=19)
	
points(apply(dm1,2,mean),col="blue",cex=2.6,pch=19)
points(apply(dm2,2,mean),col="blue",cex=2.6,pch=19)

points(apply(tp1,2,mean),col="green",cex=2.6,pch=19)
points(apply(tp2,2,mean),col="green",cex=2.6,pch=19)
Axis(side=1, labels=c('0h','24h','48h','72h'),at=c(1,2,3,4))
if(legend==TRUE){
	legend("topright", legend=c("Wt","TP53KO","DNMT1KO"), fill=c('salmon','olivedrab1','darkblue'),cex=.6, bty = "n")}
}
pdf("birc5_points.pdf")
plot.gene.trajectories.points("BIRC5")
dev.off()

barmat = rbind(as.numeric(wt1),as.numeric(wt2),as.numeric(tp1),as.numeric(tp2),as.numeric(dm1),as.numeric(dm2))
rownames(barmat)<-c("WT_1","WT_2","TP53KO_1","TP53KO_2","DNMT1_1","DNMT1_2")
colnames(barmat)<-c("0h","24h","48h","72h")
pdf("birc5_barplot.pdf")
barplot(barmat,beside=T,col=c("darkred","darkred","darkgreen","darkgreen","darkblue","darkblue"),ylim=c(0,13),main="BIRC5")
dev.off()


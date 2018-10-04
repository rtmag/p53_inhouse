
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

options(bitmapType="cairo")
options(scipen=999)
data<-read.table("Illu-Quant.txt",row.names=1,header=T,sep="\t")

pdf("gene_OMNI.pdf")
par(mfrow=c(3,3))
plot.gene.trajectories("CCNB1",legend=F)
plot.gene.trajectories("BRCA1",legend=F)
plot.gene.trajectories("BRCA2",legend=F)

plot.gene.trajectories("RAD51",legend=F)
plot.gene.trajectories("RAD54L",legend=F)
plot.gene.trajectories("BIRC5",legend=F)

plot.gene.trajectories("PLK1",legend=F)
plot.gene.trajectories("CDC6",legend=F)
plot.gene.trajectories("FAS",legend=F)
dev.off()


pdf("general_sns2.pdf")
par(mfrow=c(3,3))
plot.gene.trajectories("KAT5",legend=T)
plot.gene.trajectories("SUV39H1",legend=F)
plot.gene.trajectories("SUV39H2",legend=F)

plot.gene.trajectories("EHMT1",legend=F)
plot.gene.trajectories("EHMT2",legend=F)
plot.gene.trajectories("TRIM28",legend=F)

plot.gene.trajectories("CBX5",legend=F)
plot.gene.trajectories("IL6",legend=F)
plot.gene.trajectories("IL8",legend=F)
#
plot.gene.trajectories("ATM",legend=T)
plot.gene.trajectories("ATR",legend=F)
plot.gene.trajectories("TP53",legend=F)

plot.gene.trajectories("CDKN1A",legend=F)
plot.gene.trajectories("MDM2",legend=F)
plot.gene.trajectories("CDKN2A",legend=F)

plot.gene.trajectories("RB1",legend=F)
plot.gene.trajectories("E2F7",legend=F)
plot.gene.trajectories("FOXM1",legend=F)
#
plot.gene.trajectories("PTEN",legend=T)
plot.gene.trajectories("AKT1",legend=F)
plot.gene.trajectories("AKT2",legend=F)

plot.gene.trajectories("KRAS",legend=F)
plot.gene.trajectories("FRAP1",legend=F)
plot.gene.trajectories("E2F3",legend=F)

plot.gene.trajectories("DNMT1",legend=F)
plot.gene.trajectories("DNMT3A",legend=F)
plot.gene.trajectories("DNMT3B",legend=F)
#
plot.gene.trajectories("MAPK14",legend=T)
plot.gene.trajectories("BCL2",legend=F)
plot.gene.trajectories("BAX",legend=F)

plot.gene.trajectories("BAK1",legend=F)
plot.gene.trajectories("BBC3",legend=F)
plot.gene.trajectories("LMNB1",legend=F)

plot.gene.trajectories("LMNB2",legend=F)
plot.gene.trajectories("CDK2",legend=F)
plot.gene.trajectories("CDK4",legend=F)
dev.off()

# Senescence and apoptosis: dueling or complementary cell fates? EMBO review 2014 Notes

-ATM recognizes stress and phosphorilates p53.
-p14Arf inhibits MDM2, while MDM2 ubiquitinates p53 for degradation.
-ATM might supress ARF in some cases, high level regulation of this interaction is not understood.
-p53 upregulates p21, which in turn triggers initial cell cycle arrest to give the cell time to repair the DNA (before S-phase entry).
-Prolonged arrest activates the CDKi p16Ink4a, which in turns activates Rb. (prolonged p16Ink4a results in permanent cell cycle arrest.




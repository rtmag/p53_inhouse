cut -f1,2 wt_48_methylated.bed|awk -F"\t" '{print $1"\t"$2"\t"$2+1}' | \
annotatePeaks.pl - hg38 -annStats wt_48_methylated.annStats > wt_48_methylated.anno &

cut -f1,2 p53_48_methylated.bed|awk -F"\t" '{print $1"\t"$2"\t"$2+1}' | \
annotatePeaks.pl - hg38 -annStats p53_48_methylated.annStats > p53_48_methylated.anno &

############################################################################################################
par(mfrow=c(1,1))
pdf("distanno_WT48.pdf")
par(mar=c(11.1,4.1,4.1,2))
res=read.table(pipe("more cpg_methylated_in_WT_48h.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>360]
barplot(sort(tdown),las=2,ylab="CpGs",col="darkred",ylim=c(0,60000))
title("CpG with high methylation in WT 48h", cex.main=.9)
dev.off()

pdf("distanno_TP5348h.pdf")
par(mar=c(11.1,4.1,4.1,2))
res=read.table(pipe("more cpg_methylated_in_TP53KO_48.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>360]
barplot(sort(tdown),las=2,ylab="CpGs",col="lightblue3",ylim=c(0,20000))
title("CpG with high methylation in TP53KO 48h", cex.main=.9)
dev.off()


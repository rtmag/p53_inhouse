cut -f1,2 wt_48_methylated.bed|awk -F"\t" '{print $1"\t"$2"\t"$2+1}' | \
annotatePeaks.pl - hg38 -annStats wt_48_methylated.annStats > wt_48_methylated.anno &

cut -f1,2 p53_48_methylated.bed|awk -F"\t" '{print $1"\t"$2"\t"$2+1}' | \
annotatePeaks.pl - hg38 -annStats p53_48_methylated.annStats > p53_48_methylated.anno &

############################################################################################################

res=read.table(pipe("more p53_48_methylated.annStats |cut -f1,2,4"), sep="\t",header=F)
i1 = which(res[,1]=="Annotation")[2]+1
i2 = dim(res)[1]
res = res[ i1:i2,]
tdown = as.numeric(as.character(res[,2]))
names(tdown) = res[,1]
names(tdown) = paste(names(tdown)," ",round(tdown/sum(tdown)*100,digits=2),"%",sep="")
tdown = tdown[tdown>200]
pie(sort(tdown), main=,cex=.8,col=NA)
title("Distribution of chromatin regions open in SH compared to NT\n(15,360 regions)", cex.main=.9)

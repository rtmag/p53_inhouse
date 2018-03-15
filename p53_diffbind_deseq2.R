library(Rsubread)

x=read.table('/root/stuff/macs2/P53_24h_48h_doxo_merged.bed',sep="\t",stringsAsFactors=F)
ann = data.frame(GeneID=paste(x[,1],x[,2],x[,3],sep="_!_"),Chr=x[,1],Start=x[,2],End=x[,3],Strand='+')
bam.files <- c('/root/stuff/bam/P53_24h_doxo_s1_rmdup.bam',
              '/root/stuff/bam/P53_48h_doxo_s1_rmdup.bam')

fc_SE <- featureCounts(bam.files,annot.ext=ann,isPairedEnd=TRUE,nthreads=20)
countData=fc_SE$counts

colnames(countData)=c("P53_24h_doxo","P53_48h_doxo")
saveRDS(countData,'p53_doxo_countdata.rds')
#####
library(DESeq2)
countData=readRDS('p53_doxo_countdata.rds')

dds <- DESeqDataSetFromMatrix(
       countData = countData,
       colData = data.frame(group=c("24h","48h")),
       design = ~ group)
rld <- rlogTransformation( dds )
res <- data.frame(rLogFC = assay(rld)[,2] - assay(rld)[,1],
                  p53_24h=countData[,1],p53_48h=countData[,2])

vst = varianceStabilizingTransformation(dds)
res_vst <- data.frame(vstFC = assay(vst)[,2] - assay(vst)[,1],
                  p53_24h=countData[,1],p53_48h=countData[,2])

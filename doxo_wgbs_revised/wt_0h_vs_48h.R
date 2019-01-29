
library(methylKit)

file.list=list( "/root/hct116_wgbs_revised/combine_hct116_WGBS_1_val_1_bismark_bt2_pe.CpG_report.txt.gz",
               "/root/hct116_wgbs_revised/P007_48h_doxo_1_val_1_bismark_bt2_pe.CpG_report.txt.gz")

# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
           sample.id=list("WT_0h","WT_48h"),
           assembly="hg38",
           treatment=c(0,1),
           context="CpG",
           pipeline="bismarkCytosineReport",
           header=FALSE,
           mincov=3)

meth_0=unite(myobj, destrand=TRUE,mc.cores=40)
pooled.meth_0=pool(meth_0,sample.ids=c("WT_0h","WT_48h"))
pooled.myDiff_0=calculateDiffMeth(pooled.meth_0,num.cores=40)
pooledData_0=getData(pooled.meth_0)

wt0=fisher_0[fisher_0$qvalue<0.1 & abs(fisher_0$meth.diff)>25,]
wt0_l=fisher_0[fisher_0$qvalue<0.1 & (fisher_0$meth.diff)>(-25),]
wt48_l=fisher_0[fisher_0$qvalue<0.1 & (fisher_0$meth.diff)>25,]

write.table(wt0,"wt_0_48_methylated.bed",quote=F,col.names=F,row.names=F,sep="\t")
write.table(wt0_l,"wt_0_methylated.bed",quote=F,col.names=F,row.names=F,sep="\t")
write.table(wt48_l,"wt_48_methylated.bed",quote=F,col.names=F,row.names=F,sep="\t")


##


file.list=list( "/root/wgbs_doxo/OriginalFastq/P007_48h_doxo_1_bismark_bt2_pe.CpG_report.txt",                
                "/root/wgbs_doxo/OriginalFastq/TP53del_48h_doxo_1_bismark_bt2_pe.CpG_report.txt")

# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
           sample.id=list("WT_48h","P53KO_48h"),
           assembly="hg38",
           treatment=c(0,1),
           context="CpG",
           pipeline="bismarkCytosineReport",
           header=FALSE,
           mincov=5)

##

meth_48=unite(myobj, destrand=TRUE,mc.cores=40)
pooled.meth_48=pool(meth_48,sample.ids=c("WT_48h","P53KO_48h"))
pooled.myDiff_48=calculateDiffMeth(pooled.meth_48,num.cores=40)
pooledData_48=getData(pooled.meth_48)
####
fisher_0=data.frame(pooled.myDiff_0,WT_0h=(pooledData_0[,6]/pooledData_0[,5]),WT_48h=(pooledData_0[,9]/pooledData_0[,8]) ) 
write.csv(fisher_0,"fisher_0.csv")
saveRDS(fisher_0,"fisher_0.rds")

fisher_48=data.frame(pooled.myDiff_48,WT_48h=(pooledData_48[,6]/pooledData_48[,5]),P53KO_48h=(pooledData_48[,9]/pooledData_48[,8]) ) 
write.csv(fisher_48,"fisher_48.csv")
saveRDS(fisher_48,"fisher_48.rds")
####
dim(fisher_0[(fisher_0$qvalue<0.05 & abs(fisher_0$meth.diff)>35),])
dim(fisher_48[(fisher_48$qvalue<0.05 & abs(fisher_48$meth.diff)>50),])


wt48 = fisher_48[fisher_48$qvalue<0.05 & fisher_48$meth.diff>50,]
p5348= fisher_48[fisher_48$qvalue<0.05 & fisher_48$meth.diff<(-50),]

wt48 = wt48[as.numeric(wt48[,1])<25,]
p5348 = p5348[as.numeric(p5348[,1])<25,]

wt48[,1] = gsub("^","chr",wt48[,1],perl=T)
wt48[,1] = gsub("chr23","chrX",wt48[,1],perl=T)
wt48[,1] = gsub("chr24","chrY",wt48[,1],perl=T)

p5348[,1] = gsub("^","chr",p5348[,1],perl=T)
p5348[,1] = gsub("chr23","chrX",p5348[,1],perl=T)
p5348[,1] = gsub("chr24","chrY",p5348[,1],perl=T)

write.table(wt48,"wt_48_methylated.bed",quote=F,col.names=F,row.names=F,sep="\t")
write.table(p5348,"p53_48_methylated.bed",quote=F,col.names=F,row.names=F,sep="\t")

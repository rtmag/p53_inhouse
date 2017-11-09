library(methylKit)

file.list=list( "P007_48h_doxo_1_bismark_bt2_pe.CpG_report.txt",                
                "TP53del_48h_doxo_1_bismark_bt2_pe.CpG_report.txt")

# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
           sample.id=list("WT","P53ko"),
           assembly="hg38",
           treatment=c(1,0),
           context="CpG",
           pipeline="bismarkCytosineReport",
           header=FALSE,
           mincov=10)

##

meth=unite(myobj, destrand=TRUE,mc.cores=40)
pooled.meth=pool(meth,sample.ids=c("WT","P53ko"))
pooled.myDiff=calculateDiffMeth(pooled.meth,num.cores=40)

pooledData=getData(pooled.meth)
#
fisher=data.frame(pooled.myDiff,WT=(pooledData[,6]/pooledData[,5]),P53ko=(pooledData[,9]/pooledData[,8]) ) 

diff=fisher[(fisher$qvalue<0.05 & abs(fisher$meth.diff)>25),]

diff=diff[as.numeric(diff[,1])<26,]
chrNames=gsub("^","chr",diff[,1],perl=T)
chrNames=gsub("chr23","chrX",chrNames,perl=T)
chrNames=gsub("chr24","chrY",chrNames,perl=T)
chrNames=gsub("chr25","chrM",chrNames,perl=T)
diff[,1]=chrNames
diff[,3]=diff[,3]+1

saveRDS(diff,"diff_wt_tp53.rds")

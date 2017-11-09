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
           mincov=3)

##

meth=unite(myobj, destrand=TRUE,mc.cores=40)
pooled.meth=pool(meth,sample.ids=c("WT","P53ko"))
pooled.myDiff=calculateDiffMeth(pooled.meth,num.cores=40)

pooledData=getData(pooled.meth)
#

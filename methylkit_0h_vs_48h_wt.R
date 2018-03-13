
library(methylKit)

file.list=list( "/root/HCT116_wgbs/combine_hct116_WGBS_1_val_1_bismark_bt2_pe.CpG_report.txt",
               "/root/wgbs_doxo/OriginalFastq/P007_48h_doxo_1_bismark_bt2_pe.CpG_report.txt")

# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
           sample.id=list("WT_0h","WT_48h"),
           assembly="hg38",
           treatment=c(0,1),
           context="CpG",
           pipeline="bismarkCytosineReport",
           header=FALSE,
           mincov=5)

meth_0=unite(myobj, destrand=TRUE,mc.cores=40)
pooled.meth_0=pool(meth_0,sample.ids=c("WT_0h","WT_48h"))
pooled.myDiff_0=calculateDiffMeth(pooled.meth_0,num.cores=40)
pooledData_0=getData(pooled.meth_0)
####

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
           mincov=10)

##

meth_48=unite(myobj, destrand=TRUE,mc.cores=40)
pooled.meth_48=pool(meth_48,sample.ids=c("WT_48h","P53KO_48h"))
pooled.myDiff_48=calculateDiffMeth(pooled.meth_48,num.cores=40)
pooledData_48=getData(pooled.meth_48)
####

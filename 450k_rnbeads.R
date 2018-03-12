library(RnBeads.hg38)
library(RnBeads)
rnb.set.norm=load.rnb.set("rnb.set.norm.RData.zip")


rnb.set.norm@pheno=data.frame(rnb.set.norm@pheno, 
                              Time=c("0","12","24","48","0","12","24","48","0","12","48","ctrl",
                                     "0","12","24","48","0","12","24","48","0","12","48","ctrl"),
                              reps=c("wt_0","wt_12","wt_24","wt_48","p53_0","p53_12","p53_24","p53_48",
                                    "dnmt1_0","dnmt1_12","dnmt1_48","ctrl",
                                    "wt_0","wt_12","wt_24","wt_48","p53_0","p53_12","p53_24","p53_48",
                                    "dnmt1_0","dnmt1_12","dnmt1_48","ctrl")
                             )

rnb.set.norm@assembly = "hg38"

select = c("wt_0","wt_12","wt_24","wt_48","p53_0","p53_12","p53_24","p53_48",
                                    "dnmt1_0","dnmt1_12","dnmt1_48","ctrl",
                                    "wt_0","wt_12","wt_24","wt_48","p53_0","p53_12","p53_24","p53_48",
                                    "dnmt1_0","dnmt1_12","dnmt1_48","ctrl")

WT24h_VS_WT48h=remove.samples(rnb.set.norm,samples(rnb.set.norm)[which(!select %in% c("wt_24","wt_48"))] )

WT24h_VS_WT48h_dmc <- rnb.execute.computeDiffMeth(WT24h_VS_WT48h,pheno.cols=c("reps"))


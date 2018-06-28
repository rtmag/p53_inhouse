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

select = c("wt_0","wt_12","wt_24","wt_48","p53_0","p53_12","p53_24","p53_48",
                                    "dnmt1_0","dnmt1_12","dnmt1_48","ctrl",
                                    "wt_0","wt_12","wt_24","wt_48","p53_0","p53_12","p53_24","p53_48",
                                    "dnmt1_0","dnmt1_12","dnmt1_48","ctrl")

#########################################
WT0h_VS_WT48h=remove.samples(rnb.set.norm,samples(rnb.set.norm)[which(!select %in% c("wt_0","wt_48"))] )

meth = meth(WT0h_VS_WT48h)

WT0h_VS_WT48h_dmc <- rnb.execute.computeDiffMeth(WT0h_VS_WT48h,pheno.cols=c("reps"))
comparison <- get.comparisons(WT0h_VS_WT48h_dmc)[1]
WT0h_VS_WT48h_dmc_table <-get.table(WT0h_VS_WT48h_dmc, comparison, "sites", return.data.frame=TRUE)
#

#
saveRDS(WT0h_VS_WT48h_dmc_table,"WT0h_VS_WT48h_dmc_table.rds")

idx=which(abs(WT0h_VS_WT48h_dmc_table$mean.diff)>.35 & WT0h_VS_WT48h_dmc_table$diffmeth.p.adj.fdr<0.05)
idx=which(abs(WT0h_VS_WT48h_dmc_table$mean.diff)>.25)

#########################################
WT24h_VS_WT48h=remove.samples(rnb.set.norm,samples(rnb.set.norm)[which(!select %in% c("wt_24","wt_48"))] )
WT24h_VS_WT48h_dmc <- rnb.execute.computeDiffMeth(WT24h_VS_WT48h,pheno.cols=c("reps"))

comparison <- get.comparisons(WT24h_VS_WT48h_dmc)[1]
WT24h_VS_WT48h_dmc_table <-get.table(WT24h_VS_WT48h_dmc, comparison, "sites", return.data.frame=TRUE)

saveRDS(WT24h_VS_WT48h_dmc_table,"WT24h_VS_WT48h_dmc_table.rds")

idx=which(abs(WT24h_VS_WT48h_dmc_table$mean.diff)>.35 & WT24h_VS_WT48h_dmc_table$diffmeth.p.adj.fdr<0.05)
#########################################
WT_VS_TP53_0h=remove.samples(rnb.set.norm,samples(rnb.set.norm)[which(!select %in% c("wt_0","p53_0"))] )
WT_VS_TP53_0h_dmc <- rnb.execute.computeDiffMeth(WT_VS_TP53_0h,pheno.cols=c("reps"))


comparison <- get.comparisons(WT_VS_TP53_0h_dmc)[1]
WT_VS_TP53_0h_dmc_dmc_table <-get.table(WT_VS_TP53_0h_dmc, comparison, "sites", return.data.frame=TRUE)

saveRDS(WT_VS_TP53_0h_dmc_table,"WT_VS_TP53_0h_dmc_table.rds")

idx=which(abs(WT_VS_TP53_0h_dmc_table$mean.diff)>.35 & WT_VS_TP53_0h_dmc_table$diffmeth.p.adj.fdr<0.05)
#########################################
#p53_24 p53_48
TP5324h_VS_TP5348h=remove.samples(rnb.set.norm,samples(rnb.set.norm)[which(!select %in% c("p53_24","p53_48"))] )
TP5324h_VS_TP5348h_dmc <- rnb.execute.computeDiffMeth(TP5324h_VS_TP5348h,pheno.cols=c("reps"))

comparison <- get.comparisons(TP5324h_VS_TP5348h_dmc)[1]
TP5324h_VS_TP5348h_dmc_table <-get.table(TP5324h_VS_TP5348h_dmc, comparison, "sites", return.data.frame=TRUE)

saveRDS(TP5324h_VS_TP5348h_dmc_table,"TP5324h_VS_TP5348h_dmc_table.rds")

idx=which(abs(TP5324h_VS_TP5348h_dmc_table$mean.diff)>.35 & TP5324h_VS_TP5348h_dmc_table$diffmeth.p.adj.fdr<0.05)

#########################################
WT_VS_TP53_48h=remove.samples(rnb.set.norm,samples(rnb.set.norm)[which(!select %in% c("wt_48","p53_48"))] )
WT_VS_TP53_48h_dmc <- rnb.execute.computeDiffMeth(WT_VS_TP53_48h,pheno.cols=c("reps"))


comparison <- get.comparisons(WT_VS_TP53_48h_dmc)[1]
WT_VS_TP53_48h_dmc_dmc_table <-get.table(WT_VS_TP53_48h_dmc, comparison, "sites", return.data.frame=TRUE)

saveRDS(WT_VS_TP53_48h_dmc_dmc_table,"WT_VS_TP53_48h_dmc_table.rds")

idx=which(abs(WT_VS_TP53_48h_dmc_dmc_table$mean.diff)>.35 & WT_VS_TP53_48h_dmc_dmc_table$diffmeth.p.adj.fdr<0.05)
#############
# p53 0h vs p53 48
TP530h_VS_TP5348h=remove.samples(rnb.set.norm,samples(rnb.set.norm)[which(!select %in% c("p53_0","p53_48"))] )
TP530h_VS_TP5348h_dmc <- rnb.execute.computeDiffMeth(TP530h_VS_TP5348h,pheno.cols=c("reps"))

comparison <- get.comparisons(TP530h_VS_TP5348h_dmc)[1]
TP530h_VS_TP5348h_dmc_table <-get.table(TP530h_VS_TP5348h_dmc, comparison, "sites", return.data.frame=TRUE)

saveRDS(TP530h_VS_TP5348h_dmc_table,"TP530h_VS_TP5348h_dmc_table.rds")

idx=which(abs(TP530h_VS_TP5348h_dmc_table$mean.diff)>.35 & TP530h_VS_TP5348h_dmc_table$diffmeth.p.adj.fdr<0.05)
#############
# wt 0h vs p53 0
WT0h_VS_p530h=remove.samples(rnb.set.norm,samples(rnb.set.norm)[which(!select %in% c("wt_0","p53_0"))] )
WT0h_VS_p530h_dmc <- rnb.execute.computeDiffMeth(WT0h_VS_p530h,pheno.cols=c("reps"))

comparison <- get.comparisons(WT0h_VS_p530h_dmc)[1]
WT0h_VS_p530h_dmc_table <-get.table(WT0h_VS_p530h_dmc, comparison, "sites", return.data.frame=TRUE)

idx=which(abs(WT0h_VS_p530h_dmc_table$mean.diff)>.35 & WT0h_VS_p530h_dmc_table$diffmeth.p.adj.fdr<0.05)

saveRDS(WT0h_VS_p530h,"WT0h_VS_p530h_dmc_table.rds")
# wt 48h vs p53 48
WT48h_VS_p5348h=remove.samples(rnb.set.norm,samples(rnb.set.norm)[which(!select %in% c("wt_48","p53_48"))] )
WT48h_VS_p5348h_dmc <- rnb.execute.computeDiffMeth(WT48h_VS_p5348h,pheno.cols=c("reps"))

comparison <- get.comparisons(WT48h_VS_p5348h_dmc)[1]
WT48h_VS_p5348h_dmc_table <-get.table(WT48h_VS_p5348h_dmc, comparison, "sites", return.data.frame=TRUE)

idx=which(abs(WT48h_VS_p5348h_dmc_table$mean.diff)>.35 & WT48h_VS_p5348h_dmc_table$diffmeth.p.adj.fdr<0.05)
smoothScatter(WT48h_VS_p5348h_dmc_table$mean.diff[idx], -log10(WT48h_VS_p5348h_dmc_table$diffmeth.p.adj.fdr[idx]))
saveRDS(WT0h_VS_p530h,"WT0h_VS_p530h_dmc_table.rds")

########
par(mfrow = c(1,2))
smoothScatter(WT0h_VS_p530h_dmc_table$mean.diff, -log10(WT0h_VS_p530h_dmc_table$diffmeth.p.adj.fdr),ylab = "-log10 FDR",xlab = "WT 0h VS TP53KO 0h")
smoothScatter(WT48h_VS_p5348h_dmc_table$mean.diff, -log10(WT48h_VS_p5348h_dmc_table$diffmeth.p.adj.fdr),ylab = "-log10 FDR",xlab = "WT 48h VS TP53KO 48h")

idx0=which(abs(WT0h_VS_p530h_dmc_table$mean.diff)>.35 & WT0h_VS_p530h_dmc_table$diffmeth.p.adj.fdr<0.05)
idx48=which(abs(WT48h_VS_p5348h_dmc_table$mean.diff)>.35 & WT48h_VS_p5348h_dmc_table$diffmeth.p.adj.fdr<0.05)

smoothScatter(WT48h_VS_p5348h_dmc_table$mean.diff[idx], -log10(WT48h_VS_p5348h_dmc_table$diffmeth.p.adj.fdr[idx]))

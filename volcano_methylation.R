library(graphics)

png("WT0h_VS_WT48h_volcano.png")
x=readRDS("WT0h_VS_WT48h_dmc_table.rds")
smoothScatter(x$mean.diff,-log10(x$diffmeth.p.adj.fdr),xlim=c(-1,1),ylim=c(0,1),
             xlab = "DNA methylation difference (WT 0h / WT 48h)",
             ylab = "-log10 FDR")
dev.off()

x=readRDS("WT24h_VS_WT48h_dmc_table.rds")
png("WT24h_VS_WT48h_volcano.png")
smoothScatter(x$mean.diff,-log10(x$diffmeth.p.adj.fdr),xlim=c(-1,1),
             xlab = "DNA methylation difference (WT 24h / WT 48h)",
             ylab = "-log10 FDR")
dev.off()

x=readRDS("WT_VS_TP53_0h_dmc_table.rds")
png("WT0h_VS_TP530h_volcano.png")
smoothScatter(x$mean.diff,-log10(x$diffmeth.p.adj.fdr),xlim=c(-1,1),
             xlab = "DNA methylation difference (WT 0h / TP53KO 0h)",
             ylab = "-log10 FDR")
dev.off()

x=readRDS("WT_VS_TP53_48h_dmc_table.rds")
png("WT48h_VS_TP5348h_volcano.png")
smoothScatter(x$mean.diff,-log10(x$diffmeth.p.adj.fdr),xlim=c(-1,1),
             xlab = "DNA methylation difference (WT 48h / TP53KO 48h)",
             ylab = "-log10 FDR")
dev.off()


x=readRDS("TP5324h_VS_TP5348h_dmc_table.rds")
png("TP5324h_VS_TP5348h_volcano.png")
smoothScatter(x$mean.diff,-log10(x$diffmeth.p.adj.fdr),xlim=c(-1,1),ylim=c(0,3.8),
             xlab = "DNA methylation difference (TP53KO 24h / TP53KO 48h)",
             ylab = "-log10 FDR")
dev.off()

##
fisher_0=readRDS("fisher_48.rds")
library(graphics)

png("WT48h_VS_WT48h_volcano_WGBS.png")
smoothScatter(fisher_0$meth.diff,-log10(fisher_0$qvalue),xlim=c(-100,100),ylim=c(0,5),
             xlab = "DNA methylation difference (WT 0h / WT 48h)",
             ylab = "-log10 FDR",nrpoints=0)
abline(h=-log10(.05),lty=2)
abline(v=50,lty=2)
abline(v=-50,lty=2)
legend("topright", paste("WT 48h:",length(which(fisher_0$meth.diff>50 & fisher_0$qvalue<0.05))), bty="n") 
legend("topleft", paste("TP53 48h:",length(which(fisher_0$meth.diff<(-50) & fisher_0$qvalue<0.05))), bty="n") 
points(fisher_0$meth.diff[abs(fisher_0$meth.diff)>50 & fisher_0$qvalue<0.05],
       -log10(fisher_0$qvalue[abs(fisher_0$meth.diff)>50 & fisher_0$qvalue<0.05]),
       col=alpha("#c0392b",.05)
dev.off()






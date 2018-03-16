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

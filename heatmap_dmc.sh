cut -f1,2 wt_48_methylated.bed|awk -F"\t" '{print $1"\t"$2"\t"$2+1}' > dmc_48h.bed
echo "# WT 48h Methylated" >> dmc_48h.bed
cut -f1,2,3 p53_48_methylated.bed|awk -F"\t" '{print $1"\t"$2"\t"$2+1}' >> dmc_48h.bed
echo "# TP53KO 48h Methylated" >> dmc_48h.bed

computeMatrix reference-point \
-S /root/p53_dnmt1/bw/P53_24h_doxo_s1.bw \
/root/p53_dnmt1/bw/P53_48h_doxo_s1.bw \
-R /root/HCT116_wgbs/methylkit/dmc_48h.bed --referencePoint center \
--sortRegions descend -bs 100 -a 2000 -b 2000 -p max \
-out /root/HCT116_wgbs/methylkit/dmc_p53_48h.mat

plotHeatmap --refPointLabel "CpG" -m /root/HCT116_wgbs/methylkit/dmc_p53_48h.mat \
--colorMap Blues -out /root/HCT116_wgbs/methylkit/dmc_p53_48h.pdf
##
####
##

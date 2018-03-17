##
cut -f1,2 wt_48_methylated.bed|awk -F"\t" '{print $1"\t"$2"\t"$2+1}' > cpg_methylated_in_TP53KO_48.bed
cut -f1,2 p53_48_methylated.bed|awk -F"\t" '{print $1"\t"$2"\t"$2+1}' > cpg_methylated_in_WT_48.bed

cut -f1,2 cpg_methylated_in_TP53KO_48.bed|awk -F"\t" '{print $1"\t"$2"\t"$2+1}' > dmc_48h.bed
echo "#CpG methylated in TP53KO 48h" >> dmc_48h.bed
cut -f1,2,3 cpg_methylated_in_WT_48.bed|awk -F"\t" '{print $1"\t"$2"\t"$2+1}' >> dmc_48h.bed
echo "#CpG methylated in WT 48h" >> dmc_48h.bed


##


computeMatrix reference-point \
-S /root/p53_dnmt1/bw/P53_24h_doxo_s1.bw \
/root/p53_dnmt1/bw/P53_48h_doxo_s1.bw \
-R /root/HCT116_wgbs/methylkit/dmc_48h.bed --referencePoint center \
--sortRegions descend -bs 100 -a 2000 -b 2000 -p max \
-out /root/HCT116_wgbs/methylkit/dmc_p53_48h.mat

plotHeatmap --refPointLabel "CpG" -m /root/HCT116_wgbs/methylkit/dmc_p53_48h.mat \
--zMax 50 -out /root/HCT116_wgbs/methylkit/dmc_p53_48h.pdf
##
####
##

computeMatrix reference-point \
-S /root/p53_dnmt1/bw/WT_0h_doxo_wgbs.bw \
/root/p53_dnmt1/bw/WT_48h_doxo_wgbs.bw \
/root/p53_dnmt1/bw/P53KO_48h_doxo_wgbs.bw \
-R /root/HCT116_wgbs/methylkit/dmc_48h.bed --referencePoint center \
--sortRegions descend -bs 100 -a 2000 -b 2000 -p max \
-out /root/p53_dnmt1/heatmap/dmc_48h_wgbs.mat

plotHeatmap --refPointLabel "CpG" -m /root/p53_dnmt1/heatmap/dmc_48h_wgbs.mat \
 -out /root/p53_dnmt1/heatmap/dmc_48h_wgbs_bicolor.pdf

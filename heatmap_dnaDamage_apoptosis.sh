cat c25_tss.bed > timecourse_genes.bed
echo "# DNA Damage" >> timecourse_genes.bed
cat c134_tss.bed >> timecourse_genes.bed
echo "# Apoptosis" >> timecourse_genes.bed

computeMatrix reference-point \
-S /root/p53_dnmt1/bw/WT_0h_doxo_wgbs.bw \
/root/p53_dnmt1/bw/WT_48h_doxo_wgbs.bw \
/root/p53_dnmt1/bw/P53KO_48h_doxo_wgbs.bw \
-R /root/p53_dnmt1/heatmap/timecourse_genes.bed --referencePoint center \
--sortRegions descend -bs 100 -a 4000 -b 4000 -p max \
-out /root/p53_dnmt1/heatmap/timecourse_genes_wgbs.mat

plotHeatmap --refPointLabel "TSS" -m /root/p53_dnmt1/heatmap/timecourse_genes_wgbs.mat \
black --colorMap Blues -out /root/p53_dnmt1/heatmap/timecourse_genes_wgbs.pdf
#
computeMatrix reference-point \
-S /root/p53_dnmt1/bw/P53_24h_doxo_s1.bw \
/root/p53_dnmt1/bw/P53_48h_doxo_s1.bw \
-R /root/p53_dnmt1/heatmap/timecourse_genes.bed --referencePoint center \
--sortRegions descend -bs 100 -a 2000 -b 2000 -p max \
-out /root/p53_dnmt1/heatmap/timecourse_genes_p53.mat

plotHeatmap --refPointLabel "TSS" -m /root/p53_dnmt1/heatmap/timecourse_genes_p53.mat \
--colorMap Blues -out /root/p53_dnmt1/heatmap/timecourse_genes_p53.pdf

############################################################################################################

cat c25_gene.bed > timecourse_genesbody.bed
echo "# DNA Damage" >> timecourse_genesbody.bed
cat c134_gene.bed >> timecourse_genesbody.bed
echo "# Apoptosis" >> timecourse_genesbody.bed


computeMatrix scale-regions \
-S /root/p53_dnmt1/bw/WT_0h_doxo_wgbs.bw \
/root/p53_dnmt1/bw/WT_48h_doxo_wgbs.bw \
/root/p53_dnmt1/bw/P53KO_48h_doxo_wgbs.bw \
-R /root/p53_dnmt1/heatmap/timecourse_genesbody.bed  \
--beforeRegionStartLength 3000 \
--regionBodyLength 5000 \
--afterRegionStartLength 3000 \
--sortRegions descend -bs 100 -p max \
-out /root/p53_dnmt1/heatmap/timecourse_genesbody_wgbs.mat

plotHeatmap -m /root/p53_dnmt1/heatmap/timecourse_genesbody_wgbs.mat \
--colorMap Blues -out /root/p53_dnmt1/heatmap/timecourse_genesbody_wgbs.pdf
#
computeMatrix scale-regions \
-S /root/p53_dnmt1/bw/P53_24h_doxo_s1.bw \
/root/p53_dnmt1/bw/P53_48h_doxo_s1.bw \
-R /root/p53_dnmt1/heatmap/timecourse_genesbody.bed  \
--beforeRegionStartLength 3000 \
--regionBodyLength 5000 \
--afterRegionStartLength 3000 \
--sortRegions descend -bs 100 -p max \
-out /root/p53_dnmt1/heatmap/timecourse_genesbody_p53.mat


plotHeatmap -m /root/p53_dnmt1/heatmap/timecourse_genesbody_p53.mat \
--colorMap Blues -out /root/p53_dnmt1/heatmap/timecourse_genesbody_p53.pdf
############################################################################################################

grep "Down" ../diffbind/diffreps/P53_48h_vs_24h_std_nsd20_modeN.bed > P53_diff.bed
echo "#24h Doxo" >> P53_diff.bed
grep "Up" ../diffbind/diffreps/P53_48h_vs_24h_std_nsd20_modeN.bed >> P53_diff.bed
echo "#48h Doxo" >> P53_diff.bed


computeMatrix reference-point \
-S /root/p53_dnmt1/bw/P53_24h_doxo_s1.bw \
/root/p53_dnmt1/bw/P53_48h_doxo_s1.bw \
-R P53_diff.bed --referencePoint center \
--sortRegions descend -bs 100 -a 2000 -b 2000 -p max \
-out P53_diff.mat

plotHeatmap -m P53_diff.mat --xAxisLabel "" --yAxisLabel "" --refPointLabel "P53 Peak" \
--colorMap Blues -out P53_diff.pdf

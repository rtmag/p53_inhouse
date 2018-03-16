bamToBed -i /root/stuff/bam/P53_24h_doxo_s1_rmdup.bam > /root/stuff/diffbind/bed/P53_24h_doxo.bed &
bamToBed -i /root/stuff/bam/P53_48h_doxo_s1_rmdup.bam > /root/stuff/diffbind/bed/P53_48h_doxo.bed &

wait

diffReps.pl --treatment /root/stuff/diffbind/bed/P53_48h_doxo.bed --control /root/stuff/diffbind/bed/P53_24h_doxo.bed \
--nsd 20 --mode n --meth gt --gname hg19 --report /root/stuff/diffbind/diffreps/P53_48h_vs_24h_nsd20_modeN.diffreps --frag 0 --nproc 30 

##

diffReps.pl --treatment /root/stuff/diffbind/bed/P53_48h_doxo.bed --control /root/stuff/diffbind/bed/P53_24h_doxo.bed \
--nsd 20 --std --mode n --meth gt --gname hg19 --report /root/stuff/diffbind/diffreps/P53_48h_vs_24h_std_nsd20_modeN.diffreps \
--frag 0 --nproc 30 

##
#
cat /root/stuff/diffbind/deseq2/p53_24h_higher_rlog.3.bed > control_p53_diffbind_rlog.3.bed
echo "# P53 24h Doxo" >> control_p53_diffbind_rlog.3.bed
cat /root/stuff/diffbind/deseq2/p53_48h_higher_rlog.3.bed >> control_p53_diffbind_rlog.3.bed
echo "# P53 48h Doxo" >> control_p53_diffbind_rlog.3.bed

cat /root/stuff/diffbind/deseq2/p53_24h_higher_rlog.4.bed > control_p53_diffbind_rlog.4.bed
echo "# P53 24h Doxo" >> control_p53_diffbind_rlog.4.bed
cat /root/stuff/diffbind/deseq2/p53_48h_higher_rlog.4.bed >> control_p53_diffbind_rlog.4.bed
echo "# P53 48h Doxo" >> control_p53_diffbind_rlog.4.bed

####
computeMatrix reference-point \
-S /root/p53_dnmt1/bw/P53_24h_doxo_s1.bw \
/root/p53_dnmt1/bw/P53_48h_doxo_s1.bw \
-R /root/stuff/diffbind/deseq2/control_p53_diffbind_rlog.3.bed --referencePoint center \
--sortRegions descend -bs 100 -a 2000 -b 2000 -p max \
-out /root/p53_dnmt1/heatmap/control_p53_diffbind_rlog.3.mat

plotHeatmap --refPointLabel "peak" -m /root/p53_dnmt1/heatmap/control_p53_diffbind_rlog.3.mat \
--colorMap Blues -out /root/p53_dnmt1/heatmap/control_p53_diffbind_rlog.3.pdf

computeMatrix reference-point \
-S /root/p53_dnmt1/bw/P53_24h_doxo_s1.bw \
/root/p53_dnmt1/bw/P53_48h_doxo_s1.bw \
-R /root/stuff/diffbind/deseq2/control_p53_diffbind_rlog.4.bed --referencePoint center \
--sortRegions descend -bs 100 -a 2000 -b 2000 -p max \
-out /root/p53_dnmt1/heatmap/control_p53_diffbind_rlog.4.mat

plotHeatmap --refPointLabel "peak" -m /root/p53_dnmt1/heatmap/control_p53_diffbind_rlog.4.mat \
--colorMap Blues -out /root/p53_dnmt1/heatmap/control_p53_diffbind_rlog.4.pdf
###
##
#

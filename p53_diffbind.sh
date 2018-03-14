bamToBed -i /root/stuff/bam/P53_24h_doxo_s1_rmdup.bam > /root/stuff/diffbind/bed/P53_24h_doxo.bed &
bamToBed -i /root/stuff/bam/P53_48h_doxo_s1_rmdup.bam > /root/stuff/diffbind/bed/P53_48h_doxo.bed &

wait

diffReps.pl --treatment /root/stuff/diffbind/bed/P53_48h_doxo.bed --control /root/stuff/diffbind/bed/P53_24h_doxo.bed \
--meth gt --gname hg19 --report /root/stuff/diffbind/diffreps/P53_48h_vs_24h.diffreps --frag 0 --nproc 30 

##
#

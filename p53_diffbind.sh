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

trim_galore --illumina --paired -o /root/hct116_wgbs_revised/trim \
/root/wgbs_doxo/OriginalFastq/P007_48h_doxo_1.fq.gz \
/root/wgbs_doxo/OriginalFastq/P007_48h_doxo_2.fq.gz &

trim_galore --illumina --paired -o /root/hct116_wgbs_revised/trim \
/root/wgbs_doxo/OriginalFastq/TP53del_48h_doxo_1.fq.gz \
/root/wgbs_doxo/OriginalFastq/TP53del_48h_doxo_2.fq.gz &

#########################################################################################

~/myPrograms/Bismark/bismark --bowtie2 --multicore 34 \
/root/resources/hg38_bismark_vanilla/ \
-1 /root/HCT116_wgbs/trim/combine_hct116_WGBS_1_val_1.fq.gz \
-2 /root/HCT116_wgbs/trim/combine_hct116_WGBS_2_val_2.fq.gz

~/myPrograms/Bismark/bismark --bowtie2 --multicore 34 \
/root/resources/hg38_bismark_vanilla/ \
-1 /root/hct116_wgbs_revised/trim/P007_48h_doxo_1_val_1.fq.gz \
-2 /root/hct116_wgbs_revised/trim/P007_48h_doxo_2_val_2.fq.gz

~/myPrograms/Bismark/bismark --bowtie2 --multicore 34 \
/root/resources/hg38_bismark_vanilla/ \
-1 /root/hct116_wgbs_revised/trim/TP53del_48h_doxo_1_val_1.fq.gz \
-2 /root/hct116_wgbs_revised/trim/TP53del_48h_doxo_2_val_2.fq.gz

###########
~/myPrograms/Bismark/bismark_methylation_extractor --multicore 32 --gzip --buffer_size 100G --paired-end --ample_memory --comprehensive --cytosine_report --genome_folder /root/resources/hg38_bismark_vanilla/ /root/hct116_wgbs_revised/combine_hct116_WGBS_1_val_1_bismark_bt2_pe.bam
~/myPrograms/Bismark/bismark_methylation_extractor --multicore 32 --gzip --buffer_size 100G --paired-end --ample_memory --comprehensive --cytosine_report --genome_folder /root/resources/hg38_bismark_vanilla/ /root/hct116_wgbs_revised/P007_48h_doxo_1_val_1_bismark_bt2_pe.bam
~/myPrograms/Bismark/bismark_methylation_extractor --multicore 32 --gzip --buffer_size 100G --paired-end --ample_memory --comprehensive --cytosine_report --genome_folder /root/resources/hg38_bismark_vanilla/ /root/hct116_wgbs_revised/TP53del_48h_doxo_1_val_1_bismark_bt2_pe.bam

####

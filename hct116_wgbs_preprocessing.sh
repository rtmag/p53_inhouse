trim_galore --illumina --paired -o /root/HCT116_wgbs/trim \
/root/HCT116_wgbs/combine_hct116_WGBS_1.fq.gz \
/root/HCT116_wgbs/combine_hct116_WGBS_2.fq.gz 

~/myPrograms/Bismark/bismark --bowtie2 --multicore 34 \
/root/resources/HCT116_bs/ \
-1 /root/HCT116_wgbs/trim/combine_hct116_WGBS_1_val_1.fq.gz \
-2 /root/HCT116_wgbs/trim/combine_hct116_WGBS_2_val_2.fq.gz

~/myPrograms/Bismark/bismark_methylation_extractor --multicore 32 --gzip --buffer_size 100G --paired-end --ample_memory --comprehensive --cytosine_report --genome_folder /root/resources/HCT116_bs/ /root/HCT116_wgbs/combine_hct116_WGBS_1_val_1_bismark_bt2_pe.bam

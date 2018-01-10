cat p53_only_10k.bedpe | tail -n +2 | awk 'BEGIN{OFS="\t"} {print "chr"$1,$2,$3,"chr"$4,$5,$6,"loop_"NR,$16,"+","-"}' > p53_only_10k_chr.bedpe
cat ../rao_hct116/merged_loops.bedpe | tail -n +2 | awk 'BEGIN{OFS="\t"} {print "chr"$1,$2,$3,"chr"$4,$5,$6,"loop_"NR,$16,"+","-"}' > rao_merged_loops.bedpe
cat wt_only_10k.bedpe | tail -n +2 | awk 'BEGIN{OFS="\t"} {print "chr"$1,$2,$3,"chr"$4,$5,$6,"loop_"NR,$16,"+","-"}' > wt_only_10k_chr.bedpe

python  ~/myPrograms/liftOverBedpe/liftOverBedpe.py --lift /root/myPrograms/kentUtils/bin/liftOver \
--chain /root/resources/hg19ToHg38.over.chain.gz --i rao_merged_loops.bedpe \
--o rao_merged_loops_hg38.bedpe  --h F

python  ~/myPrograms/liftOverBedpe/liftOverBedpe.py --lift /root/myPrograms/kentUtils/bin/liftOver \
--chain /root/resources/hg19ToHg38.over.chain.gz --i p53_only_10k_chr.bedpe \
--o p53_only_10k_hg38.bedpe  --h F

python  ~/myPrograms/liftOverBedpe/liftOverBedpe.py --lift /root/myPrograms/kentUtils/bin/liftOver \
--chain /root/resources/hg19ToHg38.over.chain.gz --i wt_only_10k_chr.bedpe \
--o wt_only_10k_hg38.bedpe  --h F

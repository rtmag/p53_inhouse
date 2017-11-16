java -jar ~/myPrograms/juicer/scripts/juicer_tools_0.7.5.jar Arrowhead hct_doxo_24_inter.hic hct_doxo_24_inter.arrowhead --ignore_sparsity


# GOLD GPU CLUSTER
java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccups /hpctmp/e0056363/GSE104333_Rao-2017-untreated_combined.hic /hpctmp/e0056363/rao_hct116
java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccups hct_dmso_24_inter.hic /hpctmp/e0056363/hct_dmso_24
java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccups hct_doxo_24_inter.hic /hpctmp/e0056363/hct_doxo_24
java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccups p53_dmso_24_inter.hic /hpctmp/e0056363/p53_dmso_24
java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccups p53_doxo_24_inter.hic /hpctmp/e0056363/p53_doxo_24

java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccups hct_doxo_6_inter.hic /hpctmp/e0056363/hct_doxo_6
java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccups p53_doxo_6_inter.hic /hpctmp/e0056363/p53_doxo_6

##
#DONE



export PATH=/usr/local/cuda-8.0/bin${PATH:+:${PATH}}
# AWS GPU CLUSTER

time ( java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccups --ignore_sparsity GSE104333_Rao-2017-untreated_combined.hic rao_hct116 ) 2> rao_hct116.time > rao_hct116.out
time ( java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccups --ignore_sparsity hct_dmso_24_inter.hic hct_dmso_24 ) 2> hct_dmso_24.time > hct_dmso_24.out
time ( java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccups --ignore_sparsity hct_doxo_24_inter.hic hct_doxo_24 ) 2> hct_doxo_24.time > hct_doxo_24.out
time ( java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccups --ignore_sparsity p53_dmso_24_inter.hic p53_dmso_24 ) 2> p53_dmso_24.time > p53_dmso_24.out
time ( java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccups --ignore_sparsity p53_doxo_24_inter.hic p53_doxo_24 ) 2> p53_doxo_24.time > p53_doxo_24.out

time ( java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccups --ignore_sparsity hct_doxo_6_inter.hic hct_doxo_6 ) 2> hct_doxo_6.time > hct_doxo_6.out
time ( java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccups --ignore_sparsity p53_doxo_6_inter.hic p53_doxo_6 ) 2> p53_doxo_6.time > p53_doxo_6.out

##
#DONE
time ( java -jar juicer_tools.1.7.6_jcuda.0.8.jar hiccupsdiff --ignore_sparsity hct_doxo_24_inter.hic p53_doxo_24_inter.hic /home/ubuntu/hic/hct_doxo_24/merged_loops.bedpe /home/ubuntu/hic/p53_doxo_24/merged_loops.bedpe hiccupsDiff_wt_p53_doxo_24h ) 2> hiccupsDiff_wt_p53_doxo_24h.time > hiccupsDiff_wt_p53_doxo_24h.out


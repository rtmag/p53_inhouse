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
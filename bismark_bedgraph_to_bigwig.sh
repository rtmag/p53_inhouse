
more P007_48h_doxo_1_bismark_bt2_pe.bedGraph | \
awk -F"\t" '{ \
  if($1=='1'){print "chr1\t"$2"\t"$3"\t"$4} \
  if($1=='2'){print "chr2\t"$2"\t"$3"\t"$4} \
  if($1=='3'){print "chr3\t"$2"\t"$3"\t"$4} \
  if($1=='4'){print "chr4\t"$2"\t"$3"\t"$4} \
  if($1=='5'){print "chr5\t"$2"\t"$3"\t"$4} \
  if($1=='6'){print "chr6\t"$2"\t"$3"\t"$4} \
  if($1=='7'){print "chr7\t"$2"\t"$3"\t"$4} \
  if($1=='8'){print "chr8\t"$2"\t"$3"\t"$4} \
  if($1=='9'){print "chr9\t"$2"\t"$3"\t"$4} \
  if($1=='10'){print "chr10\t"$2"\t"$3"\t"$4} \
  if($1=='11'){print "chr11\t"$2"\t"$3"\t"$4} \
  if($1=='12'){print "chr12\t"$2"\t"$3"\t"$4} \
  if($1=='13'){print "chr13\t"$2"\t"$3"\t"$4} \
  if($1=='14'){print "chr14\t"$2"\t"$3"\t"$4} \
  if($1=='15'){print "chr15\t"$2"\t"$3"\t"$4} \
  if($1=='16'){print "chr16\t"$2"\t"$3"\t"$4} \
  if($1=='17'){print "chr17\t"$2"\t"$3"\t"$4} \
  if($1=='18'){print "chr18\t"$2"\t"$3"\t"$4} \
  if($1=='19'){print "chr19\t"$2"\t"$3"\t"$4} \
  if($1=='20'){print "chr20\t"$2"\t"$3"\t"$4} \
  if($1=='21'){print "chr21\t"$2"\t"$3"\t"$4} \
  if($1=='22'){print "chr22\t"$2"\t"$3"\t"$4} \
  if($1=='23'){print "chrX\t"$2"\t"$3"\t"$4} \
  if($1=='24'){print "chrY\t"$2"\t"$3"\t"$4} \
  if($1=='25'){print "chrM\t"$2"\t"$3"\t"$4} \
  if($1=='26'){print "chr1_KI270706v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='27'){print "chr1_KI270707v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='28'){print "chr1_KI270708v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='29'){print "chr1_KI270709v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='30'){print "chr1_KI270710v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='31'){print "chr1_KI270711v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='32'){print "chr1_KI270712v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='33'){print "chr1_KI270713v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='34'){print "chr1_KI270714v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='35'){print "chr2_KI270715v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='36'){print "chr2_KI270716v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='37'){print "chr3_GL000221v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='38'){print "chr4_GL000008v2_random\t"$2"\t"$3"\t"$4} \
  if($1=='39'){print "chr5_GL000208v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='40'){print "chr9_KI270717v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='41'){print "chr9_KI270718v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='42'){print "chr9_KI270719v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='43'){print "chr9_KI270720v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='44'){print "chr11_KI270721v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='45'){print "chr14_GL000009v2\t"$2"\t"$3"\t"$4} \
  if($1=='46'){print "chr14_GL000225v1\t"$2"\t"$3"\t"$4} \
  if($1=='47'){print "chr14_KI270722v1\t"$2"\t"$3"\t"$4} \
  if($1=='48'){print "chr14_GL000194v1\t"$2"\t"$3"\t"$4} \
  if($1=='49'){print "chr14_KI270723v1\t"$2"\t"$3"\t"$4} \
}' > WT_48h_doxo_wgbs.bedGraph


more TP53del_48h_doxo_1_bismark_bt2_pe.bedGraph | \
awk -F"\t" '{ \
  if($1=='1'){print "chr1\t"$2"\t"$3"\t"$4} \
  if($1=='2'){print "chr2\t"$2"\t"$3"\t"$4} \
  if($1=='3'){print "chr3\t"$2"\t"$3"\t"$4} \
  if($1=='4'){print "chr4\t"$2"\t"$3"\t"$4} \
  if($1=='5'){print "chr5\t"$2"\t"$3"\t"$4} \
  if($1=='6'){print "chr6\t"$2"\t"$3"\t"$4} \
  if($1=='7'){print "chr7\t"$2"\t"$3"\t"$4} \
  if($1=='8'){print "chr8\t"$2"\t"$3"\t"$4} \
  if($1=='9'){print "chr9\t"$2"\t"$3"\t"$4} \
  if($1=='10'){print "chr10\t"$2"\t"$3"\t"$4} \
  if($1=='11'){print "chr11\t"$2"\t"$3"\t"$4} \
  if($1=='12'){print "chr12\t"$2"\t"$3"\t"$4} \
  if($1=='13'){print "chr13\t"$2"\t"$3"\t"$4} \
  if($1=='14'){print "chr14\t"$2"\t"$3"\t"$4} \
  if($1=='15'){print "chr15\t"$2"\t"$3"\t"$4} \
  if($1=='16'){print "chr16\t"$2"\t"$3"\t"$4} \
  if($1=='17'){print "chr17\t"$2"\t"$3"\t"$4} \
  if($1=='18'){print "chr18\t"$2"\t"$3"\t"$4} \
  if($1=='19'){print "chr19\t"$2"\t"$3"\t"$4} \
  if($1=='20'){print "chr20\t"$2"\t"$3"\t"$4} \
  if($1=='21'){print "chr21\t"$2"\t"$3"\t"$4} \
  if($1=='22'){print "chr22\t"$2"\t"$3"\t"$4} \
  if($1=='23'){print "chrX\t"$2"\t"$3"\t"$4} \
  if($1=='24'){print "chrY\t"$2"\t"$3"\t"$4} \
  if($1=='25'){print "chrM\t"$2"\t"$3"\t"$4} \
  if($1=='26'){print "chr1_KI270706v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='27'){print "chr1_KI270707v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='28'){print "chr1_KI270708v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='29'){print "chr1_KI270709v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='30'){print "chr1_KI270710v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='31'){print "chr1_KI270711v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='32'){print "chr1_KI270712v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='33'){print "chr1_KI270713v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='34'){print "chr1_KI270714v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='35'){print "chr2_KI270715v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='36'){print "chr2_KI270716v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='37'){print "chr3_GL000221v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='38'){print "chr4_GL000008v2_random\t"$2"\t"$3"\t"$4} \
  if($1=='39'){print "chr5_GL000208v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='40'){print "chr9_KI270717v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='41'){print "chr9_KI270718v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='42'){print "chr9_KI270719v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='43'){print "chr9_KI270720v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='44'){print "chr11_KI270721v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='45'){print "chr14_GL000009v2\t"$2"\t"$3"\t"$4} \
  if($1=='46'){print "chr14_GL000225v1\t"$2"\t"$3"\t"$4} \
  if($1=='47'){print "chr14_KI270722v1\t"$2"\t"$3"\t"$4} \
  if($1=='48'){print "chr14_GL000194v1\t"$2"\t"$3"\t"$4} \
  if($1=='49'){print "chr14_KI270723v1\t"$2"\t"$3"\t"$4} \
}' > P53_48h_doxo_wgbs.bedGraph



more combine_hct116_WGBS_1_val_1_bismark_bt2_pe.bedGraph | \
awk -F"\t" '{ \
  if($1=='1'){print "chr1\t"$2"\t"$3"\t"$4} \
  if($1=='2'){print "chr2\t"$2"\t"$3"\t"$4} \
  if($1=='3'){print "chr3\t"$2"\t"$3"\t"$4} \
  if($1=='4'){print "chr4\t"$2"\t"$3"\t"$4} \
  if($1=='5'){print "chr5\t"$2"\t"$3"\t"$4} \
  if($1=='6'){print "chr6\t"$2"\t"$3"\t"$4} \
  if($1=='7'){print "chr7\t"$2"\t"$3"\t"$4} \
  if($1=='8'){print "chr8\t"$2"\t"$3"\t"$4} \
  if($1=='9'){print "chr9\t"$2"\t"$3"\t"$4} \
  if($1=='10'){print "chr10\t"$2"\t"$3"\t"$4} \
  if($1=='11'){print "chr11\t"$2"\t"$3"\t"$4} \
  if($1=='12'){print "chr12\t"$2"\t"$3"\t"$4} \
  if($1=='13'){print "chr13\t"$2"\t"$3"\t"$4} \
  if($1=='14'){print "chr14\t"$2"\t"$3"\t"$4} \
  if($1=='15'){print "chr15\t"$2"\t"$3"\t"$4} \
  if($1=='16'){print "chr16\t"$2"\t"$3"\t"$4} \
  if($1=='17'){print "chr17\t"$2"\t"$3"\t"$4} \
  if($1=='18'){print "chr18\t"$2"\t"$3"\t"$4} \
  if($1=='19'){print "chr19\t"$2"\t"$3"\t"$4} \
  if($1=='20'){print "chr20\t"$2"\t"$3"\t"$4} \
  if($1=='21'){print "chr21\t"$2"\t"$3"\t"$4} \
  if($1=='22'){print "chr22\t"$2"\t"$3"\t"$4} \
  if($1=='23'){print "chrX\t"$2"\t"$3"\t"$4} \
  if($1=='24'){print "chrY\t"$2"\t"$3"\t"$4} \
  if($1=='25'){print "chrM\t"$2"\t"$3"\t"$4} \
  if($1=='26'){print "chr1_KI270706v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='27'){print "chr1_KI270707v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='28'){print "chr1_KI270708v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='29'){print "chr1_KI270709v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='30'){print "chr1_KI270710v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='31'){print "chr1_KI270711v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='32'){print "chr1_KI270712v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='33'){print "chr1_KI270713v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='34'){print "chr1_KI270714v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='35'){print "chr2_KI270715v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='36'){print "chr2_KI270716v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='37'){print "chr3_GL000221v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='38'){print "chr4_GL000008v2_random\t"$2"\t"$3"\t"$4} \
  if($1=='39'){print "chr5_GL000208v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='40'){print "chr9_KI270717v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='41'){print "chr9_KI270718v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='42'){print "chr9_KI270719v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='43'){print "chr9_KI270720v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='44'){print "chr11_KI270721v1_random\t"$2"\t"$3"\t"$4} \
  if($1=='45'){print "chr14_GL000009v2\t"$2"\t"$3"\t"$4} \
  if($1=='46'){print "chr14_GL000225v1\t"$2"\t"$3"\t"$4} \
  if($1=='47'){print "chr14_KI270722v1\t"$2"\t"$3"\t"$4} \
  if($1=='48'){print "chr14_GL000194v1\t"$2"\t"$3"\t"$4} \
  if($1=='49'){print "chr14_KI270723v1\t"$2"\t"$3"\t"$4} \
}' > WT_0h_doxo_wgbs.bedGraph

~/myPrograms/kentUtils/bin/bedGraphToBigWig WT_0h_doxo_wgbs.bedGraph ~/resources/hg38.chrom.sizes WT_0h_doxo_wgbs.bw &

~/myPrograms/kentUtils/bin/bedGraphToBigWig WT_48h_doxo_wgbs.bedGraph ~/resources/hg38.chrom.sizes WT_48h_doxo_wgbs.bw &

~/myPrograms/kentUtils/bin/bedGraphToBigWig P53_48h_doxo_wgbs.bedGraph ~/resources/hg38.chrom.sizes P53KO_48h_doxo_wgbs.bw &


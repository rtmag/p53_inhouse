
tss = read.table("~/resources/hg38_tss.bed",sep="\t",header=F,stringsAsFactors=F)

c25 = read.table("cluster2_5_lsgenes_names.txt",sep="\t",header=F,stringsAsFactors=F)

c134 = read.table("cluster1_3_4_names.txt",sep="\t",header=F,stringsAsFactors=F)
#

table(tss[,4] %in% c25$V1)
table(tss[,4] %in% c134$V1)
#table(c25$V1 %in% tss[,4])
#c25[,1][!c25$V1 %in% tss[,4]]

write.table(tss[tss[,4] %in% c25$V1,], "c25_tss.bed", quote=F,col.names=F,row.names=F,sep="\t")
write.table(tss[tss[,4] %in% c134$V1,], "c134_tss.bed",quote=F,col.names=F,row.names=F,sep="\t")


####



gene = read.table("~/resources/hg38_genes.bed",sep="\t",header=F,stringsAsFactors=F)

c25 = read.table("cluster2_5_lsgenes_names.txt",sep="\t",header=F,stringsAsFactors=F)

c134 = read.table("cluster1_3_4_names.txt",sep="\t",header=F,stringsAsFactors=F)
#

table(gene[,4] %in% c25$V1)
table(gene[,4] %in% c134$V1)
#table(c25$V1 %in% tss[,4])
#c25[,1][!c25$V1 %in% tss[,4]]

write.table(gene[gene[,4] %in% c25$V1,], "c25_gene.bed", quote=F,col.names=F,row.names=F,sep="\t")
write.table(gene[gene[,4] %in% c134$V1,], "c134_gene.bed",quote=F,col.names=F,row.names=F,sep="\t")

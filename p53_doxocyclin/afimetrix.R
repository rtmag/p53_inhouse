# followed https://github.com/icnn/Microarray-Tutorials/wiki/Affymetrix#4
## Load packages
library(affy)   # Affymetrix pre-processing
library(limma)  # two-color pre-processing; differential
                  # expression
         
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

#gset <- getGEO("GSE30753", GSEMatrix =TRUE, getGPL=F)
#if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
#gset <- gset[[idx]]
#pData(gset)

data.affy <- ReadAffy(celfile.path = "/home/rtm/Downloads/p53_doxocyclin/GSM", filenames = list.files())
datExpr <- exprs(data.affy)
datExpr <- log2(datExpr)

###################################################################################################################
boxplot(datExpr,range=0, xaxt='n', xlab = "Array", main = "Boxplot", ylab = "Intensity")

i=1 
plot(density((datExpr[,i]),na.rm=T),
     main = "Histogram", xlab="log2 exp")
for(i in 2:dim(datExpr)[2]){
  lines(density((datExpr[,i]),na.rm=T))
}
###################################################################################################################
datExpr_norm <- rma(data.affy, background=T, normalize=T, verbose=T)
datExpr <- exprs(datExpr_norm)
datExpr <- log2(datExpr)

library( biomaRt )
library(annotate)
library(hugene10sttranscriptcluster.db)
annodb <- "hugene10sttranscriptcluster.db"
ID     <- featureNames(data.affy)
Symbol <- as.character(lookUp(ID, annodb, "SYMBOL"))
Name   <- as.character(lookUp(ID, annodb, "GENENAME"))
Entrez <- as.character(lookUp(ID, annodb, "ENTREZID"))

phenoTable = pData(data.affy)
rownames(datExpr) <- Symbol

datExpr = datExpr[Symbol!="NA",]
#
sample_names = colnames(datExpr)
sample_names = gsub("0ng","NONE",sample_names)
sample_names = gsub("3ng","LOW",sample_names)
sample_names = gsub("5ng","HIGH",sample_names)
sample_names = gsub(".+NONE","NONE",sample_names,perl=T)
sample_names = gsub(".+LOW","LOW",sample_names,perl=T)
sample_names = gsub(".+HIGH","HIGH",sample_names,perl=T)
sample_names = gsub("h\\-.+","h",sample_names,perl=T)
sample_names = paste(sample_names,1:4,sep="_")
########################################################
colnames(datExpr) = sample_names
saveRDS(datExpr,"p53_doxocyclin.RDS")
###################################################################################################################
###################################################################################################################
###################################################################################################################
options(scipen=999)
library(gplots)
library(pheatmap)
library(factoextra)
library(RColorBrewer)
colors <- c("green","white","red")
 colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(20))

data = readRDS("p53_doxocyclin.RDS")
data = data[,c(9:12,1:4,13:16,5:8,17:20)]

c25 = read.table("cluster2_5_lsgenes_names.txt")
c134 = read.table("cluster1_3_4_names.txt")

c25 = as.character(c25[c25[,1] %in% rownames(data),1])
c134 = as.character(c134[c134[,1] %in% rownames(data),1])  
        
c12345 = c(c25,c134)
matrix = data[rownames(data) %in% c12345,]
rlab = data.frame(genes = c(rep("down",length(c25)),rep("up",length(c134))))


clab_col = list(genes = c(down="darkred",up="darkgreen"))


clab_col = list( OTHER_tr=c(white="white",black="black"),NUTILIN=c(white="white",black="black"),
                DOXO=c(white="white",black="black"),p53OE=c(white="white",black="black"),
           P21=c(white="white",blue="blue"),
           OTHER_cell=c(white="white",red="red"),MCF7=c(white="white",red="red"),
                U2OS=c(white="white",red="red"),HCT116=c(white="white",red="red"),
genes = c(down="darkred",up="darkgreen"))


png("p53_doxocyclin_SelectedDream.png",width= 5.25,
  height=10,units="in",
res=1200,pointsize=4)

pheatmap(matrix,col=colors,show_rownames=F,annotation_row=data.frame(rlab),cellwidth=10,
        annotation_legend = F,scale="row",cluster_cols=F,annotation_colors = clab_col)
dev.off()

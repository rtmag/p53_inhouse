# followed https://github.com/icnn/Microarray-Tutorials/wiki/Affymetrix#4
## Load packages
library(affy)   # Affymetrix pre-processing
library(limma)  # two-color pre-processing; differential
                  # expression
         
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE30753", GSEMatrix =TRUE, getGPL=F)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
pData(gset)

data.affy <- ReadAffy(celfile.path = "/home/rtm/Downloads/p53_doxocyclin/GSM", filenames = list.files())
datExpr <- exprs(data.affy)

datExpr <- exprs(gset)
datExpr <- log2(datExpr)

boxplot(datExpr,range=0, xaxt='n', xlab = "Array", main = "Boxplot", ylab = "Intensity")
#

i=1 
plot(density((datExpr[,i]),na.rm=T),
     main = "Histogram", xlab="log2 exp")
for(i in 2:dim(datExpr)[2]){
  lines(density((datExpr[,i]),na.rm=T))
}

datExpr_norm <- rma(data.affy, background=T, normalize=T, verbose=T)
datExpr <- exprs(datExpr)
datExpr <- log2(datExpr)

library( biomaRt )

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")

f <- listFilters(ensembl)
a <- listAttributes(ensembl)

identifier <- "affy_hg_u133a"
getinfo <- c("affy_hg_u133a", "ensembl_gene_id", "entrezgene", "external_gene_name")
geneDat <- getBM(attributes = getinfo, filters=identifier, values = rownames(exprs(data.affy)),mart=ensembl)


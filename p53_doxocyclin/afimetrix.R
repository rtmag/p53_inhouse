## Load packages
library(affy)   # Affymetrix pre-processing
library(limma)  # two-color pre-processing; differential
                  # expression
                  
## import "phenotype" data, describing the experimental design
phenoData <- 
    read.AnnotatedDataFrame(system.file("extdata", "pdata.txt",
    package="arrays"))

## RMA normalization
celfiles <- system.file("extdata", package="arrays")
eset <- justRMA(phenoData=phenoData,
    celfile.path=celfiles)
    
    
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE30753", GSEMatrix =TRUE, getGPL=F)
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
pdata(gset)
# set parameters and draw the plot

dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE30753", '/', annotation(gset), " selected samples", sep ='')
boxplot(exprs(gset), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)


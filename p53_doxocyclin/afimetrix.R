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
    
    
    

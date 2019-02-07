x = read.csv("p53_meta_analysis_edited.csv")
data = x[4:dim(x)[1],2:dim(x)[2]]
rownames(data) = as.character(x[4:dim(x)[1],1])


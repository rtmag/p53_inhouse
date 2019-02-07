source('https://raw.githubusercontent.com/rtmag/tumor-meth-pipe/master/heatmap3.R')
options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)
colors <- rev(colorRampPalette( (brewer.pal(3, "RdBu")) )(3))
colors <- c("green","black","red")

x = read.csv("p53_meta_analysis_edited.csv")
data = x[4:dim(x)[1],2:dim(x)[2]]
data = matrix(as.numeric(as.matrix(data)),ncol=22)
rownames(data) = as.character(x[4:dim(x)[1],1])
colnames(data) = colnames(x[2:dim(x)[2]])
data=data+2

treatment = as.character(unlist(x[1,2:dim(x)[2]]))
p21 = as.character(unlist(x[2,2:dim(x)[2]]))
cell = as.character(unlist(x[3,2:dim(x)[2]]))

OTHER_tr = treatment
OTHER_tr[!OTHER_tr %in% c("NUTILIN","DOXO","p53OE")]=2
OTHER_tr[OTHER_tr %in% c("NUTILIN","DOXO","p53OE")]=1
OTHER_tr=as.numeric(OTHER_tr)

NUTILIN = treatment
NUTILIN[NUTILIN!='NUTILIN']=1
NUTILIN[NUTILIN=='NUTILIN']=2
NUTILIN=as.numeric(NUTILIN)

DOXO = treatment
DOXO[DOXO!='DOXO']=1
DOXO[DOXO=='DOXO']=2
DOXO=as.numeric(DOXO)

p53OE = treatment
p53OE[p53OE!='p53OE']=1
p53OE[p53OE=='p53OE']=2
p53OE=as.numeric(p53OE)

p21[p21=='WT']=1
p21[p21=='KO']=2
p21=as.numeric(p21)

OTHER_cell = cell
OTHER_cell[!OTHER_cell %in% c("MCF7","U2OS","HCT116")]=2
OTHER_cell[OTHER_cell %in% c("MCF7","U2OS","HCT116")]=1
OTHER_cell=as.numeric(OTHER_cell)

MCF7 = cell
MCF7[MCF7!='MCF7']=1
MCF7[MCF7=='MCF7']=2
MCF7=as.numeric(MCF7)

U2OS = cell
U2OS[U2OS!='U2OS']=1
U2OS[U2OS=='U2OS']=2
U2OS=as.numeric(U2OS)

HCT116 = cell
HCT116[HCT116!='HCT116']=1
HCT116[HCT116=='HCT116']=2
HCT116=as.numeric(HCT116)

col1 = c('white','red')
col2 = c('white','blue')
col3 = c('white','black')
clab=cbind(OTHER_tr=col3[OTHER_tr],NUTILIN=col3[NUTILIN],DOXO=col3[DOXO],p53OE=col3[p53OE],
           P21=col2[p21],
           OTHER_cell=col1[OTHER_cell],MCF7=col1[MCF7],U2OS=col1[U2OS],HCT116=col1[HCT116])

rownames(clab) = colnames(data)

clab_col = list( OTHER_tr=col3,NUTILIN=col3,DOXO=col3,p53OE=col3,
           P21=col2,
           OTHER_cell=col1,MCF7=col1,U2OS=col1,HCT116=col1 )


library(pheatmap)
pheatmap(matrix,col=colors,show_rownames=F,annotation_col=data.frame(clab),annotation_colors = clab_col)
##########################################################################################

dream = read.csv("dream_targets.csv")
dream = as.character(dream[,1])

matrix = data[rownames(data) %in% dream,]
matrix = matrix[abs(rowSums(matrix))>0,]

heatmap.2(matrix,col=colors,scale="none", trace="none",srtCol=90,
labRow = FALSE,xlab="",dendogram="none")

heatmap(matrix,labRow = FALSE,col=colors,scale="row",distfun = function(x) get_dist(x,method="pearson"))
        
c25 = read.table("/home/rtm/Downloads/cluster2_5_lsgenes_names.txt")
c25 = as.character(c25[,1])
           
c134 = read.table("/home/rtm/Downloads/cluster1_3_4_names.txt")
c134 = as.character(c134[,1])  
        
c12345 = c(c25,c134)
        
matrix = data[rownames(data) %in% c12345,]
matrix = matrix[abs(rowSums(matrix))>0,]

clab=cbind(colores[tp53],colores[braf],colores[kras],col.stage[stage])
        
heatmap.3(matrix,labRow = FALSE,col=colors,scale="none",ColSideColors=clab,margins=c(2,9))
         )

ColSideColors=clab

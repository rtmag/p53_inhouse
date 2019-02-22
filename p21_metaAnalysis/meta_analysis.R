source('https://raw.githubusercontent.com/rtmag/tumor-meth-pipe/master/heatmap3.R')
options(scipen=999)
library(gplots)
library(pheatmap)
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
#####################################
#####################################
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
#####################################
#####################################
col1 = c('white','red')
col2 = c('white','blue')
col3 = c('white','black')
clab=cbind(OTHER_tr=col3[OTHER_tr],NUTILIN=col3[NUTILIN],DOXO=col3[DOXO],p53OE=col3[p53OE],
           P21=col2[p21],
           OTHER_cell=col1[OTHER_cell],MCF7=col1[MCF7],U2OS=col1[U2OS],HCT116=col1[HCT116])

clab_col = list( OTHER_tr=c(white="white",black="black"),NUTILIN=c(white="white",black="black"),
                DOXO=c(white="white",black="black"),p53OE=c(white="white",black="black"),
           P21=c(white="white",blue="blue"),
           OTHER_cell=c(white="white",red="red"),MCF7=c(white="white",red="red"),
                U2OS=c(white="white",red="red"),HCT116=c(white="white",red="red"),
               genes = c(down="darkred",up="darkgreen"))

rownames(clab) = colnames(data)


library(pheatmap)

c25 = read.table("cluster2_5_lsgenes_names.txt")
c134 = read.table("cluster1_3_4_names.txt")

c25 = as.character(c25[c25[,1] %in% rownames(data),1])
c134 = as.character(c134[c134[,1] %in% rownames(data),1])  
        
c12345 = c(c25,c134)
        
matrix = data[rownames(data) %in% c12345,]

rlab = data.frame(genes = c(rep("down",length(c25)),rep("up",length(c134))))
rownames(rlab) = rownames(matrix)

png("p53_metaAnalysis_stephProjected.png",width= 5.25,
  height=10,units="in",
  res=1600,pointsize=4)
pheatmap(matrix,col=colors,show_rownames=F,annotation_col=data.frame(clab),annotation_row=data.frame(rlab),
        annotation_legend = FALSE,annotation_colors = clab_col,cellheight=.50,cellwidth=10,cluster_cols=F)
library(grid)
grid.ls(grid.force()) # "col_annotation" looks like it's the one to edit
grid.gedit("col_annotation", gp = gpar(col="black"))
dev.off()
##########################################################################################

selected_dream = c("EIF2C2",
"ATAD2",
"ATAD5",
"AURKA",
"AURKB",
"BIRC5",
"BLM",
"BUB1",
"CCNA2",
"CCNB2",
"CDC20",
"CDC25A",
"CDC25C",
"CDC2",
"CDKN1A",
"CDKN2A",
"CENPA",
"CENPE",
"CHEK1",
"DSCC1",
"E2F1",
"E2F7",
"E2F8",
"EXO1",
"FANCA",
"FEN1",
"FOXM1",
"GINS1",
"KIF14",
"KIF20A",
"KIF2C",
"KIF4A",
"KIF4B",
"MCM2",
"MCM4",
"MELK",
"MTBP",
"FAM54A",
"MYBL2",
"NUF2",
"PLK1",
"PRC1",
"PTTG1",
"PTTG2",
"RAD51",
"RAD51AP1",
"RAD54L",
"RBL1",
"RFC4",
"SKA1",
"SLC7A5",
"TPX2",
"TRIP13",
"TROAP",
"TTK",
"UBE2C")

matrix = data[rownames(data) %in% selected_dream,]

png("p53_metaAnalysis_SelectedDream.png",width= 5.25,
  height=10,units="in",
  res=1200,pointsize=4)
pheatmap(matrix,col=colors,show_rownames=T,annotation_col=data.frame(clab),
        annotation_legend = FALSE,annotation_colors = clab_col,cellheight=9,cellwidth=8)
library(grid)
grid.ls(grid.force()) # "col_annotation" looks like it's the one to edit
grid.gedit("col_annotation", gp = gpar(col="black"))
dev.off()

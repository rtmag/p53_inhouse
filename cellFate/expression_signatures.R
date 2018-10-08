# SENESENCE 
# Clasical markers
sns_class = c("CDKN1A","CDKN2A","CDKN2B","CDK4","CDK6","CDK2","RB1","NOTCH3","LMNB1",
              "CXCL1","CXCL2","IL6","IL8","IL10","SIRT1","MAPK11","GLB1","NOTCH1")

# 10 gene Marker Althubiti, et al. Cell Death and Disease 2014.
SNS_10 = c("DEP10","NTAL","EBP50","STX4","VAMP3","ARMCX3","B2MG","LANCL1","VPS26A","PLD3")

# 20 gene Marker Nagano, et al. Scientific Reports 2016.
SNS_20 = c("PVRL4","GPR172B","DAO","CCDC74B","LOXL4","EVL","PRODH","E2F7","LY6D","IGFBP2","CRABP2","EPN3","APOBEC3B","IER5",
          "ANGPTL2","SLC48AI","WBSCR27","E2F2","NXPH4","PPMID","CDKN1A","BTG2","SULF2")

SNS_20_p53_dependent = c("PVRL4","PRODH","LY6D","DAO","EPN3","GPR172B")

# Apoptosis RT2 Profiler PCR Array https://www.qiagen.com/ca/shop/pcr/primer-sets/rt2-profiler-pcr-arrays/?catno=PAHS-012Z#geneglobe
DNA_damage_repair = c("ABL1", "CIDEA", "CIDEB", "TP53", "TP73")
extracellular_signal = c("CFLAR", "DAPK1", "TNFRSF25")
pro_apoptotic = c("BAD", "BAK1", "BAX", "BCL10", "BCL2L11", "BID", "BIK", "BNIP3", "BNIP3L", "CASP1", "CASP10", "CASP14", "CASP2",
 "CASP3", "CASP4", "CASP6", "CASP8", "CD27", "CD70", "CYCS", "DFFA", "DIABLO", "FAS", "FASLG", "GADD45A", "HRK", "LTA", "NOD1", 
 "PYCARD", "TNFRSF10A", "TNFRSF9", "TNFSF10", "TNFSF8", "TP53BP2", "TRADD", "TRAF3")

anti-apoptotic = c("AKT1", "BAG1", "BAG3", "BAX", "BCL2", "BCL2A1", "BCL2L1", "BCL2L10", "BCL2L2", "BFAR", "BIRC3", "BIRC5", 
"BIRC6", "BNIP2", "BNIP3", "BNIP3L", "BRAF", "CD27", "CD40LG", "CFLAR", "DAPK1", "FAS", "HRK", "IGF1R", "IL10", "MCL1", "NAIP", 
"NFKB1", "NOL3", "RIPK2", "TNF", "XIAP")

negative_regulation_apoptosis = c("BAG1", "BAG3", "BCL10", "BCL2", "BCL2A1", "BCL2L1", "BCL2L10", "BCL2L2",
 "BFAR", "BIRC2", "BIRC3", "BIRC6", "BNIP2", "BNIP3", "BNIP3L", "BRAF", "CASP3", "CD27", "CD40LG", "CFLAR", 
 "CIDEA", "DAPK1", "DFFA", "FAS", "IGF1R", "MCL1", "NAIP", "NOL3", "TP53", "TP73", "XIAP")
positive_regulation_apoptosis = c("ABL1", "AKT1", "BAD", "BAK1", "BAX", "BCL2L11", "BID", "BIK", "BNIP3", "BNIP3L",
 "CASP1", "CASP10", "CASP14", "CASP2", "CASP4", "CASP6", "CASP8", "CD40", "CD70", "CIDEB", "CRADD", "FADD", "FASLG",
  "HRK", "LTA", "LTBR", "NOD1", "PYCARD", "RIPK2", "TNF", "TNFRSF10A", "TNFRSF10B", "TNFRSF25", "TNFRSF9", "TNFSF10",
  "TNFSF8", "TP53", "TP53BP2", "TRADD", "TRAF2", "TRAF3")

death_domain_receptor = c("CRADD", "DAPK1", "FADD", "TNFRSF10A", "TNFRSF10B", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", 
"TNFRSF21", "TNFRSF25", "TRADD","TNF")

caspases = c( "CASP1", "CASP10", "CASP14", "CASP2", "CASP3", "CASP4", "CASP5", "CASP6", "CASP7", "CASP8", "CASP9",
 "CFLAR", "CRADD", "PYCARD")
caspase_activation = c("AIFM1", "APAF1", "BAX", "BCL2L10", "CASP1", "CASP9", "NOD1", "PYCARD", "TNFRSF10A", "TNFRSF10B", "TP53")
caspase_inhibition = c("CD27", "XIAP")

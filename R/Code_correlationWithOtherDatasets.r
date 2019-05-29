# This code makes correlation plots comparing MTG with other data sets.

# FIRST NEED TO SET THE WORKING DIRECTORY TO THE MAIN FOLDER
# This is the folder where "Start.RData" is located

# SECOND, create a "comparison" directory, download the following files, and save them there:
# Lake et al 2016: NOTE, YOU NEED TO ASK PERMISSION TO ACCESS THESE DATA: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/molecular.cgi?study_id=phs000833.v7.p1&phv=219211&phd=4779&pha=&pht=4360&phvf=&phdf=&phaf=&phtf=1&dssp=1&consent=&temp=1
# Lake et al 2018:  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97930&format=file&file=GSE97930%5FFrontalCortex%5FsnDrop%2Dseq%5FUMI%5FCount%5FMatrix%5F08%2D01%2D2017%2Etxt%2Egz"
# Habib et al: https://storage.googleapis.com/gtex_additional_datasets/single_cell_data/GTEx_droncseq_hip_pcf.tar
# Boldog et al: https://github.com/AllenInstitute/L1_rosehip/tree/master/data/human 

#-----------------------------------------------------------------------------------------------
print("Read in the data and load relevant libraries")
outputFolder  = "output/"
scriptsFolder = "R/"
dataFolder    = "data/"
mtgFolder     = "MTG/"
compFolder    = "comparison/"

#-----------------------------------------------------------------------------------------------
# Load these libraries
library(stringr)
library(lowcat)
library(feather)
library(dplyr)
library(scrattch.io)
library(matrixStats)
library(ggplot2)
options(stringsAsFactors = F)

# Read in the extra scripts
source(paste0(scriptsFolder,"Support_extraFunctions.r"))

#-----------------------------------------------------------------------------------------------
# Read in the MTG data
anno       <- read_feather(paste(mtgFolder,"anno.feather",sep="")) 
Expr.dat   <- feather(paste(mtgFolder,"data.feather",sep=""))   # FPKM
Expr.dat   <- Expr.dat[match(anno$sample_id,Expr.dat$sample_id),] # Make sure the expression matches the sample information
datIn      <- as.matrix(Expr.dat[,names(Expr.dat)!="sample_id"])
rownames(datIn) <- Expr.dat$sample_id
datIn      <- t(datIn)
norm.dat   <- log2(datIn+1)


# Get median expression per cluster and the proportions
clm        = anno$cluster_id
names(clm)  = Expr.dat$sample_id
exprThresh = 1
medianExpr = do.call("cbind", tapply(names(clm), clm, function(x) rowMedians(datIn[,x]))) 
propExpr   = do.call("cbind", tapply(names(clm), clm, function(x) rowMeans(norm.dat[,x]>exprThresh))) 
rownames(medianExpr) <- rownames(propExpr) <- rownames(norm.dat)  
colnames(medianExpr) <- colnames(propExpr) <- anno$cluster_label[match(colnames(medianExpr),as.character(anno$cluster_id))]


# Select features based on beta score >0.4
betaScore  = getBetaScore(propExpr)
marker_scores_genes_subset <- names(betaScore)[betaScore>0.4]



#-----------------------------------------------------------------------------------------------
# Define the study data

clusters <- list(Lake2018="Lake2018_clusters.csv",
                 Habib2017="Habib2017_clusters.csv",
				 Boldog="Boldog2018_clusters.csv")
datas    <- list(Lake2018="GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt",
                 Habib2017="GTEx_droncseq_hip_pcf.umi_counts.txt",
				 Boldog2018="data_boldog.txt")
corVal   <- list(Lake2018=0.4, Habib2017=0.3, Boldog2018=0)
				 
# Habib is actually PFC and hippocampus and has N=11317 if we exclude cells mapping to hippocampus clusters
# Boldog has 872 non-outlier clusters
# There are different correlation cutoffs for the different data sets as currently presented
				 
for (s in 1:length(datas)){
  # Read in the study data
  dataIn     <- read.csv(paste0(compFolder,datas[[s]]), sep = "\t", row.names=1)
  dataIn     <- as.matrix(dataIn)
  annotation <- read.csv(paste0(dataFolder,clusters[[s]]),row.names=1)
  rownames(annotation) <- make.names(rownames(annotation))
  data       <- dataIn[,rownames(annotation)] # make the order match
  
  # subsetting comparison (comp) data
  new_comp_frame <- log2(data[rownames(data) %in% marker_scores_genes_subset,]+1)
  new_ref_count  <- log2(medianExpr[rownames(new_comp_frame),]+1) 

  # Run Pearson correlation (from lowcat)
  comp_annotated <- max_column_correlation(new_comp_frame, new_ref_count, method = "pearson")
  names(comp_annotated) <- c("sample_id","cluster_cor","cluster_label")

  # Plot the results
  kpSamp  <- comp_annotated$cluster_cor >= corVal[[s]]
  comp.cl <- droplevels(factor(comp_annotated$cluster_label[kpSamp],levels <- colnames(medianExpr))) 
  ref.cl  <- droplevels(factor(annotation$cluster[kpSamp],levels=sort(unique(annotation$cluster),decreasing =TRUE)))
  plotOut <- plot_annotation_comparison(comp.cl,ref.cl,anno,names(datas)[s]) 
  plotOut

  # Save the plots
  ggsave(paste0(outputFolder,"FigED5_",names(datas)[s],"_correlation.pdf"),plotOut)
}








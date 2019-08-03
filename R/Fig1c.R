print("Read in the data and load relevant libraries")
outputFolder  = "output/"
scriptsFolder = "R/"
inputFolder   = "data/"
dataFolder    = "MTG/"

library(pheatmap)
library(RColorBrewer)
library(feather)


anno     <- read_feather(paste(dataFolder,"anno.feather",sep="")) 
anno.df <- unique(anno[, c("cluster_id", "cluster_label")])
anno$cluster_label <- factor(anno$cluster_label, 
                             levels = anno.df$cluster_label[order(anno.df$cluster_id)])

layer.cl <- as.matrix(table(anno$brain_subregion_label, anno$cluster_label))
layer.cl.prop <- sweep(layer.cl, 2, colSums(layer.cl), "/")


#### Figure 1c ####
hm.colors <- colorRampPalette(c("white", brewer.pal(9, "YlOrRd")))(100)
pheatmap(layer.cl.prop, cluster_rows = FALSE, cluster_cols = FALSE, color = hm.colors)

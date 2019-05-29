# This code identifies and plots genes with common patterns in mouse and human

# FIRST NEED TO SET THE WORKING DIRECTORY TO THE MAIN FOLDER
# This is the folder where "Start.RData" is located

outputFolder  = "output/"
scriptsFolder = "R/"
inputFolder   = "data/"
humanFolder   = "MTG/"
mouseFolder   = "mouse/"

library(Hmisc)
library(gplots)
library(dplyr);
library(feather)

source(paste0(scriptsFolder,"Support_heatmap_functions.R")) 
source(paste0(scriptsFolder,"Support_violin_functions.R"))

########################################################################################
print("Read in the annotations and homology conversions.")

annoH <- read_feather(paste(humanFolder,"anno.feather",sep="")) 
annoM <- read_feather(paste(mouseFolder,"anno.feather",sep="")) 
mhConvert <- read.csv(paste(inputFolder,"human_mouse_homology.csv",sep="")) 


########################################################################################
print("Update anno.feather files with homology information above.")

mhOut <- mhConvert[,c("homol_type","tree_order","cluster_color")]
colnames(mhOut) <- c("homology_cluster_label","homology_cluster_id","homology_cluster_color")
mhHuman <- as.data.frame(mhOut[match(annoH$cluster_label,mhConvert$cluster),])
mhMouse <- as.data.frame(mhOut[match(annoM$cluster_label,mhConvert$cluster),])

## Update human

if(!file.exists(paste0(humanFolder,"anno_original.feather"))) 
  write_feather(annoH,paste0(humanFolder,"anno_original.feather"))

annoH2 <- cbind(annoH,as.data.frame(mhHuman))  

file.remove(paste0(humanFolder,"anno.feather"))
write_feather(annoH2,paste0(humanFolder,"anno.feather"))

## Update mouse

if(!file.exists(paste0(mouseFolder,"anno_original.feather"))) 
  write_feather(annoM,paste0(mouseFolder,"anno_original.feather"))

annoM2 <- cbind(annoM,as.data.frame(mhMouse))  

file.remove(paste0(mouseFolder,"anno.feather"))
write_feather(annoM2,paste0(mouseFolder,"anno.feather"))



########################################################################################
print("Plot serotonin receptors for Fig 6e.")

serotoninH = c("HTR1D","HTR1F","HTR2A","HTR2C","HTR4","HTR5A","HTR6","HTR7","HTR3A","HTR3B")
  
human_plot <- group_heatmap_plot(data_source = humanFolder,
                   genes = serotoninH,
                   group_by = "homology_cluster",
                   clusters = 1:37,
                   calculation = "trimmed_mean",
                   labelheight = 50,
                   showcounts = TRUE)
ggsave(paste0(outputFolder,"Fig6e_serotonin_genes_human.pdf"),human_plot,height = 4, width = 6)
 
mouse_plot <- group_heatmap_plot(data_source = mouseFolder,
                   genes = capitalize(tolower(serotoninH)),
                   group_by = "homology_cluster",
                   clusters = 1:37,
                   calculation = "trimmed_mean",
                   labelheight = 50,
                   showcounts = TRUE)
ggsave(paste0(outputFolder,"Fig6e_serotonin_genes_mouse.pdf"),mouse_plot,height = 4, width = 6)




########################################################################################
print("Plot common marker genes in mouse and human (identified elsewhere) for Fig ED13.")

common_genes = c("SV2C", "KIRREL", "NOS1", "HCRTR2", "TGFB2", "SP8", "PAX6", "FOXO1", "ST3GAL4",
                 "MYLK", "NDNF", "NR2E1", "RXRG", "NDST4", "SALL1", "ANGPT1", "SHISA8", "PPAPDC1A", 
				 "PTGDS", "STAC", "CDHR1", "NDRG1", "LDLR", "RNF152", "NTN1", "CHODL", "CUEDC1", "EOGT", 
				 "TACR1", "CRHBP", "TMCC3", "ARHGEF40", "MOXD1", "B3GAT2", "PYGO1", "KLHL14", "GXYLT2", 
				 "TMTC4", "CBLN4", "PVALB", "TAC1", "MYO5B", "WIF1", "HHATL", "BCAN", "MYO1E", "BACH1", 
				 "FRMPD1", "SLC39A14", "C1QL1", "PTHLH", "MME", "PLEKHH2", "UNC5B", "GFRA1", "ZAK", 
				 "TNNT2", "PRDM8", "RMST", "RORB", "S100A6", "PTER", "PTGFRN", "SNTG2", "PLXND1", "MET", 
				 "RTKN2", "CLMN", "MYLIP", "LIG4", "CROT", "PLCG1", "DAPK2", "RMDN2", "ADAMTS3", "ACVR1C", 
				 "SLC24A4", "MTHFR", "NR4A2", "FAP", "FGF10", "OPRK1", "CCDC134", "MEI1", "AMZ1", 
				 "PLEKHA2", "MAPK12", "CACNA1H", "RGMA", "CHST8", "SCN4B", "AGPAT9", "APEH", "SYT2", 
				 "KCNN1", "HYOU1", "TRPC3", "SHB", "NXPH4", "NXPH3", "CTGF", "MAOB", "RRP1B", "OLIG2", 
				 "CSPG4", "ZCCHC24", "LRP4", "VCAN", "PDGFRA", "EMID1", "CALCRL", "AFAP1L2", "FGFR3", 
				 "GJA1", "F3", "MLC1", "ETNPPL", "GRIN2C", "AQP4", "PRDM16", "SLC25A18", "OAF", "PPP1R14A", 
				 "ERMN", "OPALIN", "MOG", "MAL", "TRIM59", "MAG", "MYRF", "FA2H", "GPR37", "XAF1", "ITIH5", 
				 "LEF1", "EMCN", "UACA", "TNS1", "APCDD1", "EBF1", "CGNL1", "ATP10A", "C1QC", "LST1", 
				 "C1QB", "FGD2", "CX3CR1", "TYROBP", "LTC4S", "LY86", "LAPTM5", "CD86")
  
human_plot <- group_heatmap_plot(data_source = humanFolder,
                   genes = common_genes,
                   group_by = "homology_cluster",
                   clusters = 1:37,
                   calculation = "median",
                   labelheight = 5,
                   showcounts = TRUE)
ggsave(paste0(outputFolder,"FigED13_common_genes_human.pdf"),human_plot,height = 20, width = 6)


 
mouse_plot <- group_heatmap_plot(data_source = mouseFolder,
                   genes = capitalize(tolower(common_genes)),
                   group_by = "homology_cluster",
                   clusters = 1:37,
                   calculation = "median",
                   labelheight = 5,
                   showcounts = TRUE)
ggsave(paste0(outputFolder,"FigED13_common_genes_mouse.pdf"),mouse_plot,height = 20, width = 6)


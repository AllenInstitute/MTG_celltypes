# This code builds the QC jitterplots

# FIRST NEED TO SET THE WORKING DIRECTORY TO THE MAIN FOLDER
# This is the folder where "Start.RData" is located

#-----------------------------------------------------------------------------------------------
print("Read in the data and load relevant libraries")
outputFolder  = "output/"
scriptsFolder = "R/"
inputFolder   = "data/"
dataFolder    = "MTG/"

library(dplyr)
library(ggplot2)
library(ggbeeswarm)

## Add a cluster order variable to the annotation feather for plotting

load(paste0(outputFolder,"dend.RData"))
clord <- 1:length(labels(dend))
names(clord) <- substr(labels(dend),3,4)
names(clord) <- gsub("_","",names(clord))
anno <- read_feather(paste0(dataFolder,"anno.feather"))
if(!file.exists(paste0(dataFolder,"anno_original.feather"))) 
  write_feather(anno,paste0(dataFolder,"anno_original.feather"))

anno$dendcluster_id    <- as.numeric(as.character(anno$cluster_id))
#anno$dendcluster_id[is.na(anno$dendcluster_id)] = max(anno$dendcluster_id)+5
anno$dendcluster_label <- paste0("cl",anno$cluster_id)
anno$dendcluster_color <- anno$cluster_color
file.remove(paste0(dataFolder,"anno.feather"))


## Read in and add the cluster separability

cc.mean <- read.csv(paste0(inputFolder,"mtg_cluster_separability.csv"))
anno$cc.mean.within_label <- cc.mean[match(anno$seq_name_label,cc.mean$sample_id),"cc.min.diff"]
write_feather(anno,paste0(dataFolder,"anno.feather"))


dendcluster_anno <- anno %>%
  select(dendcluster_id, dendcluster_label, dendcluster_color) %>%
  unique()

## Required functions
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE, roundall = F) {
  require(dplyr)
  # This does the summary. For each group return a vector with
  # N, mean, and sd
  
  names(data)[names(data) == measurevar] <- "measurevar"
  
  datac <- data %>%
    select(one_of(groupvars,"measurevar")) %>%
    filter(ifelse(na.rm == T, !is.na(measurevar), T)) %>%
    mutate(measurevar = as.numeric(measurevar)) %>%
    group_by_(c(groupvars)) %>%
    summarise(N = n(),
              median = median(measurevar),
              mean = mean(measurevar),
              max = max(measurevar),
              sd = ifelse(N == 1, 0, sd(measurevar)),
              q25 = as.numeric(quantile(measurevar, 0.25, na.rm=na.rm)),
              q75 = as.numeric(quantile(measurevar, 0.75, na.rm=na.rm))) %>%
    mutate(se = sd/sqrt(N))
  #%>%
  #  mutate(ci =  se * qt(conf.interval/2 + 0.5, N-1))
  
  
  if(roundall) {
    roundcols <- c("median","mean","max","sd","q25","q75","se","ci")
    datac[roundcols] <- round(datac[roundcols],3)
  }
  
  # datac <- datac %>%
  #   mutate(xpos = 1:n())
  
  return(datac)
}


# This function makes the actual plots!
qcPlot <- function(name,scaleLimits = c(-5000, 12000), 
                     scaleBreaks = seq(0, 12000, 2000),
                     scaleLabels = seq(0,12,2),
					 ylab = "value",
					 fileName = gsub("\\.","_",gsub("_label","",name)),
					 na.rm = FALSE)
{

# dendcluster_id is the annotation for cluster ordering based on the current, bootstrapped dendrogram
stats <- summarySE(data = anno,
                         measurevar = name,
                         groupvars = "dendcluster_id",
						 na.rm = na.rm)


genes_plot <- ggplot() +
  # geom_quasirandom from the ggbeeswarm package
  # makes violin-shaped jittered point plots
  geom_quasirandom(data = anno,
                   aes(x = dendcluster_id,
                       y = eval(parse(text=name)),  # might need eval(parse(text=name))
                       color = dendcluster_color),  # "skyblue"
                   # Need to set position_jitter height = 0 to prevent
                   # jitter on the y-axis, which changes data representation
                   position = position_jitter(width = .3,height = 0),
                   size = 0.1) +
  # Errorbars built using stats values
  geom_errorbar(data = stats,
                aes(x = dendcluster_id,
                    ymin = q25,
                    ymax = q75),
                size = 0.2) +
  # Median points from stats
  geom_point(data = stats,
             aes(x = dendcluster_id,
                 y = median),
             color = "red",
             size = 0.5) +
  # Cluster labels as text objects
  geom_text(data = dendcluster_anno,
            aes(x = dendcluster_id,
                y = 0,
                label = dendcluster_label,
                color = dendcluster_color),
            angle = 90,
            hjust = 2,
            vjust = 0.3,
            size = 2*5/6) +
  scale_color_identity() +
  # Expand the y scale so that the labels are visible
  scale_y_continuous(ylab,
                     limits = scaleLimits, 
                     breaks = scaleBreaks,
                     labels = scaleLabels) +
  # Remove X-axis title
  scale_x_continuous("") +
  theme_bw() +
  # Theme tuning
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggsave(paste0(outputFolder,fileName,"_QC.pdf"), genes_plot, width = 8, height = 2, useDingbats = F)

}



## Make the QC plots for numeric variables!

qcPlot("genes_detected_cpm_criterion_label",scaleLimits = c(-2000, 14000), scaleBreaks = seq(0, 14000, 2000), 
  scaleLabels = seq(0,14,2),ylab="Genes Detected (thousands)",fileName="FigED1d_genes_detected")
qcPlot("total_reads_label",scaleLimits = c(-2000000, 7000000), scaleBreaks = seq(0, 7000000, 1000000), 
  scaleLabels = seq(0,7,1),ylab="Total reads (millions)",fileName="FigED1c_total_reads")
qcPlot("cc.mean.within_label",scaleLimits = c(-0.25,1), scaleBreaks = seq(0, 1, 0.25), na.rm=TRUE, 
  scaleLabels = seq(0,1,0.25),ylab="Mean w/in cluster CC",fileName="FigED3a_cluster_separability")


  
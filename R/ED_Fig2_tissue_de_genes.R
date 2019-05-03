library(feather)
library(limma)

shiny.path <- "//molgen-shiny-server/share/shiny/human/MTG_paper_rev/"
anno <- as.data.frame(read_feather(paste0(shiny.path, "anno.feather")))
data <- as.data.frame(read_feather(paste0(shiny.path, "data_t.feather")))
row.names(data) <- data$gene
data <- log2(data[, -grep("gene", colnames(data))] + 1)

surg.donors <- c("H16.06.002", "H16.06.008", "H16.06.009", "H16.03.004")
pm.donors <- c("H200.1023", "H200.1025", "H200.1030", "H16.24.010")
female.donors <- c("H200.1023", "H16.06.002", "H16.06.008", "H16.06.009")
anno$donor_type <- factor(ifelse(anno$external_donor_name_label %in% surg.donors, 
                          "surg", "pm"))
anno$donor_sex <- factor(ifelse(anno$external_donor_name_label %in% female.donors, 
                                "female", "male"))
cl.cnt <- table(anno$cluster_label, anno$donor_type)
keep.cl <- row.names(cl.cnt)[apply(cl.cnt, 1, min) >= 10]
anno.subset <- droplevels(subset(anno, cluster_label %in% keep.cl))
keep.samp <- match(anno.subset$sample_id, colnames(data))
data.subset <- data[, keep.samp]


# Calc DE
cl <- make.names(anno.subset$cluster_label)

design <- model.matrix(~0 + anno.subset$donor_type + anno.subset$donor_sex + cl)
colnames(design) <- sub("anno.subset$", "", colnames(design), fixed = TRUE)
fit1 <- lmFit(data.subset, design = design)
cont.matrix <- makeContrasts(tissue_type = "donor.typesurg - donor.typepm", 
                             levels = design) 
fit2 <- contrasts.fit(fit1, cont.matrix)
fit3 <- eBayes(fit2)
top1 <- topTable(fit3, p.value = 0.01, lfc = 1, number = Inf)
top2 <- data.frame(gene = row.names(top1), top1)

# Save tables
de.surg <- subset(top2, logFC > 0)
de.pm <- subset(top2, logFC <= 0)
de.pm$logFC <- -de.pm$logFC
write.csv(de.surg, "output/surg_de_genes.csv", row.names = FALSE)
write.csv(de.pm, "output/pm_de_genes.csv", row.names = FALSE)



# QC metrics by tissue type
with(anno.subset, 
     data.frame(reads = tapply(total_reads_label, donor_type, median),
                aligned = tapply(star_percent_unique_reads_label, donor_type, median),
                genes = tapply(Genes.With.FPKM_label, donor_type, median))
)


# Use GEOquery to read the local series matrix file (GSE67705)

#______________loading packages__________________#
library(TCGAbiolinks)
install.packages()
library(GEOquery)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylation450kmanifest)
library(minfi)
library(limma)
library(missMethyl)
library(DMRcate)
library(Gviz)
library(ggplot2)
library(RColorBrewer)
library(edgeR)

#________________________Import Data____________________________#
options(timeout = 6000)  # 10 minutes

gse2 <- getGEO(filename = "D:/Project 42/Methylation Profiling by array/DMR/FINAL_RUN_JUNE_28_2025/GSE67705_series_matrix.txt.gz", GSEMatrix = TRUE)
#gsm_list <- gse[[1]]  # If only one GSE object
class(gse2) # Should have an expression set

###############____________________Beta value matrix & Metadata______________________#

# Expression (beta values), typically in the assay data or GSM tables
beta_vals_GSE6770 <- as.data.frame(exprs(gse2))
beta_vals_GSE67705 <- beta_vals_GSE6770 # Using a different variable
#beta_vals_GSE67705 <- beta_vals_GSE67705[ , !(names(beta_vals_GSE67705) %in% c("GSM1653324", "GSM1653325")) ]

# Metadata (clinical/phenotype data)
clinical_GSE67705 <- pData(gse2)

#__________________Inititial Cleaning____________________#

# Identify columns (samples) with NA in beta values
na_columns <- names(beta_vals_GSE67705)[colSums(is.na(beta_vals_GSE67705)) > 0]


# Remove those columns from beta values
beta_vals_GSE67705 <- beta_vals_GSE67705[ , !(names(beta_vals_GSE67705) %in% na_columns) ]

# Remove corresponding rows from clinical data
clinical_GSE67705 <- clinical_GSE67705[!(rownames(clinical_GSE67705) %in% na_columns), ]

#_________________Check if Row and Colnames are similar_____#
#Check if all row and col names are similar and correctly ordered
all(colnames(beta_vals_GSE67705) %in% rownames(clinical_GSE67705))  # Should return TRUE
all(rownames(clinical_GSE67705) %in% colnames(beta_vals_GSE67705))  # Should also return TRUE
all(colnames(beta_vals_GSE67705) == rownames(clinical_GSE67705))  # TRUE = already aligned

# If not correctly ordered match it
# Ensure sample names match
# colnames(beta_vals_GSE67705) <- rownames(clinical_GSE67705) #RUN IF FALSE

#_________________CLEANING Beta Values_____________________##
# Get annotation
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Remove probes with NA
probe.na_5 <- rowSums(is.na(beta_vals_GSE67705))
probe_5 <- probe.na_5[probe.na_5 == 0]
beta_vals_GSE67705 <- beta_vals_GSE67705[row.names(beta_vals_GSE67705) %in% names(probe_5), ]

# Remove chrX/Y
keep <- !(row.names(beta_vals_GSE67705) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
beta_vals_GSE67705 <- beta_vals_GSE67705[keep, ]

# Remove SNPs
no.snp.probe <- ann450k$Name[is.na(ann450k$Probe_rs)]
snp.probe <- ann450k[!is.na(ann450k$Probe_rs), ]
snp5.probe <- snp.probe$Name[snp.probe$Probe_maf <= 0.05]
beta_vals_GSE67705 <- beta_vals_GSE67705[row.names(beta_vals_GSE67705) %in% c(no.snp.probe, snp5.probe), ]

# Remove cross-reactive probes
crs.reac <- read.csv("cross_reactive_probe.chen2013.csv") #Cite this paper
crs.reac <- crs.reac$TargetID[-1]
beta_vals_GSE67705 <- beta_vals_GSE67705[ -which(row.names(beta_vals_GSE67705) %in% crs.reac), ]

#________________m_value Caluculation_________________#
# Save filtered beta matrix
bval_GSE67705 <- beta_vals_GSE67705

#Convert to m value
mval_GSE67705 <- t(apply(bval_GSE67705, 1, function(x) log2(x/(1 - x))))

#_______________ORGANIZING METADATA_________________#
# Choose appropriate columns and prepare groups
clinical_GSE67705 <- clinical_GSE67705[, c("geo_accession", "source_name_ch1")]  # modify as needed
colnames(clinical_GSE67705) <- c("Sample", "Condition")             # rename if needed
clinical_GSE67705$Condition <- gsub(".*:", "", clinical_GSE67705$Condition)  # extract label
clinical_GSE67705$Condition <- trimws(clinical_GSE67705$Condition)

# Remove Texts other than HIV+ or HIV- from each row
clinical_GSE67705$Condition <- gsub("^(HIV[+-]).*$", "\\1", clinical_GSE67705$Condition)

# Make sample names row names
rownames(clinical_GSE67705) <- clinical_GSE67705$Sample

#__________________ORGANIZING & DEFINING REFERENCE_____________#
# Match order
bval_GSE67705 <- bval_GSE67705[, rownames(clinical_GSE67705)]
mval_GSE67705 <- mval_GSE67705[, rownames(clinical_GSE67705)]

# Grouping variable
clinical_GSE67705$Condition <- as.factor(clinical_GSE67705$Condition)
clinical_GSE67705$Condition <- relevel(clinical_GSE67705$Condition, ref = "HIV-")  # modify as appropriate

#_______________________FIT LINEAR MODEL_____________#
design_GSE67705 <- model.matrix(~ Condition, data = clinical_GSE67705)
fit_GSE67705 <- lmFit(mval_GSE67705, design_GSE67705) #should have done with mval
fit2_GSE67705 <- eBayes(fit_GSE67705)

# Top Differentially Methylated Probes (DMPs)
dmp_results_GSE67705 <- topTable(fit2_GSE67705, coef = 2, number = Inf, adjust.method = "BY")
write.csv(dmp_results_GSE67705,"04_DMPs_GSE67705.csv")

#_____________Volcano Plot______________________#
# Define thresholds for significance
Hyper_ID <- dmp_results_GSE67705[dmp_results_GSE67705$logFC > 0.4 & dmp_results_GSE67705$adj.P.Val < 0.05, ]
Hypo_ID  <- dmp_results_GSE67705[dmp_results_GSE67705$logFC < -0.4 & dmp_results_GSE67705$adj.P.Val < 0.05, ]

#MANIPULATING
# NotSig <- 
# write.csv(Hyper_ID,"04_Hyper_GSE67705.csv")
# write.csv(Hypo_ID,"04_Hypo_GSE67705.csv")

# Create combined data with methylation status
combined_data_adj <- rbind(
  transform(Hyper_ID, Methylation = "Hyper-methylated"),
  transform(Hypo_ID, Methylation = "Hypo-methylated"),
  transform(
    dmp_results_GSE67705[
      !(rownames(dmp_results_GSE67705) %in% c(rownames(Hyper_ID), rownames(Hypo_ID))),
    ],
    Methylation = "Not Significant"
  )
)

# Ensure Methylation is treated as a factor
combined_data_adj$Methylation <- factor(combined_data_adj$Methylation,
                                        levels = c("Hyper-methylated", "Hypo-methylated", "Not Significant"))

# Set plot title (replace with your desired string)
plot_title <- "Volcano Plot of Differential Methylation"

# Create the volcano plot
library(ggplot2)
volcano_plot <- ggplot(combined_data_adj, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Methylation), size = 2, alpha = 0.7) +
  geom_vline(xintercept = c(-0.4, 0.4), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = plot_title,
    x = "Delta Beta Value",
    y = "-log10(Adjusted P-Value)",
    color = "Methylation"
  ) +
  scale_color_manual(values = c(
    "Hyper-methylated" = "#E36A5D",
    "Hypo-methylated" = "#722C6E",
    "Not Significant" = "grey"
  )) +
  theme_minimal()


# Updated volcano plot without grid and heading
# Volcano plot with no grid but with x and y axes
volcano_plot <- ggplot(combined_data_adj, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Methylation), size = 2, alpha = 0.7) +
  geom_vline(xintercept = c(-0.4, 0.4), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Delta Beta Value",
    y = "-log10(Adjusted P-Value)",
    color = "Methylation"
  ) +
  scale_color_manual(values = c(
    "Hyper-methylated" = "#E36A5D",
    "Hypo-methylated" = "#722C6E",
    "Not Significant" = "grey"
  )) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),       # Remove grid
    plot.title = element_blank(),       # Remove title
    axis.line = element_line(color = "black"), # Add axis lines
    axis.ticks = element_line(color = "black") # Show axis ticks
  )

# To display
print(volcano_plot)


# Print the plot
print(volcano_plot)


# Filter only Hyper- and Hypo-methylated probes
significant_data <- subset(combined_data_adj, Methylation != "Not Significant")

# Volcano plot without grey points
volcano_plot <- ggplot(significant_data, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Methylation), size = 2, alpha = 0.7) +
  geom_vline(xintercept = c(-0.4, 0.4), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Delta Beta Value",
    y = "-log10(Adjusted P-Value)",
    color = "Methylation"
  ) +
  scale_color_manual(values = c(
    "Hyper-methylated" = "#E36A5D",
    "Hypo-methylated" = "#722C6E"
  )) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# Display the plot
print(volcano_plot)


#_____________ANNOTATION OF DMPs___________#

annotate_dmps <- function(dmp_results, annotation_df, output_file = NULL) {
  # Ensure probe IDs are a column
  dmp_results$Probe_ID <- rownames(dmp_results)
  
  # Check required columns in annotation
  required_cols <- c("Name", "chr", "pos", "strand", 
                     "UCSC_RefGene_Name", "UCSC_RefGene_Group", 
                     "Relation_to_Island")
  missing_cols <- setdiff(required_cols, colnames(annotation_df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in annotation:", paste(missing_cols, collapse = ", ")))
  }
  
  # Subset annotation
  ann_subset <- annotation_df[, required_cols]
  
  # Merge DMPs with annotation
  dmp_annotated <- merge(dmp_results, ann_subset, by.x = "Probe_ID", by.y = "Name")
  
  # Sort by adjusted p-value
  if ("adj.P.Val" %in% colnames(dmp_annotated)) {
    dmp_annotated <- dmp_annotated[order(dmp_annotated$adj.P.Val), ]
  }
  
  # Adjust the chromosome labels to be numeric (e.g., chr1, chr2 -> 1, 2)
  dmp_annotated$chr_numeric <- as.numeric(gsub("chr", "", dmp_annotated$chr))
  
  #bp to Mbp
  dmp_annotated$pos_mbp <- (dmp_annotated$pos)/1000000
  
  # Optionally save to file
  if (!is.null(output_file)) {
    write.csv(dmp_annotated, output_file, row.names = FALSE)
  }
  
  return(dmp_annotated)
}


# Annotate and save results---USING FUNCTION
annotated_dmps <- annotate_dmps(dmp_results_GSE67705, ann450k, "DMPs_ann_mann_GSE67705.csv")

check_ann_dmp <- as.data.frame(annotated_dmps)


#______DMR FINDING______________#

# setting some annotation
myAnnotation <- cpg.annotate(object = mval_GSE67705,
                             datatype = "array", 
                             what = "M", 
                             analysis.type = "differential", 
                             design = design_GSE67705, 
                             contrasts = FALSE, 
                             coef = "ConditionHIV+", 
                             arraytype = "450K",
                             fdr = 0.01)

str(myAnnotation)

# DMR analysis
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges

# DMR CSV File
dmr.table <- data.frame(results.ranges)
write.csv(dmr.table,"001_DMRs_GSE67705.csv")
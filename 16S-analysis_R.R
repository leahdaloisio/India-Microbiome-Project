################################################################################
################################   R STUDIO    #################################
################################################################################

# Load required packages 
library(BiocManager)
library(phyloseq)
library(ANCOMBC)
library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library (vegan)
library(dunn.test)
library(qiime2R)
library(stats)  # for cmdscale function

# Use feature table derived using QIIME2 to construct your OTU table in R
otu <- read.table(file = "feature-table-gg2.tsv", 
                  sep = "\t", header = T, row.names = 1, 
                  skip = 0, comment.char = "")

# Use taxonomy table derived using QIIME2 to construct your taxonomy table in R
taxonomy <- read.table(file = "taxonomy-gg2.tsv", sep = "\t", 
                       header = T ,row.names = 1)

# Clean the taxonomy table -- make sure the tidyverse package is loaded for this
tax <- taxonomy %>%
  select(Taxon) %>%
  separate(Taxon, c("Kingdom", "Phylum", "Class", "
                    Order", "Family", "Genus", "Species"), "; ")

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "k__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)

tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean=="__"] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], 
                                  sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = " ")
  }
}

################################################################################
############################ CREATE PHYLOSEQ FILE ##############################
################################################################################

# Read in metadata for ps object
metadata <- read.csv("sample-metadata_V5.csv")

# Handling empty sample names by assigning them a unique identifier
metadata$sample_name[metadata$sample_name == ""] <- paste("EmptyName", seq_along(metadata$sample_name[metadata$sample_name == ""]), sep="_")

# Handling duplicates by ensuring unique sample names
metadata$sample_name <- make.unique(metadata$sample_name)

# Set row names
rownames(metadata) <- metadata$sample_name
metadata <- metadata[, !names(metadata) %in% "sample_name"]

# Convert to sample_data
SAMPLE <- sample_data(metadata)

# Read in OTU table for ps object
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE) 

# Read in Taxonomy table for ps object
TAX = tax_table(as.matrix(tax.clean))
# Re-order the tax table to be the same order as OTU (need to have them matched to create ps)
TAX <- TAX[match(row.names(OTU), row.names(TAX)), ]

# Read in tree for ps object
TREE = read_tree("tree-gg2.nwk")

# Merge the data into a file called "ps"
ps <- phyloseq(OTU, TAX, SAMPLE, TREE)

################################################################################
########################## ALPHA DIVERSITY ANALYSIS ############################
################################################################################

# Ensure that 'group' is a factor and set the levels in the desired order
sample_data(ps)$group <- factor(sample_data(ps)$group, levels = c("Indian", "Indo-Immigr", "Indo-Can", "Euro-Can", "Euro-Immigr"))

### Plot Shannon
plot_shannon <- plot_richness(ps, x = "group", measures = c("Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    plot.title = element_blank()
  ) +
  labs(y = "Shannon Index", title = "")

print(plot_shannon)

### Plot Pielou's Evenness
# Calculate Pielou's Evenness
richness_estimates <- estimate_richness(ps, measures = c("Shannon"), split = TRUE)
species_richness <- estimate_richness(ps, measures = c("Observed"), split = TRUE)
richness_estimates$Pielou <- richness_estimates$Shannon / log(species_richness$Observed)

# Convert to long format for ggplot2
df_pielou <- data.frame(sample_data(ps), Pielou = richness_estimates$Pielou)

# Now plot Pielou's Evenness
plot_pielou <- ggplot(df_pielou, aes(x = group, y = Pielou)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(y = "Pielou's Evenness", title = NULL)

print(plot_pielou)

## Both Shannon and Pielou's were tested with rarefied data as well, but no major differences were observed


################################################################################
################################ RAREFACTION ###################################
################################################################################

set.seed(111) # keep result reproductive
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=12053, replace=F) 

################################################################################
########################  BETA DIVERSITY WITH ANALYSIS  ########################
################################################################################

#### Using QIIME2 results here ####

####### Bray Curtis #######

# Import the qza file
distance_matrix_qza <- read_qza("core-metrics-results_gg2b/bray_curtis_distance_matrix.qza")

# Convert to a matrix
distance_matrix <- as.matrix(distance_matrix_qza$data)

# Compute classical MDS (equivalent to PCoA)
pcoa_result <- cmdscale(distance_matrix, eig = TRUE, k = 2)

# Convert to a data frame
pcoa_data <- as.data.frame(pcoa_result$points) %>%
  tibble::rownames_to_column("sample_name")

# Flip the Y-axis
pcoa_data$V2 <- pcoa_data$V2 * -1

# Merge metadata
metadata <- read.delim("sample-metadata_V5.txt", header = TRUE)
combined_df <- left_join(metadata, pcoa_data, by = "sample_name")

# Custom RGB colors
custom_colors <- c(
  "Euro-Can" = rgb(255/255, 215/255, 0/255),
  "Indian" = rgb(3/255, 192/255, 74/255),
  "Euro-Immigr" = rgb(0/255, 204/255, 204/255),
  "Indo-Can" = rgb(242/255, 80/255, 34/255),  
  "Indo-Immigr" = rgb(191/255, 64/255, 191/255)  
)

# Placeholder percentages of variation
percent_variation_axis1 <- "21.38%"
percent_variation_axis2 <- "7.235%"
percent_variation_axis3 <- "4.523%"

# Plot using ggplot
p <- ggplot(combined_df, aes(x = V1, y = V2, color = group)) +
  geom_point() +
  stat_ellipse(aes(fill = group), geom = "polygon", level = 0.95, alpha = 0.1) + 
  labs(
    title = "Bray-Curtis PCoA with Ellipses", 
    x = paste("PCoA 1 (", percent_variation_axis1, ")", sep = ""),
    y = paste("PCoA 2 (", percent_variation_axis2, ")", sep = "")
  ) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5)
  )

print(p)


####### Weighted Unifrac #######

# Load Weighted UniFrac distance matrix
weighted_unifrac_qza <- read_qza("core-metrics-results_gg2/weighted_unifrac_distance_matrix.qza")
weighted_unifrac_matrix <- weighted_unifrac_qza$data

# PCoA analysis
pcoa_results2 <- vegan::capscale(weighted_unifrac_matrix ~ 1)

# Extract PCoA data
pcoa_data2 <- pcoa_results2$data$Vectors %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_name")

pcoa_data2 <- as.data.frame(pcoa_results2$CA$u) %>%
  rownames_to_column(var = "sample_name") %>%
  select(sample_name, MDS1, MDS2)

# Load metadata
metadata <- read.delim("sample-metadata_V5.txt", header = TRUE)

# Join metadata with PCoA data
combined_df2 <- left_join(metadata, pcoa_data2, by = "sample_name")

# Custom RGB colors
custom_colors <- c(
  "Euro-Can" = rgb(255/255, 215/255, 0/255),
  "Indian" = rgb(3/255, 192/255, 74/255),
  "Euro-Immigr" = rgb(0/255, 204/255, 204/255),
  "Indo-Can" = rgb(242/255, 80/255, 34/255),  
  "Indo-Immigr" = rgb(191/255, 64/255, 191/255)  
)

# Placeholder percentages of variation
percent_variation_axis1 <- sprintf("%.1f%%", 100 * pcoa_results2$CA$eig[1] / sum(pcoa_results2$CA$eig))
percent_variation_axis2 <- sprintf("%.1f%%", 100 * pcoa_results2$CA$eig[2] / sum(pcoa_results2$CA$eig))

# Plot using ggplot

p2 <- ggplot(combined_df2, aes(x = MDS1, y = -MDS2, color = group)) + # negative sign before MDS2 is to just flip the points (a preference)
  geom_point() +
  stat_ellipse(aes(fill = group), geom = "polygon", level = 0.95, alpha = 0.2) +
  labs(
    title = "Weighted UniFrac PCoA with Ellipses", 
    x = paste("PCoA 1 (", percent_variation_axis1, ")", sep = ""),
    y = paste("PCoA 2 (", percent_variation_axis2, ")", sep = "")
  ) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.5)
  )


print(p2)

################################################################################
#################################  LEFSE  ######################################
################################################################################

# Collapse the table.qza to a level (taxonomic level you want)
qiime taxa collapse \
--i-table table-filtered-taxa.qza  \
--o-collapsed-table table-L7.qza \
--p-level 7 \
--i-taxonomy taxonomy-gg2.qza

# Calculate relative frequency for the collapsed table (instead of counts you get relative abundance)
qiime feature-table relative-frequency \
--i-table table-L7.qza \
--o-relative-frequency-table lefse-table-L7.qza \
--output-dir lefse/
  
# Export biom file
qiime tools export \
--input-path lefse-table-L7.qza \
--output-path lefse-L7/
  
# Convert biom file to a text file (for LEFSE comparison)
biom convert \
-i lefse-L7/feature-table.biom \
-o lefse-table-L7-LAST.txt --header-key "taxonomy" --to-tsv

# One against all 
lefse_format_input.py lefse_table-4.txt lefse_table-4.in -c 1 -u 2 -o 1000000

# Run lefse analysis
lefse_run.py  lefse_table-4.in  lefse_table-4.res -l 3.5
# edit the names in the .in file to make it as you'd like .. import into excel and change the info to ".in" after

# Create visualization 
lefse_plot_res.py  NIJ-filtered-9.res  NIJ-filtered.png --left_space 0.4 --dpi 1200 

# Couldn't get the legend to appear for cladogram so did it in Galaxy instead

################################################################################
##############################  RDA PLOTS   ####################################
################################################################################ 

# Reference file: dbRDA-weighted-unifrac.R

################################################################################
##########################  SAMPLE SIZE CALCULATIONS  ##########################
################################################################################

############### Bray Curtis Distance for Sample size Calculation ###############

# Calculate the distance matrix
dist <- phyloseq::distance(ps.rarefied, method = "bray")

# Extract cohort information (replace "group" with your actual column name)
cohort <- sample_data(ps.rarefied)$group

# Convert distance matrix to a data frame
dist_df <- as.data.frame(as.matrix(dist))

# Add cohort information to the distance data frame
dist_df$cohort <- cohort

# Initialize vectors to store results
cohorts <- unique(cohort)
iqr_distance <- numeric(length(cohorts))
median_distance <- numeric(length(cohorts))

# Calculate median distance and IQR for each cohort
for (i in 1:length(cohorts)) {
  current_cohort <- cohorts[i]
  
# Extract distances for pairs of samples within the current cohort
cohort_samples <- which(cohort == current_cohort)
cohort_distances <- dist_df[cohort_samples, cohort_samples]
  
# Convert the upper triangle of the distance matrix to a vector
cohort_dist_vector <- as.vector(cohort_distances[upper.tri(cohort_distances)])
  
# Calculate median and IQR for the cohort
median_distance[i] <- median(cohort_dist_vector, na.rm = TRUE)
iqr_distance[i] <- IQR(cohort_dist_vector, na.rm = TRUE)
}

# Create a data frame for median and IQR results
summary_table <- data.frame(cohort = cohorts, median_distance = median_distance, IQR = iqr_distance)

# Print the results
print(summary_table)

#        cohort median_distance       IQR
# 1    Euro-Can       0.6218784 0.1154484
# 2      Indian       0.8273044 0.1767817
# 3 Indo-Immigr       0.6088526 0.1285987
# 4    Indo-Can       0.6330789 0.1206961
# 5 Euro-Immigr       0.6669709 0.1226251


## Indian vs. Euro-Canadians
# (0.8273044 - 0.6218784)  / 0.1767817
# 1.162032^2 = 1.350318
# (7.84/1.350318)*2
# = 11.6
# = 12 ppl per group

## Indian vs. Indo-Immigrants
# (0.8273044 - 0.6088526)  / 0.1767817
# 1.235715^2 = 1.526992
# (7.84/1.526992)*2
# = 10.3
# = 10 ppl per group

## Indian vs. Indo-Canadians
# (0.8273044 - 0.6330789) / 0.1767817
# 1.098674^2 = 1.207085
# (7.84/1.207085)*2
# = 12.9
# = 13 ppl per group

## Indo-Immigrants vs. Indo-Canadians
# (0.6330789 - 0.6088526) / 0.1285987
# 0.1883868^2 = 0.03548959
# (7.84/0.03548959)*2
# = 441.8
# = 442 ppl per group

## Indo-Immigrants vs. Euro-Canadians
# (0.6218784 - 0.6088526) / 0.1285987
# 0.1012903^2 = 0.01025972
# (7.84/0.01025972)*2
# = 1528.307
# = 1528 ppl per group





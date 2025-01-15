First, taxonomic classification was conducted using MetaPhlAn4 using script from: https://github.com/biobakery/biobakery/wiki/metaphlan4
MetaPhlAn taxonomic profile table output was then imported into R for further analysis. 

################################################################################
###############################   R STUDIO  ####################################
################################################################################

library(phyloseq)
library(ggplot2)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ape)
library(dplyr)
library(tidyr)
library(pheatmap)

# Load data
# "#CLASS" row removed in table
data <- read.csv("MPA-raw-6.csv", header=TRUE, row.names=1)

# Convert to phyloseq object (assuming OTU table for simplicity)
otu_table <- as.matrix(data)

taxa_table <- data.frame(Taxonomy = rownames(data))
rownames(taxa_table) <- rownames(data)

# Clean the taxonomy table -- make sure the tidyverse package is loaded for this
tax <- taxa_table %>%
  select(Taxonomy) %>%
  separate(Taxonomy, c("k__", "p__", "c__", "
                    o__", "f__", "g__", "s__"), "\\|")

# Add in sample metadata 
metadata <- read.csv("sample-metadata_V6.csv", row.names = 1)

# Create new phyloseq with cleaned taxonomy table 
ps <- phyloseq(otu_table(otu_table, taxa_are_rows = TRUE), tax_table(as.matrix(tax)), sample_data(metadata))

# Tree 

tree = read.tree('mpa_vJan21_CHOCOPhlAnSGB_202103.nwk')

tree_meta = read.delim('mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt', header=F)

tree_meta$V1=gsub('SGB|_group','',tree_meta$V1)

tree_meta=separate(tree_meta,V2,into=c('V2','Other1','Other2','Other3','Other4','Other5','Other6'),sep=',')

tree_meta=tree_meta[!duplicated(tree_meta[c('V2')]), ]

tree_meta2=subset(tree_meta,V1 %in% tree$tip.label)
tree2=drop.tip(tree,setdiff(tree$tip.label,tree_meta2$V1))
tree_meta2=tree_meta2[match(tree2$tip.label,tree_meta2$V1),]
identical(tree2$tip.label,tree_meta2$V1)
tree2$tip.label=tree_meta2$V2

phy_tree(ps)=tree2

ps

# Export the tree to use for Microbiome Analyst
phy_tree <- phy_tree(ps)
write.tree(phy_tree, file = "output_tree.nwk")


############################## BETA DIVERSITY ##################################

# NOTE: rarefaction not applied for shotgun sequence data because MetaPhlAn4 output is already in relative abundances

############# Bray Curtis ##############

dist_matrix <- distance(ps, method="bray")
pcoa_res <- ordinate(ps, method="PCoA", distance="bray")

# Custom colors 
custom_colors <- c(
  "Euro-Can" = rgb(255/255, 215/255, 0/255),
  "Indian" = rgb(3/255, 192/255, 74/255),
  "Euro-Immigr" = rgb(0/255, 204/255, 204/255),
  "Indo-Can" = rgb(242/255, 80/255, 34/255),  
  "Indo-Immigr" = rgb(191/255, 64/255, 191/255)  
)

### Bray-Curtis Analysis ###

# Calculate Bray-Curtis distance and perform PCoA
dist_matrix_bray <- phyloseq::distance(ps, method = "bray")
pcoa_res_bray <- ordinate(ps, method = "PCoA", distance = dist_matrix_bray)

# Plotting Bray-Curtis PCoA
p_bray <- plot_ordination(ps, pcoa_res_bray, color = "group") +
  theme_minimal() +
  labs(title = "PCoA using Bray-Curtis Distance") +
  scale_color_manual(values = custom_colors) +
  stat_ellipse(aes(fill = group), geom = "polygon", level = 0.95, alpha = 0.1) +
  geom_point(size = 2)+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black") 
  )

print(p_bray)


# Pairwise PERMANOVA for Bray-Curtis
sample_data <- data.frame(sample_data(ps))
dist_matrix_bray <- as.matrix(dist_matrix_bray)
results_bray <- data.frame()
group_combinations <- combn(unique(sample_data$group), 2, simplify = FALSE)

for(combination in group_combinations) {
  group1 <- combination[1]
  group2 <- combination[2]
  
  subset_indices <- which(sample_data$group %in% c(group1, group2))
  sub_dist_matrix <- dist_matrix_bray[subset_indices, subset_indices]
  sub_sample_data <- sample_data[subset_indices, , drop = FALSE]
  
  permanova_res <- adonis(sub_dist_matrix ~ group, data = sub_sample_data, permutations = 999)
  
  pseudo_F <- permanova_res$aov.tab$'F.Model'[1]
  p_value <- permanova_res$aov.tab$'Pr(>F)'[1]
  
  results_bray <- rbind(results_bray, data.frame(Group1 = group1, Group2 = group2, 
                                                 pseudo_F = pseudo_F, p_value = p_value))
}

results_bray$p_adjusted <- p.adjust(results_bray$p_value, method = "bonferroni")
print(results_bray)
#Group1      Group2  pseudo_F p_value p_adjusted
#1  Indo-Immigr    Indo-Can  4.215003   0.001       0.01
#2  Indo-Immigr    Euro-Can  8.976951   0.001       0.01
#3  Indo-Immigr Euro-Immigr  5.313757   0.001       0.01
#4  Indo-Immigr      Indian  9.132465   0.001       0.01
#5     Indo-Can    Euro-Can  3.084292   0.001       0.01
#6     Indo-Can Euro-Immigr  2.166793   0.001       0.01
#7     Indo-Can      Indian  8.660676   0.001       0.01
#8     Euro-Can Euro-Immigr  1.163757   0.199       1.00
#9     Euro-Can      Indian 18.284343   0.001       0.01
#10 Euro-Immigr      Indian 11.182855   0.001       0.01


### Weighted UniFrac Analysis ###

# Prune tree and merge with phyloseq object
tree <- read_tree("mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt")
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, rownames(otu_table(ps))))
physeq <- merge_phyloseq(ps, pruned_tree)

# Calculate Weighted UniFrac distance and perform PCoA
dist_matrix_unifrac <- phyloseq::distance(ps, method = "wunifrac")
pcoa_res_unifrac <- ordinate(ps, method = "PCoA", distance = dist_matrix_unifrac)

# Plotting Weighted UniFrac PCoA
p_unifrac <- plot_ordination(ps, pcoa_res_unifrac, color = "group") +
  theme_minimal() +
  labs(title = "PCoA using Weighted-Unifrac Distance") +
  scale_color_manual(values = custom_colors) +
  stat_ellipse(aes(fill = group), geom = "polygon", level = 0.95, alpha = 0.1) +
  geom_point(size = 2) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black") 
  )

print(p_unifrac)


# Pairwise PERMANOVA for Weighted UniFrac
sample_data <- data.frame(sample_data(ps))
dist_matrix_unifrac <- as.matrix(dist_matrix_unifrac)
results_unifrac <- data.frame()

for(combination in group_combinations) {
  group1 <- combination[1]
  group2 <- combination[2]
  
  subset_indices <- which(sample_data$group %in% c(group1, group2))
  sub_dist_matrix <- dist_matrix_unifrac[subset_indices, subset_indices]
  sub_sample_data <- sample_data[subset_indices, , drop = FALSE]
  
  permanova_res <- adonis(sub_dist_matrix ~ group, data = sub_sample_data, permutations = 999)
  
  pseudo_F <- permanova_res$aov.tab$'F.Model'[1]
  p_value <- permanova_res$aov.tab$'Pr(>F)'[1]
  
  results_unifrac <- rbind(results_unifrac, data.frame(Group1 = group1, Group2 = group2, 
                                                       pseudo_F = pseudo_F, p_value = p_value))
}

results_unifrac$p_adjusted <- p.adjust(results_unifrac$p_value, method = "bonferroni")
print(results_unifrac)
#      Group1      Group2  pseudo_F p_value p_adjusted
# 1  Indo-Immigr    Indo-Can  4.670629   0.001       0.01
# 2  Indo-Immigr    Euro-Can 11.932016   0.001       0.01
# 3  Indo-Immigr Euro-Immigr  6.903393   0.001       0.01
# 4  Indo-Immigr      Indian 25.984365   0.001       0.01
# 5     Indo-Can    Euro-Can  5.409267   0.001       0.01
# 6     Indo-Can Euro-Immigr  2.360739   0.010       0.10
# 7     Indo-Can      Indian 25.206253   0.001       0.01
# 8     Euro-Can Euro-Immigr  1.116192   0.348       1.00
# 9     Euro-Can      Indian 49.611215   0.001       0.01
# 10 Euro-Immigr      Indian 31.565399   0.001       0.01



################################################################################
################## ------ PREVOTELLA CLADE ABUNDANCES--------- #################
################################################################################

###### CALCULATING MEAN ABUNDANCE OF P.COPRI CLADES

# Step 1: Load the Data with suppressed warnings
data <- suppressWarnings(read.csv("MPA-prevotella.txt", sep="\t", row.names = 1, fill = TRUE, stringsAsFactors = FALSE))

# Step 2: Make Column Names Unique
unique_names <- make.unique(colnames(data))
colnames(data) <- unique_names

# Step 3: Verify and Group Columns
# Adjust the patterns based on the actual column names
cohort_groups <- list(
  Indian = grep("^Indian", unique_names, ignore.case = TRUE),
  Indo_Immigr = grep("^Indo.Immigr", unique_names, ignore.case = TRUE),
  Indo_Can = grep("^Indo.Can", unique_names, ignore.case = TRUE),
  Euro_Can = grep("^Euro.Can", unique_names, ignore.case = TRUE),
  Euro_Immigr = grep("^Euro.Immigr", unique_names, ignore.case = TRUE)
)

# Verify if groups are correct
print(sapply(cohort_groups, function(x) unique_names[x]))

# Step 4: Calculate the average abundance for each cohort
data_avg <- sapply(cohort_groups, function(cols) rowMeans(data[, cols], na.rm = TRUE))

# Convert to a data frame and set appropriate row and column names
data_avg <- as.data.frame(data_avg)
colnames(data_avg) <- c("Indian", "Indo_Immigr", "Indo_Can", "Euro_Can", "Euro_Immigr")

# Step 5: Transpose the data for heatmap (flip the orientation)
data_transposed <- data_avg

# Step 6: Create the Heatmap with adjusted parameters
pheatmap(data_transposed, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Abundances of Prevotella copri Clades in Different Cohorts",
         fontsize_row = 10, fontsize_col = 10, angle_col = "45")


# Step 7: Print the average abundances for reporting
print("Average abundances of Prevotella copri clades in different cohorts:")
print(data_avg)











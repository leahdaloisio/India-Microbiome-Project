################################################################################
#######################-#-#-#-#-# - RDA PLOTS - #-#-#-#-#-######################
################################################################################

setwd("/Users/leahdaloisio/Desktop/16S-analysis")

# Load packages
library(qiime2R)
library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(ggrepel)

# Import data
metadata <- read.csv(file = "sample-metadata_V5.csv")

# Remove the last row of metadata
metadata <- metadata[-nrow(metadata),]

baseline_var <- c("age", "alcohol", "processed_food_percentage", "vegetarian", "pescetarian", "fiber_1000kcal")

# Convert variables to appropriate data types
metadata$sample_name <- as.factor(metadata$sample_name)
metadata$group <- as.factor(metadata$group)
metadata$alcohol <- as.factor(metadata$alcohol)
metadata$processed_food_percentage <- as.numeric(metadata$processed_food_percentage)
metadata$pescetarian <- as.factor(metadata$pescetarian)
metadata$vegetarian <- as.factor(metadata$vegetarian)
metadata$fiber_1000kcal <- as.numeric(metadata$fiber_1000kcal)

# Read the distance matrix data
distance_matrix_data <- qiime2R::read_qza("weighted_unifrac_distance_matrix.qza")$data

# Convert to matrix
distance_matrix <- as.matrix(distance_matrix_data)

# Extract the sample names from the distance matrix
distance_matrix_samples <- rownames(distance_matrix)

# Print out the row names for debugging
print(head(distance_matrix_samples))

# Find the names in metadata that match the distance matrix
matching_samples <- metadata$sample_name %in% distance_matrix_samples

# Run these so the sample names match
metadata$sample_name <- as.character(metadata$sample_name)
distance_matrix_samples <- as.character(distance_matrix_samples)

# Import distance matrix QIIME 2
wunifrac_dist <- qiime2R::read_qza("weighted_unifrac_distance_matrix.qza")$data %>% 
  as.matrix() %>%
  .[metadata$sample_name, metadata$sample_name] %>% 
  as.dist()

### Check for colinearity ###

baseline_model <- paste(baseline_var, collapse="+")

# Run Distance-Based Redundancy Analysis (i.e. "dbrda")
set.seed(2022)
wunifrac_dbrda_base <- vegan::dbrda(as.formula(paste("wunifrac_dist ~", baseline_model)),
                                data=metadata, scale=TRUE, na.action = na.exclude, parallel=6)
vegan::vif.cca(wunifrac_dbrda_base) #checks for collinearity, values <10 are ok.

# ANOVA
set.seed(2022)
(bray_dbrda_base_anova2 <- vegan::anova.cca(wunifrac_dbrda_base, parallel=6, step=1000)) 

# Calculating R^2 value
bray_dbrda_base_adjR2 <- vegan::RsquareAdj(wunifrac_dbrda_base)$adj.r.squared %>% 
  round(3)

print(wunifrac_dbrda_base_adjR2) # R^2 value of 17%


# Number of permutations: 999

#Model: vegan::dbrda(formula = wunifrac_dist ~ age + alcohol + processed_food_percentage + vegetarian + pescetarian + fiber_1000kcal, data = metadata, na.action = na.exclude, scale = TRUE, parallel = 6)
#Df SumOfSqs      F Pr(>F)    
#Model      6   1.8728 5.9923  0.001 ***
#  Residual 165   8.5946                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Triplot

perc <- round(100*(summary(wunifrac_dbrda_base)$cont$importance[2, 1:2]), 2) # axis explaining %

sc_si <- scores(wunifrac_dbrda_base, display="sites", choices=c(1,2), scaling=2) %>% # grab group coordinates
  as_tibble(rownames="sites") %>% drop_na(dbRDA1)

sc_bp <- scores(wunifrac_dbrda_base, display="bp", choices=c(1, 2), scaling=2) %>% # grab environmental coordinates
  as_tibble(rownames="vars")

sc_bp_num <- sc_bp %>% filter(vars %in% c("age", "processed_food_percentage", "fiber_1000kcal")) %>%
  mutate(vars=str_replace_all(vars, c("age" = "age","processed_food_percentage" = "ultra-processed food", "fiber_1000kcal" = "fiber"))) 

sc_bp_cat <- sc_bp %>% 
  filter(!vars %in% c("age", "processed_food_percentage", "alcohol", "vegetarian", "pescetarian", "fiber_1000kcal")) %>%
  mutate(vars=str_replace_all(vars, c("alcohol1"="alcohol", "vegetarian1" = "vegetarian", 
                                      "pescetarian1" = "pescetarian")))


# Set colours for cohorts
custom_colors <- c(
  "Euro-Can" = rgb(255/255, 215/255, 0/255),
  "Indian" = rgb(3/255, 192/255, 74/255),
  "Euro-Immigr" = rgb(0/255, 204/255, 204/255),
  "Indo-Can" = rgb(242/255, 80/255, 34/255),  
  "Indo-Immigr" = rgb(191/255, 64/255, 191/255)  
)

# Merge the dataframes
sc_si2 <- merge(sc_si, metadata[, c("sample_name", "group")], 
                by.x="sites", by.y="sample_name", all.x=TRUE)

# Plot RDA
dbrda_triplot <- (
  ggplot(sc_si2, aes(x=dbRDA1, y=dbRDA2)) + 
    labs(x=paste("dbRDA-1 (", perc[1], "%)", sep=""), 
         y=paste("dbRDA-2 (", perc[2], "%)", sep="")) +
    
    # Add horizontal and vertical reference lines
    geom_hline(yintercept = 0, color = "#666666", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "#666666", linetype = "dashed") +
    
    geom_point(aes(fill = group), color = "NA", shape=21, alpha=1, size=3) + 
    scale_fill_manual(values=custom_colors, guide=guide_legend(override.aes=list(size=4))) +
    
    # Environmental Vectors 
    geom_segment(
      data= dplyr::mutate(sc_bp_cat, dbRDA1 = dbRDA1 * 2, dbRDA2 = dbRDA2 * 2),
      size=1,
      color="steelblue",
      inherit.aes=F,
      show.legend=F,
      aes(x=0, xend=dbRDA1, y=0, yend=dbRDA2), 
      arrow = arrow(length = unit(0.3,"cm"))
    ) +
    
    # Labels for Environmental Vectors
    geom_label_repel(
      data=sc_bp_cat, 
      aes(label=vars, x=dbRDA1, y=dbRDA2),
      box.padding = 2.5,
      max.overlaps = Inf,
      fill="steelblue",
      color="white",
      segment.color = "steelblue",
      show.legend=FALSE,
      force = 10,
      max.iter = 5000,
      direction = "both"
    ) +
    
    # Environmental Vectors for Numeric Variables
    geom_segment(
      data= dplyr::mutate(sc_bp_num, dbRDA1 = dbRDA1 * 2, dbRDA2 = dbRDA2 * 2),
      size=1,
      color="steelblue",
      inherit.aes=F,
      show.legend=F,
      aes(x=0, xend=dbRDA1, y=0, yend=dbRDA2), 
      arrow = arrow(length = unit(0.3,"cm"))
    ) +
    
    # Labels for Environmental Vectors for Numeric Variables
    geom_label_repel(
      data=sc_bp_num, 
      aes(label=vars, x=dbRDA1, y=dbRDA2),
      box.padding = 2.5,
      max.overlaps = Inf,
      fill="steelblue",
      color="white",
      segment.color = "steelblue",
      show.legend=FALSE,
      force = 10,
      nudge_y = ifelse(sc_bp_cat$vars == "pescetarian", 0.5, 0),
      max.iter = 5000,
      direction = "both"
    ) +
    
    theme_bw() +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(colour = "#333333"),
      axis.line.y = element_line(colour = "#333333"),
      legend.text = element_text(size = 8)
    )
)

print(dbrda_triplot)


### If you're interested to get individual adjusted p values

# Test the significance of each term with permutation tests
wunifrac_dbrda_base_anova_terms <- vegan::anova.cca(wunifrac_dbrda_base, by = "term", permutations = 999)
print(wunifrac_dbrda_base_anova_terms)

# Extract p-values from the result of anova.cca
p_values <- wunifrac_dbrda_base_anova_terms$"Pr(>F)"

# Adjust p-values for multiple testing using FDR correction
adjusted_p_values <- p.adjust(p_values, method = "fdr")
print(adjusted_p_values)

#              Variable Adjusted_p_value
#                 age        0.0440
#             alcohol        0.0030
# processed_food_percentage  0.0030
#          vegetarian        0.3710
#         pescetarian        0.0705
#               fiber        0.2268



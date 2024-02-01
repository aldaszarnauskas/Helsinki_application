# Author: Aldas Žarnauskas
# Title: Analysis of Taxonomic Abundance in Infant Samples
# Date: 21/007/2023
# Last time updated: 21/007/2023
# Objective: investigate if previous antibiotic exposure affects alpha diversity and
# microbial composition of 12 month old infant samples

### ------------------------- Cleaning Environment -------------------------- ##
################################################################################
rm(list = ls(all.names = T))
gc()

### ------------------------------- Libraries ------------------------------- ##
################################################################################
library(tidyverse)
library(phyloseq)
library(vegan)
library('knitr')

set.seed(199402)

### ----------------------- Data Import + Wrangling ------------------------- ##
################################################################################
# Loading metadata and taxonomic abundance table
load("candidate_assignment.Rdata.rdata")
# Rounding taxonomic abundance values rounded upwards to the nearest integer
taxtable <- apply(taxtable %>% as.matrix(), 2, ceiling) %>% t()
# Convert taxtable into the otu_table object
otu <- otu_table(taxtable, taxa_are_rows=T)

# Set sample IDs to rownames
metadata <- meta %>% column_to_rownames(var = "sample_name") 

#--------------------------------- Antibiotic Exposure ------------------------#
# Select rows where Age is 1 month
metadata_1 <- metadata %>% .[.$Age == "1 month",] %>% 
  {.$Antibiotic_1 <- .$Antibiotic; .} %>% 
  dplyr::select(Participant_id, Antibiotic_1)
# Select rows where Age is 6 months
metadata_6 <- metadata %>% .[.$Age == "6 months",] %>% 
  {.$Antibiotic_6 <- .$Antibiotic; .} %>% 
  dplyr::select(Participant_id, Antibiotic_6)
# Select rows where Age is 12 months
metadata_12 <- metadata %>% .[.$Age == "12 months",] %>% 
  {.$Antibiotic_12 <- .$Antibiotic; .} %>% 
  dplyr::select(Participant_id, Antibiotic_12)

# Quality inspection. Check if participant_id is present in above three objects
all(metadata_1$Participant_id == metadata_6$Participant_id)
all(metadata_1$Participant_id %in% metadata_6$Participant_id)
all(metadata_1$Participant_id == metadata_12$Participant_id)
all(metadata_1$Participant_id %in% metadata_6$Participant_id)

# Join all three data frames to have patient's antibiotic exposure history on a single row 
metadata_join <- full_join(metadata_1, metadata_6, by = "Participant_id") %>% 
  full_join(., metadata_12, by = "Participant_id") %>% 
  {.$antibiotic_exposure <- ifelse(.$Antibiotic_1 == 1, "Exposed_at_1m", ifelse(
    .$Antibiotic_6 == 1, "Exposed_at_6m", ifelse(
      .$Antibiotic_12 == 1, "Exposed_at_12m", "No_exposure")
  )); .} %>% 
  dplyr::select(Participant_id, antibiotic_exposure)

# Add new columns showing patients antibiotic exposure to a metadata
metadata_new <- metadata %>% left_join(., metadata_join, by = "Participant_id") %>%
  {rownames(.) <- rownames(metadata); .} %>% 
  sample_data() # convert into sample_data object

###---------------------- Phyloseq Construction and Analysis ----------------###
################################################################################
physeq <- phyloseq(otu, metadata_new)
rm(taxtable)
# Subset rows where Age = 12 months
physeq_12m <- subset_samples(physeq, Age == "12 months")

# Estimate samples' alpha diversity
richness <- estimate_richness(physeq_12m)

# Plot selected richness metrics
plot_richness(physeq_12m, x = "antibiotic_exposure",
              measures = c("Observed", "Chao1", "Shannon"),
              color = "antibiotic_exposure") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
  
# Transform data into relative abundance
physeq_12m.transf <- transform_sample_counts(physeq_12m, function(x) x/sum(x))
# Perform PCoA using Bray-Curtis dissimilarities
ord.meas <- ordinate(physeq_12m.transf, method = "PCoA", distance = "bray")
# Plot a PCoA of Bray-Curtis dissimilarity meassures
p <- plot_ordination(physeq_12m.transf, ord.meas, color = "antibiotic_exposure") +
  geom_point(size=2) + 
  theme_classic() + 
  theme(text = element_text(size=18), axis.text = element_text(size=16), 
        legend.position = "right") +
  ggtitle("PCoA: Bray-Curtis")

###------------------------------- PERMANOVA: adonis2 -----------------------###
################################################################################

#------------------------------------- All Groups PERMANOVA -------------------#
# Calculate Bray-Curtis dissimilarities for infant samples
dist.meas <- phyloseq::distance(physeq_12m.transf, method = "bray")

# Use PERMANOVA to indentify if at least one group differes in Bray-Curtis dissimilarities
bray <- adonis2(dist.meas ~ 
                       sample_data(physeq_12m.transf)$antibiotic_exposure, 
                     permutations = 999, method = "bray")
# Use PERMANOVA to indentify if at least one group differes in observed diversity
Observed <- adonis2(richness[["Observed"]] ~ sample_data(physeq_12m.transf)$antibiotic_exposure,
                 data = metadata,
                 permutations = 999)
# Use PERMANOVA to indentify if at least one group differes in Chao1 richness
Chao1 <- adonis2(richness[["Chao1"]] ~ sample_data(physeq_12m.transf)$antibiotic_exposure,
        data = metadata,
        permutations = 999)
# Use PERMANOVA to indentify if at least one group differes in Shannon diversity
Shannon <- adonis2(richness[["Shannon"]] ~ sample_data(physeq_12m.transf)$antibiotic_exposure,
                 data = metadata,
                 permutations = 999)

#-------------------------------- Pairwise Comparison -------------------------#
# Select rows where antibiotic_exposure is Ex1m or Ex6m
Ex1m_Ex6m <- sample_data(physeq_12m.transf) %>% as.data.frame() %>% 
  .[.$antibiotic_exposure == "Exposed_at_1m" | .$antibiotic_exposure == "Exposed_at_6m",]
# Select rows where antibiotic_exposure is Ex1m or Ex12m
Ex1m_Ex12m <- sample_data(physeq_12m.transf) %>% as.data.frame() %>%
  .[.$antibiotic_exposure == "Exposed_at_1m" | .$antibiotic_exposure == "Exposed_at_12m",]
# Select rows where antibiotic_exposure is Ex1m or ExNo
Ex1m_ExNo <- sample_data(physeq_12m.transf) %>% as.data.frame() %>%
  .[.$antibiotic_exposure == "Exposed_at_1m" | .$antibiotic_exposure == "No_exposure",]
# Select rows where antibiotic_exposure is Ex6m or Ex12m
Ex6m_Ex12m <- sample_data(physeq_12m.transf) %>% as.data.frame() %>%
  .[.$antibiotic_exposure == "Exposed_at_6m" | .$antibiotic_exposure == "Exposed_at_12m",]
# Select rows where antibiotic_exposure is Ex6m or ExNo
Ex6m_ExNo <- sample_data(physeq_12m.transf) %>% as.data.frame() %>%
  .[.$antibiotic_exposure == "Exposed_at_6m" | .$antibiotic_exposure == "No_exposure",]
# Select rows where antibiotic_exposure is Ex12m or ExNo
Ex12m_ExNo <- sample_data(physeq_12m.transf) %>% as.data.frame() %>%
  .[.$antibiotic_exposure == "Exposed_at_12m" | .$antibiotic_exposure == "No_exposure",]

# Store F values from the pairwise comparisons
pairwise_adonis <- numeric()

# Store the pair names
pairs <- ls()[4:9]

# pairwise_comparison performs pairwise analysis of variance
pairwise_comparison <- function(value_type, values, pair_list){
  
  # Stores pairwise comparison F values
  pairwise_adonis <- numeric()
  
  # Iteratively performs pairwise comparison for each pair
  for (pair in pair_list){
    # Stores the value
    name <- paste(value_type, pair, sep = "_")
    pair <- get(pair)
    
    if (value_type == "Bray"){
      
      values <- values %>% as.matrix()
      pairwise_adonis[name] <- adonis2(values[rownames(
        pair),rownames(pair)] ~ pair$antibiotic_exposure, permutations = 999) %>% 
        .[["Pr(>F)"]] 
      
    } else{
      
      pairwise_adonis[name] <- adonis2(values[value_type][rownames(
        pair),] ~ pair$antibiotic_exposure, permutations = 999) %>% 
        .[["Pr(>F)"]] 
      
    }
  } 

  return(pairwise_adonis)
  
}

Chao1_pairwise <- pairwise_comparison("Chao1", richness, pairs)
Bray_pairwise <- pairwise_comparison("Bray", dist.meas, pairs)
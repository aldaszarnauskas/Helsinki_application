rm(list = ls(all.names = T))
gc()

library(tidyverse)
library(phyloseq)
library(vegan)

set.seed(199402)

load("candidate_assignment.Rdata.rdata")
taxtable <- apply(taxtable %>% as.matrix(), 2, ceiling) %>% t()
otu <- otu_table(taxtable, taxa_are_rows=T)

metadata <- meta %>% column_to_rownames(var = "sample_name") 
metadata_filt <- metadata %>% 
  {.$Age_months <- ifelse(.$Age == "1 month", 1,
                          ifelse(.$Age == "6 months", 6, 12));.} %>% 
  sample_data()

#------------------------- Antibiotic Exposure
metadata_1 <- metadata %>% .[.$Age == "1 month",] %>% 
  {.$Antibiotic_1 <- .$Antibiotic; .} %>% 
  dplyr::select(Participant_id, Antibiotic_1)
metadata_6 <- metadata %>% .[.$Age == "6 months",] %>% 
  {.$Antibiotic_6 <- .$Antibiotic; .} %>% 
  dplyr::select(Participant_id, Antibiotic_6)
metadata_12 <- metadata %>% .[.$Age == "12 months",] %>% 
  {.$Antibiotic_12 <- .$Antibiotic; .} %>% 
  dplyr::select(Participant_id, Antibiotic_12)

metadata_1$Participant_id == metadata_6$Participant_id
metadata_1$Participant_id == metadata_12$Participant_id

metadata_join <- full_join(metadata_1, metadata_6, by = "Participant_id") %>% 
  full_join(., metadata_12, by = "Participant_id") %>% 
  {.$antibiotic_exposure <- ifelse(.$Antibiotic_1 == 1, "Exposed_at_1m", ifelse(
    .$Antibiotic_6 == 1, "Exposed_at_6m", ifelse(
      .$Antibiotic_12 == 1, "Exposed_at_12m", "No_exposure")
  )); .} %>% 
  dplyr::select(Participant_id, antibiotic_exposure)

metadata_new <- metadata %>% left_join(., metadata_join, by = "Participant_id") %>%
  {rownames(.) <- rownames(metadata); .} %>% 
  sample_data()

#------------------------ Phyloseq object construction and analysis
physeq <- phyloseq(otu, metadata_new)
physeq_12m <- subset_samples(physeq, Age == "12 months")

richness <- estimate_richness(physeq_12m)

plot_richness(physeq_12m, x = "antibiotic_exposure"
              measures = c("Observed", "Chao1", "Shannon"),
              color = "antibiotic_exposure") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
  


physeq_12m.transf <- transform_sample_counts(physeq_12m, function(x) x/sum(x))
ord.meas <- ordinate(physeq_12m.transf, method = "PCoA", distance = "bray")
p <- plot_ordination(physeq_12m.transf, ord.meas, color = "antibiotic_exposure") +
  geom_point(size=2) + 
  theme_classic() + 
  theme(text = element_text(size=18), axis.text = element_text(size=16), 
        legend.position = "right") +
  ggtitle("PCoA: Bray-Curtis")

#------------------------ Adonis
library('knitr')



dist.meas <- phyloseq::distance(physeq_12m.transf, method = "bray")
bray <- adonis2(dist.meas ~ 
                       sample_data(physeq_12m.transf)$antibiotic_exposure, 
                     permutations = 999, method = "bray")
adonis.res





# Observed, Chao1 (richness), and Shannon (diversity)
Observed <- adonis2(richness[["Observed"]] ~ sample_data(physeq_12m.transf)$antibiotic_exposure,
                 data = metadata,
                 permutations = 999)
Chao1 <- adonis2(richness[["Chao1"]] ~ sample_data(physeq_12m.transf)$antibiotic_exposure,
        data = metadata,
        permutations = 999)
Shannon <- adonis2(richness[["Shannon"]] ~ sample_data(physeq_12m.transf)$antibiotic_exposure,
                 data = metadata,
                 permutations = 999)


pairwise_adonis <- numeric()

# Ex1m vs Ex6m
Ex1m_Ex6m <- sample_data(physeq_12m.transf) %>% as.data.frame() %>% 
  .[.$antibiotic_exposure == "Exposed_at_1m" | .$antibiotic_exposure == "Exposed_at_6m",]
# Ex1m vs Ex12m
Ex1m_Ex12m <- sample_data(physeq_12m.transf) %>% as.data.frame() %>%
  .[.$antibiotic_exposure == "Exposed_at_1m" | .$antibiotic_exposure == "Exposed_at_12m",]
# Ex1m vs ExNo
Ex1m_ExNo <- sample_data(physeq_12m.transf) %>% as.data.frame() %>%
  .[.$antibiotic_exposure == "Exposed_at_1m" | .$antibiotic_exposure == "No_exposure",]
# Ex6m vs Ex12m
Ex6m_Ex12m <- sample_data(physeq_12m.transf) %>% as.data.frame() %>%
  .[.$antibiotic_exposure == "Exposed_at_6m" | .$antibiotic_exposure == "Exposed_at_12m",]
# Ex6m vs ExNo
Ex6m_ExNo <- sample_data(physeq_12m.transf) %>% as.data.frame() %>%
  .[.$antibiotic_exposure == "Exposed_at_6m" | .$antibiotic_exposure == "No_exposure",]
# Ex12m vs ExNo
Ex12m_ExNo <- sample_data(physeq_12m.transf) %>% as.data.frame() %>%
  .[.$antibiotic_exposure == "Exposed_at_12m" | .$antibiotic_exposure == "No_exposure",]


pairwise_adonis["Chao1_Ex1m_Ex6m"] <- adonis2(
  richness["Chao1"][rownames(Ex1m_Ex6m),] ~ Ex1m_Ex6m$antibiotic_exposure,
        permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Chao1_Ex1m_Ex12m"] <- adonis2(
  richness["Chao1"][rownames(Ex1m_Ex12m),] ~ Ex1m_Ex12m$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Chao1_Ex1m_ExNo"] <- adonis2(
  richness["Chao1"][rownames(Ex1m_ExNo),] ~ Ex1m_ExNo$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Chao1_Ex6m_Ex12m"] <- adonis2(
  richness["Chao1"][rownames(Ex6m_Ex12m),] ~ Ex6m_Ex12m$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Chao1_Ex6m_ExNo"] <- adonis2(
  richness["Chao1"][rownames(Ex6m_ExNo),] ~ Ex6m_ExNo$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Chao1_Ex12m_ExNo"] <- adonis2(
  richness["Chao1"][rownames(Ex12m_ExNo),] ~ Ex12m_ExNo$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]




pairwise_adonis["Observed_Ex1m_Ex6m"] <- adonis2(
  richness["Observed"][rownames(Ex1m_Ex6m),] ~ Ex1m_Ex6m$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Observed_Ex1m_Ex12m"] <- adonis2(
  richness["Observed"][rownames(Ex1m_Ex12m),] ~ Ex1m_Ex12m$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Observed_Ex1m_ExNo"] <- adonis2(
  richness["Observed"][rownames(Ex1m_ExNo),] ~ Ex1m_ExNo$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Observed_Ex6m_Ex12m"] <- adonis2(
  richness["Observed"][rownames(Ex6m_Ex12m),] ~ Ex6m_Ex12m$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Observed_Ex6m_ExNo"] <- adonis2(
  richness["Observed"][rownames(Ex6m_ExNo),] ~ Ex6m_ExNo$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Observed_Ex12m_ExNo"] <- adonis2(
  richness["Observed"][rownames(Ex12m_ExNo),] ~ Ex12m_ExNo$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]




# Shannon
pairwise_adonis["Shannon_Ex1m_Ex6m"] <- adonis2(
  richness["Shannon"][rownames(Ex1m_Ex6m),] ~ Ex1m_Ex6m$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Shannon_Ex1m_Ex12m"] <- adonis2(
  richness["Shannon"][rownames(Ex1m_Ex12m),] ~ Ex1m_Ex12m$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Shannon_Ex1m_ExNo"] <- adonis2(
  richness["Shannon"][rownames(Ex1m_ExNo),] ~ Ex1m_ExNo$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Shannon_Ex6m_Ex12m"] <- adonis2(
  richness["Shannon"][rownames(Ex6m_Ex12m),] ~ Ex6m_Ex12m$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Shannon_Ex6m_ExNo"] <- adonis2(
  richness["Shannon"][rownames(Ex6m_ExNo),] ~ Ex6m_ExNo$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Shannon_Ex12m_ExNo"] <- adonis2(
  richness["Shannon"][rownames(Ex12m_ExNo),] ~ Ex12m_ExNo$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]


dist.meas <- dist.meas %>% as.matrix()
# Bray
bray <- adonis2(dist.meas ~ 
                  sample_data(physeq_12m.transf)$antibiotic_exposure, 
                permutations = 999, method = "bray")

pairwise_adonis["Bray_Ex1m_Ex6m"] <- adonis2(
  dist.meas[rownames(Ex1m_Ex6m),rownames(Ex1m_Ex6m)] ~ Ex1m_Ex6m$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Bray_Ex1m_Ex12m"] <- adonis2(
  richness["Shannon"][rownames(Ex1m_Ex12m),] ~ Ex1m_Ex12m$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Bray_Ex1m_ExNo"] <- adonis2(
  richness["Shannon"][rownames(Ex1m_ExNo),] ~ Ex1m_ExNo$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Bray_Ex6m_Ex12m"] <- adonis2(
  richness["Shannon"][rownames(Ex6m_Ex12m),] ~ Ex6m_Ex12m$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Bray_Ex6m_ExNo"] <- adonis2(
  richness["Shannon"][rownames(Ex6m_ExNo),] ~ Ex6m_ExNo$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

pairwise_adonis["Bray_Ex12m_ExNo"] <- adonis2(
  richness["Shannon"][rownames(Ex12m_ExNo),] ~ Ex12m_ExNo$antibiotic_exposure,
  permutations = 999) %>% .[["Pr(>F)"]]

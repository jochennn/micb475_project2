##### load libraries #####
library(readr)
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

# load fish phyloseq 
load("./create_phyloseq_obj/outputs/fish_phyloseq.RData")
fish_phyloseq_wdiv <- data.frame(sample_data(fish_phyloseq), estimate_richness(fish_phyloseq))

# PERMANOVA 
# weighted unifrac
dm_wunifrac <- UniFrac(fish_phyloseq, weighted=TRUE)

adonis2(dm_wunifrac ~ anatomical_location*spatial_region, data=fish_phyloseq_wdiv)

adonis2(dm_wunifrac ~ sample_type*habitat_depth_level3, data=fish_phyloseq_wdiv)

adonis2(dm_wunifrac ~ habitat_depth_level3, data=fish_phyloseq_wdiv)

adonis2(dm_wunifrac ~ sample_type, data=fish_phyloseq_wdiv)

# unweighted unifrac
dm_unifrac <- UniFrac(fish_phyloseq, weighted=FALSE)

adonis2(dm_unifrac ~ anatomical_location*spatial_region, data=fish_phyloseq_wdiv)

adonis2(dm_unifrac ~ sample_type*habitat_depth_level3, data=fish_phyloseq_wdiv)

adonis2(dm_unifrac ~ habitat_depth_level3, data=fish_phyloseq_wdiv)

adonis2(dm_unifrac ~ sample_type, data=fish_phyloseq_wdiv)

# bray
dm_bray <- vegdist(t(otu_table(fish_phyloseq)), method="bray")

adonis2(dm_bray ~ anatomical_location*spatial_region, data=fish_phyloseq_wdiv)

adonis2(dm_bray ~ sample_type*habitat_depth_level3, data=fish_phyloseq_wdiv)

adonis2(dm_bray ~ habitat_depth_level3, data=fish_phyloseq_wdiv)

# jaccard

dm_jaccard <- vegdist(t(otu_table(fish_phyloseq)), method="jaccard")

adonis2(dm_jaccard ~ anatomical_location*spatial_region, data=fish_phyloseq_wdiv)

adonis2(dm_jaccard ~ sample_type*habitat_depth_level3, data=fish_phyloseq_wdiv)

adonis2(dm_jaccard ~ habitat_depth_level3, data=fish_phyloseq_wdiv)




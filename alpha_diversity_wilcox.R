#### Load libraries ####
library(tidyverse)
library(phyloseq)

#### Load phyloseq object data ####
load("create_phyloseq_obj/outputs/fish_phyloseq.Rdata")

alphadiv <- estimate_richness(fish_phyloseq)
samp_dat <- sample_data(fish_phyloseq)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

#### Wilcox rank sum exact test ####
wilcox.test(Shannon ~ anatomical_location, data=samp_dat_wdiv)
wilcox.test(Shannon ~ spatial_region, data=samp_dat_wdiv)



## load packages
library(phyloseq)
library(indicspecies)
library(ggVennDiagram)
library(tidyverse)

## load phyloseq object
load("./create_phyloseq_obj/outputs/fish_phyloseq.Rdata")

## Convert to relative abundance
fish_phyloseq_genus <- 
  tax_glom(fish_phyloseq, "Genus", NArm = FALSE)

fish_phyloseq_genus_RA <- 
  transform_sample_counts(fish_phyloseq_genus, fun=function(x) x/sum(x))

## Indicator Species Analysis
isa_fish <- 
  multipatt(
    t(otu_table(fish_phyloseq_genus_RA)),
    cluster = sample_data(fish_phyloseq_genus_RA)$anat_space_combine
    )

summary(isa_fish)

taxtable <- 
  tax_table(fish_phyloseq) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_fish$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View()











# load packages
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

# Convert to relative abundance
fish_phyloseq_RA <- 
  transform_sample_counts(fish_phyloseq, fun=function(x) x/sum(x))

# Filter dataset by "anat_space_combine"
c_e_phyloseq <- 
  subset_samples(
    fish_phyloseq_RA, 
    anat_space_combine=="coastal/external")

c_i_phyloseq <- 
  subset_samples(
    fish_phyloseq_RA, 
    anat_space_combine=="coastal/internal")

o_e_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_space_combine=="open ocean/external")

o_i_phyloseq <-
  subset_samples(
    fish_phyloseq_RA, 
    anat_space_combine=="open ocean/internal")

# What ASVs are found in more than 50% of samples in each category?
c_e_core <- core_members(c_e_phyloseq, detection=0.001, prevalence=0.5)
c_i_core <- core_members(c_i_phyloseq, detection=0.001, prevalence=0.5)
o_e_core <- core_members(o_e_phyloseq, detection=0.001, prevalence=0.5)
o_i_core <- core_members(o_i_phyloseq, detection=0.001, prevalence=0.5)

# venn diagram
ggVennDiagram(x=list(c_e_core, c_i_core, o_e_core, o_i_core), 
              category.names=list(
                "c/e", 
                "c/i", 
                "o/e", 
                "o/i"
                ),
              label="count")


# What are these ASVs?
prune_taxa(o_e_core,fish_phyloseq) %>%
  plot_bar(, fill="Genus") +
    facet_wrap(.~anat_space_combine, scales="free")

prune_taxa(o_i_core,fish_phyloseq) %>%
  plot_bar(, fill="Genus") +
  facet_wrap(.~anat_space_combine, scales="free")





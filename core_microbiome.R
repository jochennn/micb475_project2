library(phyloseq)
library(microbiome)


phyloseq_object <- 
  transform_sample_counts(phyloseq_rawabundance, fun=function(x) x/sum(x))


group1_phyloseq <- subset_samples(phyloseq_object, Group=="group1")
group2_phyloseq <- subset_samples(phyloseq_object, Group=="group2")


group1_core <- core_members(group1_phyloseq, detection=0, prevalence=0.8)
group2_core <- core_members(group2_phyloseq, detection=0, prevalence=0.8)

library(ggVennDiagram)
ggVennDiagram(x=list(group1_core, group2_core))















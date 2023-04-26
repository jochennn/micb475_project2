#### Load libraries ####
library(ggplot2)
library(phyloseq)
library(ggpubr)
# Note: ggpubr package is used for running Wilxocon test on faceted graph at the end

#### all alpha diversity metrics ####
plot_richness(fish_phyloseq)

#### Selecting alpha diversity metrics ####
plot_richness(fish_phyloseq, measures = c("Shannon", "Observed", "Chao1", "ACE"))

#### Add ggplot layers for spatial region + Wilcoxon rank sum test ####
plot_richness(fish_phyloseq, x = "spatial_region", measures = c("Shannon", "Observed", "Chao1", "ACE")) +
  xlab("Spatial region") +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test")

gg_richness_shannonobservedchao1ace_spatial_region <- plot_richness(fish_phyloseq, x = "spatial_region", measures = c("Shannon", "Observed", "Chao1", "ACE")) +
  xlab("Spatial region") +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test")

gg_richness_shannonobservedchao1ace_spatial_region 

ggsave(filename = "alpha_diversity/ggplot_shannonobservedchao1ace_spatial_region.png"
       , gg_richness_shannonobservedchao1ace_spatial_region
       , height=6, width=8)

#### Add ggplot layers for anatomical location + Wilcoxon rank sum test ####
plot_richness(fish_phyloseq, x = "anatomical_location", measures = c("Shannon", "Observed", "Chao1", "ACE")) +
  xlab("Anatomical sample location") +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test")

gg_richness_shannonobservedchao1ace_anatomical_location <- plot_richness(fish_phyloseq, x = "anatomical_location", measures = c("Shannon", "Observed", "Chao1", "ACE")) +
  xlab("Anatomical sample location") +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test")

gg_richness_shannonobservedchao1ace_anatomical_location 

ggsave(filename = "alpha_diversity/ggplot_shannonobservedchao1ace_anatomical_location.png"
       , gg_richness_shannonobservedchao1ace_anatomical_location
       , height=6, width=8)

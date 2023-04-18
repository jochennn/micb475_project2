#### Load libraries ####
library(ggplot2)
library(phyloseq)
library(ggpubr)
# Note: ggpubr package is used for running Wilxocon test on faceted graph at the end

## load phyloseq object
load("./create_phyloseq_obj/outputs/fish_phyloseq.Rdata")

#### all alpha diversity metrics ####
plot_richness(fish_phyloseq)

#### select certain alpha diversity metrics ####
plot_richness(fish_phyloseq, measures = c("Shannon"))

#### Add ggplot layers for spatial region ####
gg_richness_spatial_region <- 
  plot_richness(fish_phyloseq, 
                x = "spatial_region", 
                measures = c("Shannon")
                ) +
  xlab("Spatial region") +
  geom_boxplot()

ggsave(filename = "./alpha_diversity/ggplot_spatial_region.png"
       , gg_richness_spatial_region
       , height=5, width=4)

#### Add ggplot layers for anatomical location ####
gg_richness_anatomical_location <- 
  plot_richness(fish_phyloseq, 
                x = "anatomical_location", 
                measures = c("Shannon")
                ) +
  xlab("Anatomical location") +
  geom_boxplot()

ggsave(filename = "./alpha_diversity/ggplot_anatomical_location.png"
       , gg_richness_anatomical_location
       , height=5, width=4)

#### Faceting by anatomical location ####
gg_richness_anatomical_location_faceted <- 
  plot_richness(fish_phyloseq, 
                x = "spatial_region", 
                measures = c("Shannon")
                ) +
  xlab("Spatial region") +
  geom_boxplot() +
  facet_grid(. ~ anatomical_location)

ggsave(filename = "./alpha_diversity/ggplot_anatomical_location_faceted.png"
       , gg_richness_anatomical_location_faceted
       , height=4, width=6)

#### Faceting by spatial region ####
gg_richness_spatial_region_faceted <- 
  plot_richness(fish_phyloseq, 
                x = "anatomical_location", 
                measures = c("Shannon")
                ) +
  xlab("Anatomical sample location") +
  geom_boxplot() + 
  facet_grid(. ~ spatial_region)

ggsave(filename = "./alpha_diversity/ggplot_spatial_region_faceted.png"
       , gg_richness_spatial_region_faceted
       , height=4, width=6)

#### Wilcox rank sum test on graph using stat_compare_means() ####
gg_richness_sr_wilcox <- 
  gg_richness_spatial_region + 
  stat_compare_means(method = "wilcox.test")

gg_richness_al_wilcox <- 
  gg_richness_anatomical_location + 
  stat_compare_means(method = "wilcox.test")

gg_richness_srf_wilcox <- 
  gg_richness_spatial_region_faceted + 
  stat_compare_means(method = "wilcox.test")

gg_richness_alf_wilcox <- 
  gg_richness_anatomical_location_faceted + 
  stat_compare_means(method = "wilcox.test")


#### Save images as png with wilcox p-values ####
ggsave(filename = "./alpha_diversity/ggplot_spatial_region_wilcox.png"
       , gg_richness_sr_wilcox
       , height=4, width=6)

ggsave(filename = "./alpha_diversity/ggplot_anatomical_location_wilcox.png"
       , gg_richness_al_wilcox
       , height=4, width=6)

ggsave(filename = "./alpha_diversity/ggplot_spatial_region_faceted_wilcox.png"
       , gg_richness_srf_wilcox
       , height=4, width=6)

ggsave(filename = "./alpha_diversity/ggplot_anatomical_location_faceted_wilcox.png"
       , gg_richness_alf_wilcox
       , height=4, width=6)

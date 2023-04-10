## load packages
library(phyloseq)
library(indicspecies)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)

## load phyloseq object
load("./create_phyloseq_obj/outputs/fish_phyloseq.Rdata")

## Convert to relative abundance, glom to genus
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
  tax_table(fish_phyloseq) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="ASV")

fish_table <- isa_fish$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05)

View(fish_table)

### Visualising data ###

#coastal external#
fish_table$Coastal_External <- rep(NA, nrow(fish_table))

for (i in 1:nrow(fish_table)) {
  if (fish_table$`s.coastal/external`[i] == 1) {
    fish_table$Coastal_External[i] <- fish_table$Genus[i]}
  else {fish_table$Coastal_External[i] <- "NA"}
}

for (i in 1:nrow(fish_table)) {
  if (is.na(fish_table$Coastal_External[i]) == TRUE|
      fish_table$Coastal_External[i] == "g__uncultured") {
    fish_table$Coastal_External[i] <- fish_table$Family[i]
  }
  if (is.na(fish_table$Coastal_External[i]) == TRUE|
      fish_table$Coastal_External[i] == "f__uncultured"|
      fish_table$Coastal_External[i] == "f__Unknown_Family") {
    fish_table$Coastal_External[i] <- fish_table$Order[i]
  }
  if (is.na(fish_table$Coastal_External[i]) == TRUE|
      fish_table$Coastal_External[i] == "o__uncultured") {
    fish_table$Coastal_External[i] <- fish_table$Class[i]
  }
  if (is.na(fish_table$Coastal_External[i]) == TRUE|
      fish_table$Coastal_External[i] == "c__uncultured") {
    fish_table$Coastal_External[i] <- fish_table$Phylum[i]
  } 
  if (is.na(fish_table$Coastal_External[i]) == TRUE|
      fish_table$Coastal_External[i] == "p__uncultured") {
    fish_table$Coastal_External[i] <- fish_table$Domain[i]
  }
}

#coastal internal#
fish_table$Coastal_Internal <- rep(NA, nrow(fish_table))

for (i in 1:nrow(fish_table)) {
  if (fish_table$`s.coastal/internal`[i] == 1) {
    fish_table$Coastal_Internal[i] <- fish_table$Genus[i]}
  else {fish_table$Coastal_Internal[i] <- "NA"}
}

for (i in 1:nrow(fish_table)) {
  if (is.na(fish_table$Coastal_Internal[i]) == TRUE|
      fish_table$Coastal_Internal[i] == "g__uncultured") {
    fish_table$Coastal_Internal[i] <- fish_table$Family[i]
  }
  if (is.na(fish_table$Coastal_Internal[i]) == TRUE|
      fish_table$Coastal_Internal[i] == "f__uncultured"|
      fish_table$Coastal_Internal[i] == "f__Unknown_Family") {
    fish_table$Coastal_Internal[i] <- fish_table$Order[i]
  }
  if (is.na(fish_table$Coastal_Internal[i]) == TRUE|
      fish_table$Coastal_Internal[i] == "o__uncultured") {
    fish_table$Coastal_Internal[i] <- fish_table$Class[i]
  }
  if (is.na(fish_table$Coastal_Internal[i]) == TRUE|
      fish_table$Coastal_Internal[i] == "c__uncultured") {
    fish_table$Coastal_Internal[i] <- fish_table$Phylum[i]
  } 
  if (is.na(fish_table$Coastal_Internal[i]) == TRUE|
      fish_table$Coastal_Internal[i] == "p__uncultured") {
    fish_table$Coastal_Internal[i] <- fish_table$Domain[i]
  }
}

#ocean external#
fish_table$Ocean_External <- rep(NA, nrow(fish_table))

for (i in 1:nrow(fish_table)) {
  if (fish_table$`s.open ocean/external`[i] == 1) {
    fish_table$Ocean_External[i] <- fish_table$Genus[i]}
  else {fish_table$Ocean_External[i] <- "NA"}
}

for (i in 1:nrow(fish_table)) {
  if (is.na(fish_table$Ocean_External[i]) == TRUE|
      fish_table$Ocean_External[i] == "g__uncultured") {
    fish_table$Ocean_External[i] <- fish_table$Family[i]
  }
  if (is.na(fish_table$Ocean_External[i]) == TRUE|
      fish_table$Ocean_External[i] == "f__uncultured"|
      fish_table$Ocean_External[i] == "f__Unknown_Family") {
    fish_table$Ocean_External[i] <- fish_table$Order[i]
  }
  if (is.na(fish_table$Ocean_External[i]) == TRUE|
      fish_table$Ocean_External[i] == "o__uncultured") {
    fish_table$Ocean_External[i] <- fish_table$Class[i]
  }
  if (is.na(fish_table$Ocean_External[i]) == TRUE|
      fish_table$Ocean_External[i] == "c__uncultured") {
    fish_table$Ocean_External[i] <- fish_table$Phylum[i]
  } 
  if (is.na(fish_table$Ocean_External[i]) == TRUE|
      fish_table$Ocean_External[i] == "p__uncultured") {
    fish_table$Ocean_External[i] <- fish_table$Domain[i]
  }
}

#ocean internal#
fish_table$Ocean_Internal <- rep(NA, nrow(fish_table))

for (i in 1:nrow(fish_table)) {
  if (fish_table$`s.open ocean/internal`[i] == 1) {
    fish_table$Ocean_Internal[i] <- fish_table$Genus[i]}
  else {fish_table$Ocean_Internal[i] <- "NA"}
}

for (i in 1:nrow(fish_table)) {
  if (is.na(fish_table$Ocean_Internal[i]) == TRUE|
      fish_table$Ocean_Internal[i] == "g__uncultured") {
    fish_table$Ocean_Internal[i] <- fish_table$Family[i]
  }
  if (is.na(fish_table$Ocean_Internal[i]) == TRUE|
      fish_table$Ocean_Internal[i] == "f__uncultured"|
      fish_table$Ocean_Internal[i] == "f__Unknown_Family") {
    fish_table$Ocean_Internal[i] <- fish_table$Order[i]
  }
  if (is.na(fish_table$Ocean_Internal[i]) == TRUE|
      fish_table$Ocean_Internal[i] == "o__uncultured") {
    fish_table$Ocean_Internal[i] <- fish_table$Class[i]
  }
  if (is.na(fish_table$Ocean_Internal[i]) == TRUE|
      fish_table$Ocean_Internal[i] == "c__uncultured") {
    fish_table$Ocean_Internal[i] <- fish_table$Phylum[i]
  } 
  if (is.na(fish_table$Ocean_Internal[i]) == TRUE|
      fish_table$Ocean_Internal[i] == "p__uncultured") {
    fish_table$Ocean_Internal[i] <- fish_table$Domain[i]
  }
}

# subseting and viewing new table
fish_table_view <- 
  subset(fish_table, 
         select = c(
           Coastal_External, 
           Coastal_Internal, 
           Ocean_External, 
           Ocean_Internal,
           index,
           p.value
         )
  )


View(fish_table_view)

## formatting the table for easier viewing ##

CE <- data.frame(
  x = "Coastal_External", 
  y = fish_table_view$Coastal_External,
  z = fish_table_view$index,
  w = fish_table_view$p.value
  )

CI <- data.frame(
  x = "Coastal_Internal", 
  y = fish_table_view$Coastal_Internal,
  z = fish_table_view$index,
  w = fish_table_view$p.value
)

OE <- data.frame(
  x = "Ocean_External", 
  y = fish_table_view$Ocean_External,
  z = fish_table_view$index,
  w = fish_table_view$p.value
)

OI <- data.frame(
  x = "Ocean_Internal", 
  y = fish_table_view$Ocean_Internal,
  z = fish_table_view$index,
  w = fish_table_view$p.value
)

fish_table2 <- data.frame(
  anat_space_combine = c(OE$x, OI$x, CE$x, CI$x),
  taxon = c(OE$y, OI$y, CE$y, CI$y),
  index = c(OE$z, OI$z, CE$z, CI$z),
  p.value = c(OE$w, OI$w, CE$w, CI$w)
)

fish_table2$taxon <- 
  ifelse(fish_table2$taxon == "NA", NA, fish_table2$taxon)

fish_table3 <- 
  fish_table2[complete.cases(fish_table2$taxon), ] %>%
  separate(col = taxon, into = c("level", "taxon"), sep = "__")

for (i in 1:nrow(fish_table3)) {
  if (fish_table3$level[i] == "g") {fish_table3$level[i] <- "Genus"}
  if (fish_table3$level[i] == "f") {fish_table3$level[i] <- "Family"}
  if (fish_table3$level[i] == "o") {fish_table3$level[i] <- "Order"}
  if (fish_table3$level[i] == "c") {fish_table3$level[i] <- "Class"}
  if (fish_table3$level[i] == "p") {fish_table3$level[i] <- "Phylum"}
  if (fish_table3$level[i] == "d") {fish_table3$level[i] <- "Domain"}
  }

View(fish_table3)

## Saving data as .tsv ##
write.table(
  fish_table3, 
  file = "./indicator_species/outputs/fish_indic_spp.tsv", 
  sep = "\t", 
  row.names = FALSE
)

## Saving Data as PDF ##
# Set the font size and style
theme_set(theme_gray(base_size = 3, base_family = "Arial"))

#OceanExternal
OceanInternalTable <- 
  subset(fish_table3, anat_space_combine == "Ocean_Internal") %>% 
  arrange(desc(index)) %>%
  subset(select = -c(anat_space_combine)) 

pdf("./indicator_species/outputs/OceanInternalTable.pdf", width = 7, height = 32)
grid.table(OceanInternalTable)
dev.off()

#OceanExternal
OceanExternalTable <- 
  subset(fish_table3, anat_space_combine == "Ocean_External") %>% 
  arrange(desc(index)) %>%
  subset(select = -c(anat_space_combine))

pdf("./indicator_species/outputs/OceanExternalTable.pdf", width = 7, height = 25)
grid.table(OceanExternalTable)
dev.off()

#CoastalInternal
CoastalInternalTable <- 
  subset(fish_table3, anat_space_combine == "Coastal_Internal") %>% 
  arrange(desc(index)) %>%
  subset(select = -c(anat_space_combine))

pdf("./indicator_species/outputs/CoastalInternalTable.pdf", width = 5, height = 8)
grid.table(CoastalInternalTable)
dev.off()

#CoastalExternal
CoastalExternalTable <- 
  subset(fish_table3, anat_space_combine == "Coastal_External") %>% 
  arrange(desc(index)) %>%
  subset(select = -c(anat_space_combine))

pdf("./indicator_species/outputs/CoastalExternalTable.pdf", width = 5, height = 5)
grid.table(CoastalExternalTable)
dev.off()



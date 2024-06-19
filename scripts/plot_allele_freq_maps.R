

#Geographic distribution of haplotypes
library(grid)
library(tidyverse)
library(ggmap)
library(gridExtra)
library(scatterpie)

clades_pallete <- c("#94B669", "#dd4124", "#0C7996")


#Prep CP data for GWAS
final_diagnostic <- read_tsv("data/wild.cp.diagnostic40.genotypes.txt" )

cp_freq <- final_diagnostic %>%
  #filter(species == "Ann") %>%
  filter(!grepl("ED",sample)) %>%
  #filter(clade == "Clade 1" | clade == "Clade 2") %>%
  group_by(species, pop, lat, long, clade) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  spread(.,clade, n,fill=0)

map_ann <- get_stamenmap(c(left = -123, bottom = 25, right = -93.5, top = 51), zoom = 5, maptype = "toner-lite")
map_petpet <- get_stamenmap(c(left = -113, bottom = 34, right = -93.5, top = 41), zoom = 6, maptype = "toner-lite")
map_petfal <- get_stamenmap(c(left = -115, bottom = 30, right = -100, top = 39), zoom = 6, maptype = "toner-lite")
map_arg <- get_stamenmap(c(left = -99, bottom = 26, right = -96, top = 30), zoom = 8, maptype = "toner-lite")

ggmap(map_ann)

pdf("figures/ann_map.v1.pdf",height=8,width=8)
ggmap(map_ann) +
  geom_scatterpie(data=cp_freq %>% filter(species == "Ann"),
                  aes(x=long, y=lat, r=0.3),
                  cols=c("Ancestral clade", "Clade 1", "Clade 2"),
                  color="black", alpha=.9,size=0.1) +
  scale_fill_manual(labels = c("Ancestral", "Clade 1", "Clade 2"),
                    values = clades_pallete)
dev.off()
      
pdf("figures/petpet_map.v1.pdf",height=5,width=8)

ggmap(map_petpet) +
  geom_scatterpie(data=cp_freq %>% filter(species == "PetPet"),
                  aes(x=long, y=lat, r=0.3),
                  cols=c("Clade 1", "Clade 2"),
                  color="black", alpha=.9,size=0.1) +
  scale_fill_manual(labels = c("Clade 1", "Clade 2"),
                    values = clades_pallete[c(2,3)])      
dev.off()

pdf("figures/petfal_map.v1.pdf",height=5,width=8)

ggmap(map_petfal) +
  geom_scatterpie(data=cp_freq %>% filter(species == "PetFal"),
                  aes(x=long, y=lat, r=0.3),
                  cols=c("Ancestral clade", "Clade 1", "Clade 2"),
                  color="black", alpha=.9,size=0.1) +
  scale_fill_manual(labels = c("Ancestral", "Clade 1", "Clade 2"),
                    values = clades_pallete) 
dev.off()

pdf("figures/arg_map.v1.pdf",height=8,width=5)

ggmap(map_arg) +
  geom_scatterpie(data=cp_freq %>% filter(species == "Arg"),
                  aes(x=long, y=lat, r=0.1),
                  cols=c("Clade 1", "Clade 2"),
                  color="black", alpha=.9,size=0.1) +
  scale_fill_manual(labels = c("Clade 1", "Clade 2"),
                    values = clades_pallete[c(2,3)])
dev.off()

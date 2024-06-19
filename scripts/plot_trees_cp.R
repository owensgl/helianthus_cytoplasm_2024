#Plotting phylogeny
#Clade_1 <- "#dd4124"
#Clade_2 <- "#0C7996"
#Ancestral <- "#94B669"

clades_pallete <- c("#94B669", "#dd4124", "#0C7996")


library(ggtree)
library(tidyverse)
library(PNWColors)
library(patchwork)
library(ape)
library(viridisLite)
library(ggtreeExtra)
library(ggthemes)

info <- read_tsv("sample_info_apr_2018.tsv") %>%
  rename(pop = population) %>%
  select(name,species, pop) %>%
  mutate(species = case_when(species == "PetCan" ~ "NivCan",
                             TRUE ~ species))
final_diagnostic <- read_tsv("data/wild.cp.diagnostic40.genotypes.txt" )

locations <- read_tsv("pop_loc_allnum.txt")
tree <- read.tree("../fasta/wild.snps.v2.cp.fa.treefile")
#tree <- drop.tip(read.tree("data/wild.snps.mt.fa.gz.treefile"), "PET0628", trim.internal = FALSE)

compressed_info <- info %>% mutate(species = case_when(species == "Gro" | species == "Gig" | species == "Div" | species == "Dec" ~ "Per",
                TRUE ~ species)) %>% as.data.frame() %>% select(-pop) %>%
  filter(name %in% tree$tip.label)

compressed_info$species <- fct_relevel(compressed_info$species, "Ann", "Arg", "Par", "PetFal", "PetPet", "Deb", "NivCan", "Ano", "Des", "Per", "Pho")


custom_colors <- calc_pal()(12)
names(custom_colors) <- c("Ann", "Arg", "Par", "PetFal", "PetPet", "Deb", "NivCan", "Ano", "Des", "Per", "Pho")



p2 <- ggtree(tree) %<+% (final_diagnostic %>% select(sample,clade)) + 
  geom_tippoint(aes(color=as.factor(clade)),size=0.6) +
  scale_color_manual(values=clades_pallete,
                    name="CP clade") +
  ggtitle("Cytoplasmic tree") +
  geom_treescale() +
  new_scale_colour()


p3 <- p2 + geom_fruit(data=compressed_info, geom=geom_point,
                mapping=aes(x=species, y=name,fill=species,color=species),
                offset = 0.01,size = 1,
                inherit.aes = FALSE) +
  scale_colour_manual(
    name="Species",
    values=custom_colors,
    guide=guide_legend(keywidth=0.5,
                       keyheight=0.5,
                       order=2,
                       override.aes=list(size=2,alpha=1))
  )

pdf("figures/CP_tree.v2test.pdf",height=9,width=9)
p3
dev.off()


mismatched <- final_diagnostic %>%
  select(sample, clade) %>%
  rename(CP_clade = clade) %>%
  inner_join(mt_clades) %>%
  filter(CP_clade != clade)

pdf("figures/CP_tree_mismatch.v2.pdf",height=9,width=9)
ggtree(tree) %<+% (mismatched) + 
  geom_tippoint(aes(color=as.factor(clade)),size=2) +
  scale_color_manual(values=clades_pallete[c(1,2,3)],
                     name="MT clade",
                     na.value = "transparent") +
  ggtitle("Chloroplast tree") +
  geom_treescale() +
  new_scale_colour()
dev.off()




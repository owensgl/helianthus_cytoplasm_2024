#Plotting phylogeny
#Clade_1 <- "#dd4124"
#Clade_2 <- "#0C7996"
#Ancestral <- "#94B669"

clades_pallete <- c("#94B669", "#dd4124", "#0C7996")
#Plotting phylogeny

library(ggtree)
library(tidyverse)
library(PNWColors)
library(patchwork)
library(ape)
library(TreeTools)
pnw_colors <-pnw_palette("Bay",5)

info <- read_tsv("sample_info_apr_2018.tsv") %>%
  rename(pop = population) %>%
  select(name,species, pop) %>%
  mutate(species = case_when(species == "PetCan" ~ "NivCan",
                             TRUE ~ species))
locations <- read_tsv("pop_loc_allnum.txt")
excluded_samples <- read_tsv("../fasta/wild.snps.mt.fa.50missingsamples.txt",col_names = "samples")
missing_info <- read_tsv("../fasta/wild.snps.mt.fa.Ncount.txt") %>%
  dplyr::rename(sample=chr) %>%
  mutate(percent_missing = Ns/bases)

tree <- read.tree("../fasta/wild.snps.v2.mt.noindels.samplepruned.fa.treefile")
#clade 1 = node 1723
#clade 2 = node 2957
clade_1_tree <- extract.clade(tree, 1723)
clade_2_tree <- extract.clade(tree, 2957)
clade_1_samples <- clade_1_tree$tip.label
clade_2_samples <- clade_2_tree$tip.label
ancestral_samples <- tree$tip.label[(!tree$tip.label %in% c(clade_1_samples, clade_2_samples))]

tmp_1 <- tibble(sample=clade_1_samples,clade="Clade 2")
tmp_2 <- tibble(sample=clade_2_samples,clade="Clade 1")
tmp_3 <- tibble(sample=ancestral_samples,clade="Ancestral clade")
mt_clades <- rbind(tmp_1, tmp_2, tmp_3)
write_tsv(mt_clades, "data/mt_phylogeny_clades.v2.txt")

compressed_info <- info %>% mutate(species = case_when(species == "Gro" | species == "Gig" | species == "Div" | species == "Dec" ~ "Per",
                                                       TRUE ~ species)) %>% as.data.frame() %>% select(-pop) %>%
  filter(name %in% tree$tip.label)

compressed_info$species <- fct_relevel(compressed_info$species, "Ann", "Arg", "Par", "PetFal", "PetPet", "Deb", "NivCan", "Ano", "Des", "Per", "Pho")


custom_colors <- calc_pal()(12)
names(custom_colors) <- c("Ann", "Arg", "Par", "PetFal", "PetPet", "Deb", "NivCan", "Ano", "Des", "Per", "Pho")

  

p2 <- ggtree(tree) %<+% (mt_clades) + 
  geom_tippoint(aes(color=as.factor(clade)),size=0.6) +
  scale_color_manual(values=clades_pallete,
                     name="MT clade") +
  ggtitle("Mitochondrial tree") +
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
pdf("figures/MT_tree.v2.pdf",height=9,width=9)
p3
dev.off()

#Mismatched CP and MT from MT phylogeny

mismatched <- final_diagnostic %>%
  select(sample, clade) %>%
  rename(CP_clade = clade) %>%
  inner_join(mt_clades) %>%
  filter(CP_clade != clade)
  
pdf("figures/MT_tree_mismatch.v2.pdf",height=9,width=9)
ggtree(tree) %<+% (mismatched) + 
  geom_tippoint(aes(color=as.factor(CP_clade)),size=2) +
  scale_color_manual(values=clades_pallete[c(2,3)],
                     name="CP clade",
                     na.value = "transparent") +
  ggtitle("Mitochondrial tree") +
  geom_treescale() +
  new_scale_colour()
dev.off()


#This script makes a plot of haplotype frequencies for each species
library(PNWColors)
library(tidyverse)
library(broom)
library(cowplot)
library(patchwork)
library(ggfortify)
library(ggrepel)

color_set_1 <- pnw_palette("Bay",6)
color_set_2 <- pnw_palette("Cascades",6)

final_diagnostic <- read_tsv("data/wild.cp.diagnostic40.genotypes.txt" )
removal_samples <- read_tsv("data/samples_to_remove.txt")
final_diagnostic <- final_diagnostic %>% anti_join(removal_samples)

compressed_info <- final_diagnostic %>% mutate(species = case_when(species == "Gro" | species == "Gig" | species == "Div" | species == "Dec" ~ "Per",
                                                                   species == "PetCan" ~ "NivCan",
                                                       TRUE ~ species)) %>% as.data.frame()

pdf("figures/cp_species_freq.v2.pdf",height=10,width=10)
compressed_info %>%
  group_by(species, clade) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n),
         total_n = sum(n)) %>%
  ggplot(.,aes(y=species,x=freq,fill=clade)) +
  geom_bar(position="stack",stat="identity") +
  scale_fill_manual(values=color_set_1[c(3,6,2)],
                    name="CP clade") +
  theme_cowplot() +
  scale_y_discrete(limits = c("Pho", "Per", "Ano","Des", "NivCan", "Deb","PetFal","PetPet","Par","Arg", "Ann")) +
  scale_x_continuous(breaks=c(0.00,0.25,0.5,0.75,1)) +
  geom_text(aes(y=species,x=1.1,label=total_n)) +
  xlab("Frequency") +
  ylab("Species")
dev.off()

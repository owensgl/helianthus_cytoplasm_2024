#This is for comparing the ABBA BABA scores to see if mismatched populations are more introgressed.

library(PNWColors)
library(tidyverse)
library(broom)
library(cowplot)
library(patchwork)
library(missMDA)
library(ggfortify)
library(ggrepel)
library(tidyr)
library(infer)

color_set_1 <- pnw_palette("Bay",6)
color_set_2 <- pnw_palette("Cascades",6)

final_diagnostic <- read_tsv("data/wild.cp.diagnostic40.genotypes.txt" )
removal_samples <- read_tsv("data/samples_to_remove.txt")
final_diagnostic <- final_diagnostic %>% anti_join(removal_samples)

compressed_info <- final_diagnostic %>% mutate(species = case_when(species == "Gro" | species == "Gig" | species == "Div" | species == "Dec" ~ "Per",
                                                                   species == "PetCan" ~ "NivCan",
                                                                   TRUE ~ species)) %>% as.data.frame()
#Samples where the cytotype does not match the species.

ann_mismatch <- compressed_info %>%
  group_by(species, pop, clade) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n),
         total_n = sum(n)) %>% 
  filter(species == "Ann") %>%
  filter(clade == "Clade 1" & n > 1) %>%
  pull(pop)

petpet_mismatch <- compressed_info %>%
  group_by(species, pop, clade) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n),
         total_n = sum(n)) %>% 
  filter(species == "PetPet") %>%
  filter(clade == "Clade 1" & n > 1) %>%
  pull(pop)

petfal_mismatch <- compressed_info %>%
  group_by(species, pop, clade) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n),
         total_n = sum(n)) %>% 
  filter(species == "PetFal") %>%
  filter(clade == "Clade 2" & n > 1) %>%
  pull(pop)

arg_mismatch <- compressed_info %>%
  group_by(species, pop, clade) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n),
         total_n = sum(n)) %>% 
  filter(species == "Arg") %>%
  filter(clade == "Clade 1" & n > 1) %>%
  pull(pop)

#Samples where the cytotype matches the species.
niv_match <- compressed_info %>%
  group_by(species, pop) %>%
  filter(species == "NivCan") %>%
  summarise(n = n()) %>%
  pull(pop)

ann_match <- compressed_info %>%
  group_by(species, pop) %>%
  summarise(n = n()) %>%
  filter(species == "Ann") %>%
  filter(!pop %in% ann_mismatch) %>%
  pull(pop)

petpet_match <- compressed_info %>%
  group_by(species, pop) %>%
  summarise(n = n()) %>%
  filter(species == "PetPet") %>%
  filter(!pop %in% petpet_mismatch) %>%
  pull(pop)

petfal_match <- compressed_info %>%
  group_by(species, pop) %>%
  summarise(n = n()) %>%
  filter(species == "PetFal") %>%
  filter(!pop %in% petfal_mismatch) %>%
  pull(pop)

arg_match <- compressed_info %>%
  group_by(species, pop) %>%
  summarise(n = n()) %>%
  filter(species == "Arg") %>%
  filter(!pop %in% arg_mismatch) %>%
  pull(pop)

petpet <- c(petpet_match,petpet_mismatch)
petfal <- c(petfal_match,petfal_mismatch)
niv <- niv_match
ann <- c(ann_match,ann_mismatch)
arg <- c(arg_match,arg_mismatch)


pop_matches <- compressed_info %>%
  filter(species == "Ann" | species == "NivCan" | species == "PetPet" | species == "PetFal" | species == "Arg") %>%
  select(species, pop) %>%
  unique() %>%
  mutate(mismatch = case_when(pop %in% c(ann_mismatch, petfal_mismatch, petpet_mismatch,arg_mismatch) ~ "mismatch",
                              TRUE ~ "match"))


#Load dstats

tmp <- read_tsv("data/w1506_outgroup.pop_5h_combined_Dmin.txt") %>%
  mutate(type = "not_flipped")

tmp_mirror <- tmp %>%
  rename(tmp_p2 = P1, tmp_p1 = P2) %>%
  rename(P1 = tmp_p1, P2 = tmp_p2) %>%
  mutate(Dstatistic = -Dstatistic, `f4-ratio` = -`f4-ratio`, type = "flip")

dstat <- rbind(tmp, tmp_mirror)  
set_1 <- dstat %>% 
  filter(P1 %in% ann) %>%
  filter(P2 %in% ann) %>%
  filter(P3 %in% petfal) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  mutate(target_species = "Ann", donor_species = "PetFal")
set_2 <- dstat %>% 
  filter(P1 %in% ann) %>%
  filter(P2 %in% ann) %>%
  filter(P3 %in% niv) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  mutate(target_species = "Ann", donor_species = "NivCan")
set_3 <- dstat %>% 
  filter(P1 %in% arg) %>%
  filter(P2 %in% arg) %>%
  filter(P3 %in% niv) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  mutate(target_species = "Arg", donor_species = "NivCan")
set_4 <- dstat %>% 
  filter(P1 %in% arg) %>%
  filter(P2 %in% arg) %>%
  filter(P3 %in% petfal) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  mutate(target_species = "Arg", donor_species = "PetFal")
set_5 <- dstat %>% 
  filter(P1 %in% petpet) %>%
  filter(P2 %in% petpet) %>%
  filter(P3 %in% petfal) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  mutate(target_species = "PetPet", donor_species = "PetFal")
set_6 <- dstat %>% 
  filter(P1 %in% petpet) %>%
  filter(P2 %in% petpet) %>%
  filter(P3 %in% niv) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  mutate(target_species = "PetPet", donor_species = "NivCan")

set_7 <- dstat %>% 
  filter(P1 %in% petfal) %>%
  filter(P2 %in% petfal) %>%
  filter(P3 %in% petpet) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  mutate(target_species = "PetFal", donor_species = "PetPet")
set_8 <- dstat %>% 
  filter(P1 %in% petfal) %>%
  filter(P2 %in% petfal) %>%
  filter(P3 %in% ann) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  mutate(target_species = "PetFal", donor_species = "Ann")
set_9 <- dstat %>% 
  filter(P1 %in% petfal) %>%
  filter(P2 %in% petfal) %>%
  filter(P3 %in% arg) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  mutate(target_species = "PetFal", donor_species = "Arg")

pdf("figures/all_populations_mean_D.v1.pdf",height=6,width=6)
rbind(set_1, set_2, set_3, set_4, set_5, set_6, set_7, set_8, set_9) %>%
  ggplot(.,aes(x=donor_species,y=mean_D,color=mismatch)) +
  geom_boxplot() +
  facet_wrap(~target_species,scales="free") +
  theme_cowplot() +
  ylab("Mean D") + xlab("Donor species") +
  scale_color_manual(values=c("#dd4124", "#0C7996"),
                     name="Cytotype",
                     labels=c("Matched", "Mismatched"))
dev.off()

##Make T tests
set_1 <- dstat %>% 
  filter(P1 %in% ann) %>%
  filter(P2 %in% ann) %>%
  filter(P3 %in% petfal) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  t_test(formula = mean_D ~ mismatch,
         alternative = "two-sided") %>%
  mutate(target_species = "Ann", donor_species = "PetFal") 
set_2 <- dstat %>% 
  filter(P1 %in% ann) %>%
  filter(P2 %in% ann) %>%
  filter(P3 %in% niv) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  t_test(formula = mean_D ~ mismatch,
         alternative = "two-sided") %>%
  mutate(target_species = "Ann", donor_species = "NivCan")
set_3 <- dstat %>% 
  filter(P1 %in% arg) %>%
  filter(P2 %in% arg) %>%
  filter(P3 %in% niv) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  t_test(formula = mean_D ~ mismatch,
         alternative = "two-sided") %>%
  mutate(target_species = "Arg", donor_species = "NivCan")
set_4 <- dstat %>% 
  filter(P1 %in% arg) %>%
  filter(P2 %in% arg) %>%
  filter(P3 %in% petfal) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  t_test(formula = mean_D ~ mismatch,
         alternative = "two-sided") %>%
  mutate(target_species = "Arg", donor_species = "PetFal") 
set_5 <- dstat %>% 
  filter(P1 %in% petpet) %>%
  filter(P2 %in% petpet) %>%
  filter(P3 %in% petfal) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  t_test(formula = mean_D ~ mismatch,
         alternative = "two-sided") %>%
  mutate(target_species = "PetPet", donor_species = "PetFal")
  
set_6 <- dstat %>% 
  filter(P1 %in% petpet) %>%
  filter(P2 %in% petpet) %>%
  filter(P3 %in% niv) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  t_test(formula = mean_D ~ mismatch,
         alternative = "two-sided") %>%
  mutate(target_species = "PetPet", donor_species = "NivCan") 

set_7 <- dstat %>% 
  filter(P1 %in% petfal) %>%
  filter(P2 %in% petfal) %>%
  filter(P3 %in% petpet) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  t_test(formula = mean_D ~ mismatch,
         alternative = "two-sided") %>%
  mutate(target_species = "PetFal", donor_species = "PetPet")
set_8 <- dstat %>% 
  filter(P1 %in% petfal) %>%
  filter(P2 %in% petfal) %>%
  filter(P3 %in% ann) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  t_test(formula = mean_D ~ mismatch,
         alternative = "two-sided") %>%
  mutate(target_species = "PetFal", donor_species = "Ann") 
set_9 <- dstat %>% 
  filter(P1 %in% petfal) %>%
  filter(P2 %in% petfal) %>%
  filter(P3 %in% arg) %>% 
  group_by(P2) %>%
  summarize(mean_D = mean(Dstatistic)) %>%
  inner_join(pop_matches %>% rename(P2 = pop) %>% select(-species)) %>%
  t_test(formula = mean_D ~ mismatch,
         alternative = "two-sided") %>%
  mutate(target_species = "PetFal", donor_species = "Arg")

rbind(set_1, set_2, set_3, set_4, set_5, set_6, set_7, set_8, set_9) %>%
  write_tsv("data/all_populations_mean_D_Ttest.txt")

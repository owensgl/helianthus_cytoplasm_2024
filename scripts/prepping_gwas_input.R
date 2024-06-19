library(tidyverse)
library(PNWColors)

#Prep CP data for GWAS
folder <- "/media/drive_5_usb/owens/working/cytoplasm"
final_diagnostic <- read_tsv("data/wild.cp.diagnostic40.genotypes.txt" )

final_diagnostic %>%
  filter(species == "Ann") %>%
  filter(!grepl("ED",sample)) %>%
  mutate(clade = case_when(clade == "clade_1" ~ 0.0,
                           clade == "clade_2" ~ 1.0,
                           TRUE ~ NA_real_)) %>%
  select(sample, clade) %>%
  mutate(sample2 = sample) %>%
  select(sample, sample2, clade) %>%
  arrange(sample) %>%
  
  write_tsv(., paste0(folder,"/gwas/annuus/cp.clade1.2.genotype.txt"),col_names = F)

final_diagnostic %>%
  filter(species == "Arg") %>%
  mutate(clade = case_when(clade == "clade_1" ~ 0.0,
                           clade == "clade_2" ~ 1.0,
                           TRUE ~ NA_real_)) %>%
  select(sample, clade) %>%
  mutate(sample2 = sample) %>%
  select(sample, sample2, clade) %>%
  arrange(sample) %>%
  write_tsv(., paste0(folder,"/gwas/argophyllus/cp.clade1.2.genotype.txt"),col_names = F)

final_diagnostic %>%
  filter(species == "PetFal") %>%
  mutate(clade = case_when(clade == "clade_1" ~ 0.0,
                           clade == "clade_2" ~ 1.0,
                           TRUE ~ NA_real_)) %>%
  select(sample, clade) %>%
  mutate(sample2 = sample) %>%
  select(sample, sample2, clade) %>%
  arrange(sample) %>%
  
  write_tsv(., paste0(folder,"/gwas/petfal/cp.clade1.2.genotype.txt"),col_names = F)

final_diagnostic %>%
  filter(species == "PetPet") %>%
  mutate(clade = case_when(clade == "clade_1" ~ 0.0,
                           clade == "clade_2" ~ 1.0,
                           TRUE ~ NA_real_)) %>%
  select(sample, clade) %>%
  mutate(sample2 = sample) %>%
  select(sample, sample2, clade) %>%
  arrange(sample) %>%
  
  write_tsv(., paste0(folder,"/gwas/petpet/cp.clade1.2.genotype.txt"),col_names = F)


####Prepping for population GWAS

final_diagnostic %>%
  filter(species == "Ann") %>%
  filter(!grepl("ED",sample)) %>%
  mutate(clade = case_when(clade == "Clade 1" ~ 0.0,
                           clade == "Clade 2" ~ 1.0,
                           TRUE ~ NA_real_)) %>%
  group_by(pop) %>%
  summarise(clade_1_freq = mean(clade,na.rm=T)) %>%
  write_tsv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/annuus/cp.clade1.2.popfreq.txt",col_names = F)


final_diagnostic %>%
  filter(species == "Arg") %>%
  filter(!grepl("ED",sample)) %>%
  mutate(clade = case_when(clade == "Clade 1" ~ 0.0,
                           clade == "Clade 2" ~ 1.0,
                           TRUE ~ NA_real_)) %>%
  group_by(pop) %>%
  summarise(clade_1_freq = mean(clade,na.rm=T)) %>%
  write_tsv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/argophyllus/cp.clade1.2.popfreq.txt",col_names = F)


final_diagnostic %>%
  filter(species == "PetPet") %>%
  filter(!grepl("ED",sample)) %>%
  mutate(clade = case_when(clade == "Clade 1" ~ 0.0,
                           clade == "Clade 2" ~ 1.0,
                           TRUE ~ NA_real_)) %>%
  group_by(pop) %>%
  summarise(clade_1_freq = mean(clade,na.rm=T)) %>%
  write_tsv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/petpet/cp.clade1.2.popfreq.txt",col_names = F)


final_diagnostic %>%
  filter(species == "PetFal") %>%
  filter(!grepl("ED",sample)) %>%
  mutate(clade = case_when(clade == "Clade 1" ~ 0.0,
                           clade == "Clade 2" ~ 1.0,
                           TRUE ~ NA_real_)) %>%
  group_by(pop) %>%
  summarise(clade_1_freq = mean(clade,na.rm=T)) %>%
  write_tsv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/petfal/cp.clade1.2.popfreq.txt",col_names = F)

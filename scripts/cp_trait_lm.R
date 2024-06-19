library(PNWColors)
library(tidyverse)
library(broom)
library(cowplot)
library(ggpubr)
library(ggpmisc)
library(patchwork)
library(nlme)
library(car)
library(ggrepel)


color_set_1 <- pnw_palette("Bay",6)
color_set_2 <- pnw_palette("Sailboat",6)

final_diagnostic <- read_tsv("data/wild.cp.diagnostic40.genotypes.txt" )
trait_categories <- read_tsv("data/all_traits_category.txt",col_names = c("trait","category"))




t.test.results <- tibble(trait=character(),clade_1 = numeric(),clade_2=numeric(),estimate=numeric(),estimate1=numeric(),estimate2=numeric(),
                         statistic=numeric(),p.value=numeric(),parameter=numeric(),conf.low=numeric(),conf.high=numeric(),method=character(),alternative=character(),
                         species=character())
lm.results <- tibble(trait=character(),term=character(),estimate=numeric(),std.error=numeric(),statistic=numeric(),
                     p.value=numeric(),species=character())
species_list <- c("annuus","argophyllus","petpet","petfal")

for (chosen_species in species_list){
  phenotype_data <- read_tsv(paste0("data/",chosen_species,".phenotypes.txt")) %>% select(-Genotype_ID) %>%
    pivot_longer(-sample, names_to = "trait",values_to = "value")
  
  tmp <- phenotype_data %>%
    inner_join(final_diagnostic) %>%
    filter(clade != "Ancestral clade") %>%
    mutate(clade = case_when(clade == "Clade 1" ~ "clade_1",
                             clade == "Clade 2" ~ "clade_2")) %>%
    group_by(trait) %>%
    mutate(mean_trait = mean(value,na.rm=T), sd_trait = sd(value,na.rm=T)) %>%
    ungroup() %>%
    mutate(value = (value-mean_trait)/sd_trait) %>%
    select(trait,value,clade) %>%
    group_by(trait, clade) %>% 
    nest() %>% 
    spread(key = clade, value = data) %>% 
    mutate(
      t_test = map2(clade_1, clade_2, ~{t.test(.x$value, .y$value) %>% tidy()}),
      clade_1 = map(clade_1, nrow),
      clade_2 = map(clade_2, nrow)
    ) %>% 
    unnest() %>%
    mutate(species = chosen_species) %>%
    ungroup()
  t.test.results <- rbind(tmp,t.test.results)
  
  tmp2 <- phenotype_data %>%
    inner_join(final_diagnostic) %>%
    filter(clade != "Ancestral clade") %>%
    mutate(clade = case_when(clade == "Clade 1" ~ "clade_1",
                             clade == "Clade 2" ~ "clade_2")) %>%
    group_by(trait) %>%
    mutate(mean_trait = mean(value,na.rm=T), sd_trait = sd(value,na.rm=T)) %>%
    ungroup() %>%
    mutate(value = (value-mean_trait)/sd_trait) %>%
    group_by(trait,pop) %>%
    mutate(total_called = n(),total_c1 = sum(clade == "clade_1")) %>% 
    filter(total_c1 < total_called & total_c1 > 0) %>%
    group_by(trait) %>%
    group_modify(~ broom::tidy(lm(value ~ pop + clade, data = .x))) %>%
    filter(term == "cladeclade_2") %>%
    ungroup()
  tmp2$species <- chosen_species
  lm.results <- rbind(tmp2, lm.results)
}




#################
##TESTING
######
niv_traits <- read_tsv("data/niv.phenotypes.txt") %>% select(-sample) %>% rename(sample = Genotype_ID) %>% pivot_longer(-sample, names_to="trait",values_to="value")
petpet_traits <- read_tsv("data/petpet.phenotypes.txt") %>% select(-sample) %>% rename(sample = Genotype_ID) %>% pivot_longer(-sample, names_to="trait",values_to="value")
petfal_traits <- read_tsv("data/petfal.phenotypes.txt") %>% select(-sample) %>% rename(sample = Genotype_ID) %>% pivot_longer(-sample, names_to="trait",values_to="value")
argophyllus_traits <- read_tsv("data/argophyllus.phenotypes.txt") %>% select(-sample) %>% rename(sample = Genotype_ID) %>% pivot_longer(-sample, names_to="trait",values_to="value")
annuus_traits <- read_tsv("data/annuus.phenotypes.txt") %>% select(-sample) %>% rename(sample = Genotype_ID) %>% pivot_longer(-sample, names_to="trait",values_to="value")

rbind(niv_traits, petpet_traits, petfal_traits, argophyllus_traits, annuus_traits) %>%
  inner_join(final_diagnostic) %>%
  filter(!is.na(value)) %>%
  group_by(trait,species) %>%
  summarize(mean=mean(value)) %>%
  group_by(trait) %>%
  summarize(n=n()) %>%
  filter(n <= 2) %>% pull(trait) -> ann_specific_traits

all_species_nlme <- rbind(niv_traits, petpet_traits, petfal_traits, argophyllus_traits, annuus_traits) %>%
  inner_join(final_diagnostic) %>%
  inner_join(trait_categories) %>%  
  filter(clade != "ancestral") %>%
  filter(!trait %in%  ann_specific_traits) %>%
  group_by(trait) %>%
  mutate(mean_trait = mean(value,na.rm=T), sd_trait = sd(value,na.rm=T)) %>%
  ungroup() %>%
  mutate(value = (value-mean_trait)/sd_trait) %>%
  ungroup() %>%
  mutate(value = (value-mean_trait)/sd_trait) %>%
  group_by(trait) %>%
  group_modify(~ broom::tidy(anova(lme(value ~ clade, random=~1|species,
                                 method="REML",
                                 na.action=na.omit,
                                  data = .x)))) %>%
  filter(column == "p-value") %>%
  select(trait,mean) %>%
  rename(pvalue=mean)

all_species_nlme %>%
  ungroup() %>%
  inner_join(trait_categories) %>%  
  arrange(pvalue) %>%
  mutate(order=row_number()) -> data
  
data %>%
  ggplot(.,aes(x=-order,y=-log10(pvalue),label=trait)) + 
  geom_point(aes(color=category),size=2) +
  theme_cowplot() + 
  geom_hline(yintercept=-log10(0.05),linetype="dotted") +
  scale_color_brewer(palette="Set1")
  

read_tsv("data/less_species_more_traits_pca.coords.txt") %>%
  pivot_longer(c(Dim.1,Dim.2,Dim.3,Dim.4),
               names_to = "trait", values_to = "value") %>%
  filter(clade != "ancestral") %>%
  group_by(trait) %>% 
  group_modify(~ broom::tidy(anova(lme(value ~ clade, random=~1|species,
                                       method="REML",
                                       na.action=na.omit,
                                       data = .x)))) %>%
  filter(column == "p-value")


###############
#TESTING
#######





#Load GEMMA results
gemma.results <- tibble(trait=character(),beta=numeric(),se=numeric(),p.value=numeric(),species=character())
capital_list <- c("Annuus","Argophyllus","Petiolaris","Petiolaris")
tag_list <- c("gwas","gwas","petpet","petfal")
directory_list <- c("annuus","argophyllus","petpet","petfal")
for (i in 1:4){
  for (chosen_trait in trait_categories$trait){
    if(file.exists(paste0("/home/owens/working/cytoplasm/gemma_out_",directory_list[i],"/",capital_list[i],
                          ".tranche90.snp.",tag_list[i],".90.bi.remappedHa412HO.beagle.",chosen_trait,".assoc.txt"))){
      tmp <- read_tsv(paste0("/home/owens/working/cytoplasm/gemma_out_",directory_list[i],"/",capital_list[i],
                             ".tranche90.snp.",tag_list[i],".90.bi.remappedHa412HO.beagle.",chosen_trait,".assoc.txt")) %>% 
        select(beta,se,p_wald) %>% mutate(species = directory_list[i], trait= chosen_trait)
      gemma.results <- rbind(gemma.results, tmp)
    }
  }
}

gemma.plot <- gemma.results %>%
  filter(!trait %in%  ann_specific_traits) %>%
  inner_join(trait_categories) %>%
  ggplot(.,aes(x=fct_reorder(trait,category),y=2*(-beta),ymin=2*(-beta-(se*1.96)),ymax=2*(-beta+(se*1.96)),color=category,shape=species,linetype=species)) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_linerange(position = position_dodge(width = 0.5)) +
  theme_cowplot() +
  scale_shape_manual(values=c(1,1,1,1)) +
  scale_linetype_manual(values=c(1,1,1,1)) +
  scale_color_manual(values=color_set_1[c(1,2,3,4,5,6)]) + 
  geom_hline(yintercept=0,linetype="dotted") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylab("Estimate") + xlab("Trait") + 
  ggtitle("Gemma per species") +
  guides(shape=FALSE,linetype=FALSE)




t.test.plot <- t.test.results %>%
  filter(!trait %in%  ann_specific_traits) %>%
  inner_join(trait_categories) %>%
  ggplot(.,aes(x=fct_reorder(trait,category),y=estimate,ymin=conf.low,ymax=conf.high,color=category,shape=species,linetype=species)) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_linerange(position = position_dodge(width = 0.5)) +
  theme_cowplot() +
  scale_shape_manual(values=c(1,1,1,1)) +
  scale_linetype_manual(values=c(1,1,1,1)) +
  scale_color_manual(values=color_set_1[c(1,2,3,4,5,6)]) + 
  geom_hline(yintercept=0,linetype="dotted") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylab("Estimate") + xlab("Trait") + ggtitle("T-test per species") +
  guides(shape=FALSE,linetype=FALSE)

all_trait_species <- tibble(species=character(),trait=character())
for (species in unique(lm.results$species)){
  for (trait in unique(lm.results$trait)){
    tmp <- tibble(species=species,trait=trait)
    all_trait_species <- rbind(all_trait_species,tmp)
  }
}

lm.results %>%
  mutate(significant = case_when(p.value < 0.05 ~ "Significant",
                                 TRUE ~ "Non-significant")) %>%
  inner_join(trait_categories) %>%
  full_join(all_trait_species) %>% 
  filter(!trait %in%  ann_specific_traits) %>%
  ungroup() %>%

  ggplot(.,aes(x=fct_reorder(trait,category),y=-estimate,ymin=-estimate-(std.error*1.96),ymax=-estimate+(std.error*1.96),color=category,shape=species,linetype=species)) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_linerange(position = position_dodge(width = 0.5)) +
  theme_cowplot() +
  scale_shape_manual(values=c(1,1,1,1)) +
  scale_linetype_manual(values=c(1,1,1,1)) +
  scale_color_manual(values=color_set_1[c(1,2,3,4,5,6)]) + 
  geom_hline(yintercept=0,linetype="dotted") +
  #facet_wrap(~category,ncol=1,scales="free") +
 # geom_vline(xintercept=seq(50)+0.5) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylab("Estimate") + xlab("Trait") +
  ggtitle("T-test estimate")

my.formula <- y ~ x

lm.results %>%
  inner_join(trait_categories) %>%
  filter(species == "petfal" | species == "annuus") %>%
  select(estimate, trait,species,category) %>%
  pivot_wider(names_from=species,values_from=estimate) %>%
  ggplot(.,aes(x=petfal,y=annuus)) + geom_point() +
  geom_smooth(method = "lm", se=FALSE, formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  theme_cowplot() +
  geom_hline(yintercept=0,linetype="dashed") +
  geom_vline(xintercept=0,linetype="dashed") +
  facet_wrap(~category)

t.test.results %>%
  inner_join(trait_categories) %>%
  filter(species == "argophyllus" | species == "annuus") %>%
  select(estimate, trait,species,category) %>%
  pivot_wider(names_from=species,values_from=estimate) %>%
  ggplot(.,aes(x=argophyllus,y=annuus)) + geom_point() +
  geom_smooth(method = "lm", se=FALSE, formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE) 

#Look at plot of some seed traits
niv_traits <- read_tsv("data/niv.phenotypes.txt") %>% select(-sample) %>% rename(sample = Genotype_ID) %>% pivot_longer(-sample, names_to="trait",values_to="value")
petpet_traits <- read_tsv("data/petpet.phenotypes.txt") %>% select(-sample) %>% rename(sample = Genotype_ID) %>% pivot_longer(-sample, names_to="trait",values_to="value")
petfal_traits <- read_tsv("data/petfal.phenotypes.txt") %>% select(-sample) %>% rename(sample = Genotype_ID) %>% pivot_longer(-sample, names_to="trait",values_to="value")
argophyllus_traits <- read_tsv("data/argophyllus.phenotypes.txt") %>% select(-sample) %>% rename(sample = Genotype_ID) %>% pivot_longer(-sample, names_to="trait",values_to="value")
annuus_traits <- read_tsv("data/annuus.phenotypes.txt") %>% select(-sample) %>% rename(sample = Genotype_ID) %>% pivot_longer(-sample, names_to="trait",values_to="value")


rbind(niv_traits, petpet_traits, petfal_traits, argophyllus_traits, annuus_traits) %>%
  inner_join(final_diagnostic) %>%
  filter(!is.na(value)) %>%
  group_by(trait,species) %>%
  summarize(mean=mean(value)) %>%
  group_by(trait) %>%
  summarize(n=n()) %>%
  filter(n <= 2) %>% pull(trait) -> ann_specific_traits

pdf("figures/cp_seed_traits.v1.pdf",height=17,width=17)
rbind(niv_traits, petpet_traits, petfal_traits, argophyllus_traits, annuus_traits) %>%
  inner_join(final_diagnostic) %>%
  inner_join(trait_categories) %>%
  filter(category== "Seed") %>%
  filter(clade != "ancestral") %>%
  ggplot(.,aes(x=species,y=value,color=clade)) + geom_point(position = position_jitterdodge()) +
  geom_boxplot(alpha=0.8) + theme_cowplot() + facet_wrap(~trait,scales="free_y") +
  scale_color_manual(values=color_set_1[c(1,6)])  

rbind(niv_traits, petpet_traits, petfal_traits, argophyllus_traits, annuus_traits) %>%
  inner_join(final_diagnostic) %>%
  inner_join(trait_categories) %>%
  filter(trait == "Distance_of_first_branching_from_ground" | 
           trait == "Primary_branches" |
           trait == "LIR" | trait == "TLN" |
           trait == "Internode_length" | 
           trait == "Peduncle_length_of_first_flower") %>%
  filter(clade != "ancestral") %>%
  ggplot(.,aes(x=species,y=value,color=clade)) + geom_point(position = position_jitterdodge()) +
  geom_boxplot(alpha=0.8) + theme_cowplot() + facet_wrap(~trait,scales="free_y") +
  scale_color_manual(values=color_set_1[c(1,6)])  
dev.off()
  


all_species_lm <- rbind(niv_traits, petpet_traits, petfal_traits, argophyllus_traits, annuus_traits) %>%
  inner_join(final_diagnostic) %>%
  inner_join(trait_categories) %>%  
  filter(clade != "ancestral") %>%
  filter(!trait %in%  ann_specific_traits) %>%
  group_by(trait) %>%
  mutate(mean_trait = mean(value,na.rm=T), sd_trait = sd(value,na.rm=T)) %>%
  ungroup() %>%
  mutate(value = (value-mean_trait)/sd_trait) %>%
  group_by(trait) %>%
  group_modify(~ broom::tidy(lm(value ~ species + clade, data = .x))) %>%
  filter(term == "cladeclade_2")

lm.plot <- all_species_lm %>% 
  inner_join(trait_categories) %>%
  ggplot(.,aes(x=fct_reorder(trait,category),y=-estimate,ymin=-estimate-(std.error*1.96),ymax=-estimate+(std.error*1.96),color=category)) + 
  geom_point(position = position_dodge(width = 1)) + 
  geom_linerange(position = position_dodge(width = 1)) +
  theme_cowplot() +
  scale_color_manual(values=color_set_1[c(1,2,3,4,5,6)]) + 
  geom_hline(yintercept=0,linetype="dotted") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylab("Estimate") + xlab("Trait") +
  ggtitle("lm(value ~ species + cp_clade)")

pdf("figures/cp.trait_associations.v1.pdf",height=15,width=15)
t.test.plot / gemma.plot / lm.plot
dev.off()

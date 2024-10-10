library(PNWColors)
library(tidyverse)
library(broom)
library(cowplot)
library(patchwork)
library(ggfortify)
library(ggrepel)
color_set_1 <- pnw_palette("Bay",6)
color_set_2 <- pnw_palette("Cascades",6)
#clades_pallete <- c("#94B669", "#dd4124", "#0C7996")
clades_pallete <- c("#dd4124", "#0C7996")


mt_calls <- read_tsv("data/mt_phylogeny_clades.v2.txt")

sample_species <- read_tsv("data/wild.cp.diagnostic40.genotypes.txt" ) %>%
  mutate(species = case_when(species == "PetCan" ~ "NivCan",
                             TRUE ~ species)) %>%
  dplyr::select(sample, species, pop)
final_diagnostic <- inner_join(mt_calls, sample_species)
trait_categories <- read_tsv("cytoplasm_WGS/data/all_traits_category.txt",col_names = c("trait","category"))
removal_samples <- read_tsv("data/samples_to_remove.txt")



species_list <- c("annuus","argophyllus","petpet","petfal","niv")

phenotype_all_data <- tibble(sample=character(), trait=character(),value=numeric())
for (chosen_species in species_list){
  phenotype_data <- read_tsv(paste0("cytoplasm_WGS/data/",chosen_species,".phenotypes.txt")) %>% select(-Genotype_ID) %>%
    pivot_longer(-sample, names_to = "trait",values_to = "value")
  phenotype_all_data <- rbind(phenotype_all_data, phenotype_data)
  
}
phenotype_all_data <- phenotype_all_data %>% anti_join(removal_samples)

phenotype_all_data %>%
  filter(!is.na(value)) %>%
  inner_join(final_diagnostic %>% select(sample,species)) %>%
  group_by(species) %>%
  mutate(total_samples=n_distinct(sample)) %>%
  group_by(trait,species,total_samples) %>%
  summarize(measured_n = sum(!is.na(value))) %>%
  ungroup() %>%
  mutate(percent_measured = measured_n/total_samples) %>%
  inner_join(trait_categories) %>%
  ggplot(.,aes(x=fct_reorder(trait,category), y=species)) +
  geom_tile(aes(fill=category,alpha=percent_measured)) +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme_cowplot()

phenotype_all_data %>%
  filter(!is.na(value)) %>%
  inner_join(final_diagnostic %>% select(sample,species)) %>%
  group_by(species) %>%
  mutate(total_samples=n_distinct(sample)) %>%
  group_by(trait,species,total_samples) %>%
  summarize(measured_n = sum(!is.na(value))) %>%
  ungroup() %>%
  mutate(percent_measured = measured_n/total_samples) %>%
  inner_join(trait_categories) %>%
  group_by(trait) %>%
  summarize(n=n()) -> species_per_trait

traits_in_all <- species_per_trait %>%filter(n >=5) 
traits_in_most <- species_per_trait %>%filter(n >=4) 

all_results <- tibble()
#LEAF
###analysis of leaf traits

phenotype_all_data %>%
  inner_join(trait_categories) %>%
  filter(trait %in% traits_in_all$trait) %>%
  filter(category == "Leaf") %>%
  select(-category) %>%
  group_by(sample) %>%
  filter(!is.na(value)) %>%
  mutate(n_traits = n()) %>% 
  ungroup() %>%
  mutate(max_traits = max(n_traits)) %>%
  filter(n_traits == max_traits) %>%
  select(-n_traits,-max_traits) %>%
  pivot_wider(names_from = trait, values_from = value)  %>%
  inner_join(final_diagnostic %>% select(sample,clade,species))%>%
  relocate(sample,species,clade)-> phenotype_wide



res.comp <- imputePCA(phenotype_wide %>% select(-sample,-species,-clade),) # iterativePCA algorithm


fit <- prcomp(phenotype_wide %>% select(-sample,-species,-clade),scale=T)
pr.var <- fit$sdev^2
pve <- pr.var / sum(pr.var)
leaf_pca_results <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  ungroup()


pc1_anova <- anova(lm(PC1 ~ species + clade,na.action=na.omit,
                      data = leaf_pca_results))

afss <- pc1_anova$"Sum Sq"
tmp_1 <- as_tibble(cbind(pc1_anova,PctExp=afss/sum(afss)*100)) %>%
  mutate(trait = "leaf", variable = "PC1")
pc2_anova <- anova(lm(PC2 ~ species + clade,na.action=na.omit,
                      data = leaf_pca_results))
afss <- pc2_anova$"Sum Sq"
tmp_2 <- as_tibble(cbind(pc2_anova,PctExp=afss/sum(afss)*100)) %>%
  mutate(trait = "leaf", variable = "PC2")

all_results <- rbind(all_results, tmp_1, tmp_2)


p1 <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  ggplot(.,aes(x=PC1,y=PC2)) + geom_point(aes(color=species,shape=clade),alpha=0.5) +
  theme_cowplot() +
  scale_color_manual(values=color_set_1[c(1,2,3,4,5)],
                     name="Species") +
  ylab(paste0("PC2 pve=",round(pve[2]*100,2),"%")) +
  xlab(paste0("PC1 pve=",round(pve[1]*100,2),"%")) +
  scale_shape_manual(name="CP clade",
                     values=c(1,2)) +
  ggtitle("Leaf traits")

summarized_pc <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  select(PC1,PC2,clade,species,sample) %>%
  pivot_longer(-c(clade,species,sample), names_to = "PC",values_to = "score") %>%
  group_by(species, clade, PC) %>%
  summarize(mean_score = mean(score), n=n(), sdev = sd(score)/sqrt(n)) %>%
  mutate(lower= mean_score- (1.96*sdev), upper=mean_score+(1.96*sdev))


p2 <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  select(PC1,PC2,clade,species,sample) %>%
  pivot_longer(-c(clade,species,sample), names_to = "PC",values_to = "score") %>%
  ggplot(.,) + 
  geom_point(aes(y=score,x=species,fill=clade),position = position_jitterdodge(jitter.width=0.2),alpha=0.05) +
  geom_hpline(data=summarized_pc,aes(x=species,y=mean_score,color=clade),size=0.1, width=0.3,position=position_dodge(width=0.7)) +
  geom_linerange(data=summarized_pc,aes(x=species, ymin=lower, ymax=upper,color=clade),
                 position=position_dodge(width=0.7),linewidth=1.2) +
  facet_wrap(~PC,scales="free_y") +
  theme_cowplot() +
  ylab("PC Score") +
  scale_color_manual(values=clades_pallete,
                     name="CP clade") +
  scale_fill_manual(values=clades_pallete,
                    name="CP clade") +
  xlab("Species") +
  ggtitle("Leaf traits") +
  scale_x_discrete(guide = guide_axis(n.dodge=2))

pdf("figures/leaf_pca.v1.pdf",height=8,width=6)
p1 / p2
dev.off()

###analysis of architecture traits

phenotype_all_data %>%
  inner_join(trait_categories) %>%
  filter(trait %in% traits_in_most$trait) %>%
  filter(category == "Architecture") %>%
  select(-category) %>%
  group_by(sample) %>%
  filter(!is.na(value)) %>%
  mutate(n_traits = n()) %>% 
  ungroup() %>%
  mutate(max_traits = max(n_traits)) %>%
  filter(n_traits == max_traits) %>%
  select(-n_traits,-max_traits) %>%
  pivot_wider(names_from = trait, values_from = value)  %>%
  inner_join(final_diagnostic %>% select(sample,clade,species))%>%
  relocate(sample,species,clade)-> phenotype_wide



res.comp <- imputePCA(phenotype_wide %>% select(-sample,-species,-clade),) # iterativePCA algorithm


fit <- prcomp(phenotype_wide %>% select(-sample,-species,-clade),scale=T)
pr.var <- fit$sdev^2
pve <- pr.var / sum(pr.var)
architecture_pca_results <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  ungroup()


pc1_anova <- anova(lm(PC1 ~ species + clade,na.action=na.omit,
                      data = architecture_pca_results))
afss <- pc1_anova$"Sum Sq"
tmp_1 <- as_tibble(cbind(pc1_anova,PctExp=afss/sum(afss)*100)) %>%
  mutate(trait = "architecture", variable = "PC1")
pc2_anova <- anova(lm(PC2 ~ species + clade,na.action=na.omit,
                      data = architecture_pca_results))
afss <- pc2_anova$"Sum Sq"
tmp_2 <- as_tibble(cbind(pc2_anova,PctExp=afss/sum(afss)*100)) %>%
  mutate(trait = "architecture", variable = "PC2")

all_results <- rbind(all_results, tmp_1, tmp_2)

p1 <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  ggplot(.,aes(x=PC1,y=PC2)) + geom_point(aes(color=species,shape=clade),alpha=0.5) +
  theme_cowplot() +
  scale_color_manual(values=color_set_1[c(1,2,3,4,5)],
                     name="Species") +
  ylab(paste0("PC2 pve=",round(pve[2]*100,2),"%")) +
  xlab(paste0("PC1 pve=",round(pve[1]*100,2),"%")) +
  scale_shape_manual(name="CP clade",
                     values=c(1,2)) +
  ggtitle("Architecture traits")

summarized_pc <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  select(PC1,PC2,clade,species,sample) %>%
  pivot_longer(-c(clade,species,sample), names_to = "PC",values_to = "score") %>%
  group_by(species, clade, PC) %>%
  summarize(mean_score = mean(score), n=n(), sdev = sd(score)/sqrt(n)) %>%
  mutate(lower= mean_score- (1.96*sdev), upper=mean_score+(1.96*sdev))


p2 <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  select(PC1,PC2,clade,species,sample) %>%
  pivot_longer(-c(clade,species,sample), names_to = "PC",values_to = "score") %>%
  ggplot(.,) + 
  geom_point(aes(y=score,x=species,fill=clade),position = position_jitterdodge(jitter.width=0.2),alpha=0.05) +
  geom_hpline(data=summarized_pc,aes(x=species,y=mean_score,color=clade),size=0.1, width=0.3,position=position_dodge(width=0.7)) +
  geom_linerange(data=summarized_pc,aes(x=species, ymin=lower, ymax=upper,color=clade),
                 position=position_dodge(width=0.7),linewidth=1.2) +
  facet_wrap(~PC,scales="free_y") +
  theme_cowplot() +
  ylab("PC Score") +
  scale_color_manual(values=clades_pallete,
                     name="CP clade") +
  scale_fill_manual(values=clades_pallete,
                    name="CP clade") +
  xlab("Species") +
  ggtitle("Architecture traits") +
  scale_x_discrete(guide = guide_axis(n.dodge=2))

pdf("figures/architecture_pca.v1.pdf",height=8,width=6)
p1 / p2
dev.off()

###analysis of Inflorescence

phenotype_all_data %>%
  inner_join(trait_categories) %>%
  filter(trait %in% traits_in_most$trait) %>%
  filter(category == "Inflorescence") %>%
  select(-category) %>%
  group_by(sample) %>%
  filter(!is.na(value)) %>%
  mutate(n_traits = n()) %>% 
  ungroup() %>%
  mutate(max_traits = max(n_traits)) %>%
  filter(n_traits == max_traits) %>%
  select(-n_traits,-max_traits) %>%
  pivot_wider(names_from = trait, values_from = value)  %>%
  inner_join(final_diagnostic %>% select(sample,clade,species))%>%
  relocate(sample,species,clade)-> phenotype_wide



res.comp <- imputePCA(phenotype_wide %>% select(-sample,-species,-clade),) # iterativePCA algorithm


fit <- prcomp(phenotype_wide %>% select(-sample,-species,-clade),scale=T)
pr.var <- fit$sdev^2
pve <- pr.var / sum(pr.var)
inflorescence_pca_results <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  ungroup()


pc1_anova <- anova(lm(PC1 ~ species + clade,na.action=na.omit,
                      data = inflorescence_pca_results))
afss <- pc1_anova$"Sum Sq"
tmp_1 <- as_tibble(cbind(pc1_anova,PctExp=afss/sum(afss)*100)) %>%
  mutate(trait = "inflorescence", variable = "PC1")
pc2_anova <- anova(lm(PC2 ~ species + clade,na.action=na.omit,
                      data = inflorescence_pca_results))
afss <- pc2_anova$"Sum Sq"
tmp_2 <- as_tibble(cbind(pc2_anova,PctExp=afss/sum(afss)*100)) %>%
  mutate(trait = "inflorescence", variable = "PC2")

all_results <- rbind(all_results, tmp_1, tmp_2)

p1 <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  ggplot(.,aes(x=PC1,y=PC2)) + geom_point(aes(color=species,shape=clade),alpha=0.5) +
  theme_cowplot() +
  scale_color_manual(values=color_set_1[c(1,2,3,4,5)],
                     name="Species") +
  ylab(paste0("PC2 pve=",round(pve[2]*100,2),"%")) +
  xlab(paste0("PC1 pve=",round(pve[1]*100,2),"%")) +
  scale_shape_manual(name="CP clade",
                     values=c(1,2)) +
  ggtitle("Inflorescence traits")

summarized_pc <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  select(PC1,PC2,clade,species,sample) %>%
  pivot_longer(-c(clade,species,sample), names_to = "PC",values_to = "score") %>%
  group_by(species, clade, PC) %>%
  summarize(mean_score = mean(score), n=n(), sdev = sd(score)/sqrt(n)) %>%
  mutate(lower= mean_score- (1.96*sdev), upper=mean_score+(1.96*sdev))


p2 <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  select(PC1,PC2,clade,species,sample) %>%
  pivot_longer(-c(clade,species,sample), names_to = "PC",values_to = "score") %>%
  ggplot(.,) + 
  geom_point(aes(y=score,x=species,fill=clade),position = position_jitterdodge(jitter.width=0.2),alpha=0.05) +
  geom_hpline(data=summarized_pc,aes(x=species,y=mean_score,color=clade),size=0.1, width=0.3,position=position_dodge(width=0.7)) +
  geom_linerange(data=summarized_pc,aes(x=species, ymin=lower, ymax=upper,color=clade),
                 position=position_dodge(width=0.7),linewidth=1.2) +
  facet_wrap(~PC,scales="free_y") +
  theme_cowplot() +
  ylab("PC Score") +
  scale_color_manual(values=clades_pallete,
                     name="CP clade") +
  scale_fill_manual(values=clades_pallete,
                    name="CP clade") +
  xlab("Species") +
  ggtitle("Inflorescence traits") +
  scale_x_discrete(guide = guide_axis(n.dodge=2))

pdf("figures/inflorescence_pca.v1.pdf",height=8,width=6)
p1 / p2
dev.off()

###analysis of seed traits

phenotype_all_data %>%
  inner_join(trait_categories) %>%
  filter(trait %in% traits_in_most$trait) %>%
  filter(category == "Seed") %>%
  select(-category) %>%
  group_by(sample) %>%
  filter(!is.na(value)) %>%
  mutate(n_traits = n()) %>% 
  ungroup() %>%
  mutate(max_traits = max(n_traits)) %>%
  filter(n_traits == max_traits) %>%
  select(-n_traits,-max_traits) %>%
  pivot_wider(names_from = trait, values_from = value)  %>%
  inner_join(final_diagnostic %>% select(sample,clade,species))%>%
  relocate(sample,species,clade)-> phenotype_wide





fit <- prcomp(phenotype_wide %>% select(-sample,-species,-clade),scale=T)
pr.var <- fit$sdev^2
pve <- pr.var / sum(pr.var)
seed_pca_results <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  ungroup()



pc1_anova <- anova(lm(PC1 ~ species + clade,na.action=na.omit,
                      data = seed_pca_results))
afss <- pc1_anova$"Sum Sq"
tmp_1 <- as_tibble(cbind(pc1_anova,PctExp=afss/sum(afss)*100)) %>%
  mutate(trait = "seed", variable = "PC1")
pc2_anova <- anova(lm(PC2 ~ species + clade,na.action=na.omit,
                      data = seed_pca_results))
afss <- pc2_anova$"Sum Sq"
tmp_2 <- as_tibble(cbind(pc2_anova,PctExp=afss/sum(afss)*100)) %>%
  mutate(trait = "seed", variable = "PC2")

all_results <- rbind(all_results, tmp_1, tmp_2)

p1 <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  ggplot(.,aes(x=PC1,y=PC2, color=species,shape=clade)) + 
  geom_point(alpha=0.5) +
  theme_cowplot() +
  scale_color_manual(values=color_set_1[c(1,2,3,4,5)],
                     name="Species") +
  ylab(paste0("PC2 pve=",round(pve[2]*100,2),"%")) +
  xlab(paste0("PC1 pve=",round(pve[1]*100,2),"%")) +
  scale_shape_manual(name="CP clade",
                     values=c(1,2)) +
  ggtitle("Seed traits")

summarized_pc <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  select(PC1,PC2,clade,species,sample) %>%
  pivot_longer(-c(clade,species,sample), names_to = "PC",values_to = "score") %>%
  group_by(species, clade, PC) %>%
  summarize(mean_score = mean(score), n=n(), sdev = sd(score)/sqrt(n)) %>%
  mutate(lower= mean_score- (1.96*sdev), upper=mean_score+(1.96*sdev))


p2 <- as_tibble(fit$x) %>%
  cbind(.,phenotype_wide$sample ) %>%
  rename(sample = `phenotype_wide$sample`) %>%
  inner_join(.,final_diagnostic) %>%
  select(PC1,PC2,clade,species,sample) %>%
  pivot_longer(-c(clade,species,sample), names_to = "PC",values_to = "score") %>%
  ggplot(.,) + 
  geom_point(aes(y=score,x=species,fill=clade),position = position_jitterdodge(jitter.width=0.2),alpha=0.05) +
  geom_hpline(data=summarized_pc,aes(x=species,y=mean_score,color=clade),size=0.1, width=0.3,position=position_dodge(width=0.7)) +
  geom_linerange(data=summarized_pc,aes(x=species, ymin=lower, ymax=upper,color=clade),
                 position=position_dodge(width=0.7),linewidth=1.2) +
  facet_wrap(~PC,scales="free_y") +
  theme_cowplot() +
  ylab("PC Score") +
  scale_color_manual(values=clades_pallete,
                     name="CP clade") +
  scale_fill_manual(values=clades_pallete,
                    name="CP clade") +
  xlab("Species") +
  ggtitle("Seed traits") +
  scale_x_discrete(guide = guide_axis(n.dodge=2))

pdf("figures/seed_pca.v1.pdf",height=8,width=6)
p1 / p2
dev.off()

pdf("figures/seed_pca_alone.v1.pdf",height=5,width=8)

p2
dev.off()
#######Plotting all anovas
all_results %>%
  filter(Df == 1) %>%
  mutate(bon.p = p.adjust(`Pr(>F)`,"bonferroni")) %>%
  write_tsv("data/all_trait_lm.v1.txt")


#Zooming in on seed traits
pdf("figures/cp_clade_seed_traits.v1.pdf",height=5,width=10)
phenotype_all_data %>%
  inner_join(trait_categories) %>%
  filter(trait %in% traits_in_most$trait) %>%
  filter(category == "Seed") %>%
  group_by(trait) %>%
  mutate(mean_trait = mean(value,na.rm=T), sd_trait = sd(value,na.rm=T)) %>%
  ungroup() %>%
  mutate(value = (value-mean_trait)/sd_trait) %>%
  inner_join(final_diagnostic %>% select(sample,clade,species))%>%
  group_by(trait) %>%
  group_modify(~ broom::tidy(lm(value ~ species + clade, data = .x))) %>%
  filter(term == "cladeClade 2") %>%
  arrange(p.value) %>% 
  ggplot(.,aes(x=fct_reorder(trait,p.value),y=-estimate,ymin=-estimate-(std.error*1.96),ymax=-estimate+(std.error*1.96))) + 
  geom_point(position = position_dodge(width = 1)) + 
  geom_linerange(position = position_dodge(width = 1)) +
  theme_cowplot() +
  scale_color_manual(values=color_set_1[c(1,2,3,4,5,6)]) + 
  geom_hline(yintercept=0,linetype="dotted") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  ylab("Estimate") + xlab("Trait") +
  ggtitle("Effect size of CP clade on seed traits")
dev.off()

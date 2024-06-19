library(PNWColors)
library(tidyverse)
library(broom)
library(cowplot)
library(MASS)


color_set_1 <- pnw_palette("Bay",6)
mt_calls <- read_tsv("data/mt_phylogeny_clades.v2.txt")

sample_species <- read_tsv("data/wild.cp.diagnostic40.genotypes.txt" ) %>%
  mutate(species = case_when(species == "PetCan" ~ "NivCan",
                             TRUE ~ species)) %>%
  dplyr::select(sample, species, pop)
final_diagnostic <- inner_join(mt_calls, sample_species)
env_data <- read_tsv("env_data.txt") %>% 
  #rename(sodium=X33) %>%
  rename(pop=`Population ID`)

variable_types <- read_tsv("environment_variables.txt") %>% rename(variable_type = type)

freq_with_env <- final_diagnostic %>%
  group_by(pop,species,clade) %>%
  filter(clade != "Ancestral clade") %>%
  summarize(count=n()) %>%
  pivot_wider(names_from = clade,values_from=count) %>%
  mutate(clade_1 = case_when(is.na(`Clade 1`) ~ as.numeric(0),
                             TRUE ~ as.numeric(`Clade 1`)),
         clade_2 = case_when(is.na(`Clade 2`) ~ as.numeric(0),
                             TRUE ~ as.numeric(`Clade 2`))) %>%
  dplyr::select(-`Clade 1`,-`Clade 2`) %>%
  mutate(clade_1_freq = clade_1/(clade_1 + clade_2)) %>%
  inner_join(env_data) %>%
  dplyr::select(-clade_1,-clade_2) %>%
  pivot_longer(-c(pop,species,clade_1_freq), names_to = "variable", values_to = "score") %>%
  group_by(variable) %>%
  mutate(mean_score = mean(score,na.rm=T),sd_score = sd(score,na.rm=T)) %>%
  ungroup() %>%
  mutate(z_score = (score - mean_score)/sd_score)


# 
# ############
# #TESTING START#
# ############
# 
freq_with_env <- final_diagnostic %>%
  group_by(pop,species,clade) %>%
  filter(clade != "Ancestral clade") %>%
  summarize(count=n()) %>%
  pivot_wider(names_from = clade,values_from=count) %>%
  mutate(clade_1 = case_when(is.na(`Clade 1`) ~ as.numeric(0),
                             TRUE ~ as.numeric(`Clade 1`)),
         clade_2 = case_when(is.na(`Clade 2`) ~ as.numeric(0),
                             TRUE ~ as.numeric(`Clade 2`))) %>%
  dplyr::select(-`Clade 1`,-`Clade 2`) %>%
  mutate(clade_1_freq = clade_1/(clade_1 + clade_2)) %>%
  inner_join(env_data) %>%
  dplyr::select(-clade_1,-clade_2) %>%
  pivot_longer(-c(pop,species,clade_1_freq), names_to = "variable", values_to = "score") %>%
  group_by(variable) %>%
  mutate(mean_score = mean(score,na.rm=T),sd_score = sd(score,na.rm=T)) %>%
  ungroup() %>%
  mutate(z_score = (score - mean_score)/sd_score)
# 
# freq_with_env %>%
#   dplyr::select(pop,species,clade_1_freq,variable,z_score) %>%
#   pivot_wider(names_from=variable,values_from=z_score) -> freq_with_env_wide
# 
# full <- lme(clade_1_freq ~Latitude, random=~1|species, data = freq_with_env_wide)   
# full_1 <- lme(clade_1_freq ~Longitude + Latitude, random=~1|species, 
#               data = freq_with_env_wide, method = "ML")   
# Anova(full_1)
# z <- stepAIC(full_1, direction="both")
# 
# freq_with_env %>%
#   write_tsv("data/cp_env_data.txt")
# 


############
#TESTING END#
############
permutation_results <- tibble(variable=character(),term=character(),
                              estimate=numeric(),std.error=numeric(),
                              statistic=numeric(),p.value=numeric())
for(i in 1:1000){
  print(i)
  tmp <- freq_with_env  %>%
    group_by(species, variable) %>%
    mutate(permuted_value = z_score[sample(row_number())]) %>%
    group_by(variable) %>%
    do(tidy(lm(clade_1_freq ~ species + permuted_value, .))) %>%
    filter(term == "permuted_value") %>%
    ungroup() 
  
  permutation_results <- rbind(permutation_results,tmp)
  
  # freq_with_env  %>%
  #   group_by(species, variable) %>%
  #   mutate(permuted_value = z_score[sample(row_number())]) %>%
  #   group_by(variable) %>%
  #   group_modify(~ broom::tidy(anova(lme(value ~ clade, random=~1|species,
  #                                        method="REML",
  #                                        na.action=na.omit,
  #                                        data = .x)))) %>%
  #   do(tidy(lm(clade_1_freq ~ species + permuted_value, .))) %>%
  #   filter(term == "permuted_value") %>%
  #   ungroup() 
  
}
permutation_results$type <- "permutation"

permutation_results %>%
  ggplot(.,aes(x=p.value)) +
  geom_histogram() +
  facet_wrap(~variable)

test_results <- freq_with_env  %>%
  group_by(variable) %>%
  do(tidy(lm(clade_1_freq ~ species + z_score, .))) %>%
  filter(term == "z_score") %>%
  ungroup() %>%
  mutate(type = "test")

test_results_no_species <- freq_with_env  %>%
  group_by(variable) %>%
  do(tidy(lm(clade_1_freq ~  z_score, .))) %>%
  filter(term == "z_score") %>%
  ungroup() %>%
  mutate(type = "test") %>%
  arrange(p.value)

rbind(permutation_results) %>%
  ggplot(.,aes(p.value))+ geom_histogram() + facet_wrap(~variable) +
  geom_vline(data=test_results, aes(xintercept=p.value),color="red")

rbind(permutation_results,test_results) %>%
  group_by(variable) %>%
  mutate(percrank=rank(p.value)/length(p.value)) %>%
  filter(type == "test") %>%
  arrange(p.value) %>% View()

 test_results %>%
   inner_join(variable_types) %>%
   write_tsv(.,"data/environment_mt.lm.results.txt")

plot_1 <- test_results %>%
  inner_join(variable_types) %>%
  ggplot(.,aes(x=fct_reorder(variable,variable_type),y=-estimate,ymin=-estimate-(std.error*1.96),ymax=-estimate+(std.error*1.96),color=variable_type)) + 
  geom_point(position = position_dodge(width = 1)) + 
  geom_linerange(position = position_dodge(width = 1)) +
  theme_cowplot() +
  scale_color_manual(values=color_set_1[c(1,2,3,4,5,6)],
                     name="Variable type") + 
  geom_hline(yintercept=0,linetype="dotted") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylab("Estimate") + xlab("Environmental variable") +
  ggtitle("Environmental associations")

plot_2 <- freq_with_env  %>%
  filter(variable == "RH") %>%
  ggplot() + 
  geom_point(aes(x=z_score,y=clade_1_freq,color=species)) + 
  geom_smooth(aes(x=z_score,y=clade_1_freq,color=species),method="lm", se=FALSE) +
  geom_smooth(aes(x=z_score,y=clade_1_freq),
              method="lm",color="black",se=FALSE,linetype="dashed") +
  theme_cowplot() +
  scale_color_manual(values=color_set_1[c(1,2,3,4,5,6)],
                     name="Species") +
  ylab("Clade 1 frequency") + 
  xlab("Normalized relative humidity")


pdf("figures/cp.env_associations.v2.pdf",height=12,width=13)
plot_1 / plot_2
dev.off()

p.adjust(sort(test_results$p.value),method="holm")
freq_with_env  %>%
  group_by(variable) %>%
  do(tidy(lm(clade_1_freq ~  z_score, .))) %>%
  filter(term == "z_score") %>%
  ungroup() %>%
  mutate(type = "test") %>% View()


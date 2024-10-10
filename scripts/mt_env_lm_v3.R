library(PNWColors)
library(tidyverse)
library(broom)
library(cowplot)
library(MASS)
library(vegan)
library(missMDA)
library(patchwork)

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
  mutate(z_score = (score - mean_score)/sd_score) %>%
  mutate(variable = case_when(variable == "NA" ~ "sodium",
                              TRUE ~ variable))

freq_with_env %>%
  dplyr::select(pop,species,clade_1_freq,variable,z_score) %>%
  pivot_wider(names_from=variable,values_from=z_score) %>%
  dplyr::select(-pop) -> freq_with_env_wide

testing_variables <- colnames(freq_with_env_wide)[3:41]

analysis_data <- freq_with_env_wide
results_df <- data.frame(
  variable = character(),
  variance_explained_pc_only = numeric(),
  variance_explained_pc_with_species = numeric(),
  variance_explained_species = numeric(),
  p_value_pc_only = numeric(),
  p_value_species = numeric(),
  p_value_pc_with_species = numeric(),
  stringsAsFactors = FALSE
)

# Run models and store results
for(i in 1:length(testing_variables)) {
  # Null model for comparison
  null_model <- glm(clade_1_freq ~ 1, family = quasibinomial(link = "logit"), 
                    data = analysis_data)
  null_deviance <- null_model$null.deviance
  
  # Model 1 (Variable only)
  formula_1 <- as.formula(paste0("clade_1_freq ~ ", testing_variables[i]))
  model_1 <- glm(formula_1, family = quasibinomial(link = "logit"), data = analysis_data)
  anova_1 <- Anova(model_1)
  
  # Calculate variance explained for PC-only model
  var_exp_pc_only <- (null_deviance - model_1$deviance) / null_deviance * 100
  
  # Model 2 (Species only)
  model_species <- glm(clade_1_freq ~ species, family = quasibinomial(link = "logit"), 
                       data = analysis_data)
  
  # Model 3 (Species + PC)
  formula_2 <- as.formula(paste0("clade_1_freq ~ species + ", testing_variables[i]))
  model_2 <- glm(formula_2, family = quasibinomial(link = "logit"), data = analysis_data)
  anova_2 <- Anova(model_2, type = 2)  # Type 2 ANOVA for independent effects
  
  # Calculate additional variance explained by PC in species+PC model
  var_exp_env_with_species <- (model_species$deviance - model_2$deviance) / null_deviance * 100
  
  # Calculate variance explained by species
  var_exp_species <- (null_deviance - model_species$deviance) / null_deviance * 100
  
  # Store all results in the data frame
  results_df <- rbind(results_df, data.frame(
    variable = testing_variables[i],
    variance_explained_env_only = var_exp_pc_only,
    variance_explained_env_with_species = var_exp_env_with_species,
    variance_explained_species = var_exp_species,
    p_value_env_only = anova_1$`Pr(>Chisq)`[1],
    p_value_species = anova_2$`Pr(>Chisq)`[1],
    p_value_env_with_species = anova_2$`Pr(>Chisq)`[2]
  ))
}

results_df <- results_df %>% inner_join(variable_types) 
results_df$p_value_env_only_bon <-   p.adjust(results_df$p_value_env_only, method="bonferroni")
results_df$p_value_env_with_species_bon <-   p.adjust(results_df$p_value_env_with_species, method="bonferroni")

results_df %>%arrange(desc(variance_explained_env_only)) %>%
  write_tsv("data/mt_env_glm.v2.txt")



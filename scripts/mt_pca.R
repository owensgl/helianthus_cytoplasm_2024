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
  dplyr::select(-pop) %>%
  filter(!is.na(SOL_SALTS)) -> freq_with_env_wide

freq_with_env_wide_noclade <- freq_with_env_wide %>%
  dplyr::select(-clade_1_freq) %>%
  dplyr::select(-species)

pca_result <- prcomp(freq_with_env_wide_noclade, scale. = FALSE)  # already scaled
factoextra::fviz_eig(pca_result, addlabels = TRUE)
plot3 <- factoextra::fviz_contrib(pca_result, choice = "var", axes = 2, top = 10, 
                                  title="Contribution of variables to PC2")

pc_scores <- as.data.frame(predict(pca_result))
names(pc_scores) <- paste0("PC", 1:ncol(pc_scores))

analysis_data <- cbind(freq_with_env_wide, pc_scores)

plot1 <- analysis_data %>%
  ggplot(.,aes(x=PC1,y=PC2, color=clade_1_freq,shape=species)) + 
  geom_point() +
  scale_color_viridis_c(name="MT Clade 1\nproportion") +
  scale_shape_manual(name="Species",values=c(15:20)) +
  theme_cowplot()

plot2 <- analysis_data %>%
  ggplot(.,aes(x=PC3,y=PC4, color=clade_1_freq,shape=species)) + 
  geom_point() +
  scale_color_viridis_c(name="MT Clade 1\nproportion") +
  scale_shape_manual(name="Species",values=c(15:20)) +
  theme_cowplot()

pdf("figures/mt_pca.plots.v1.pdf",height=5,width=6)
plot1 + plot2 +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect') 
dev.off()
for(i in 1:5) {
  formula <- as.formula(paste0("clade_1_freq ~ PC", i))
  model <- glm(formula, family = quasibinomial(link = "logit"), data = analysis_data)
  anova <- Anova(model)
  formula_2 <- as.formula(paste0("clade_1_freq ~ species + PC", i))
  model_2 <- glm(formula, family = quasibinomial(link = "logit"), data = analysis_data)
  anova_2 <- Anova(model)
}

# Create empty data frame to store all results
results_df <- data.frame(
  PC = numeric(),
  variance_explained_pc_only = numeric(),
  variance_explained_pc_with_species = numeric(),
  variance_explained_species = numeric(),
  p_value_pc_only = numeric(),
  p_value_species = numeric(),
  p_value_pc_with_species = numeric(),
  stringsAsFactors = FALSE
)

# Run models and store results
for(i in 1:5) {
  # Null model for comparison
  null_model <- glm(clade_1_freq ~ 1, family = quasibinomial(link = "logit"), 
                    data = analysis_data)
  null_deviance <- null_model$null.deviance
  
  # Model 1 (PC only)
  formula_1 <- as.formula(paste0("clade_1_freq ~ PC", i))
  model_1 <- glm(formula_1, family = quasibinomial(link = "logit"), data = analysis_data)
  anova_1 <- Anova(model_1)
  
  # Calculate variance explained for PC-only model
  var_exp_pc_only <- (null_deviance - model_1$deviance) / null_deviance * 100
  
  # Model 2 (Species only)
  model_species <- glm(clade_1_freq ~ species, family = quasibinomial(link = "logit"), 
                       data = analysis_data)
  
  # Model 3 (Species + PC)
  formula_2 <- as.formula(paste0("clade_1_freq ~ species + PC", i))
  model_2 <- glm(formula_2, family = quasibinomial(link = "logit"), data = analysis_data)
  anova_2 <- Anova(model_2, type = 2)  # Type 2 ANOVA for independent effects
  
  # Calculate additional variance explained by PC in species+PC model
  var_exp_pc_with_species <- (model_species$deviance - model_2$deviance) / null_deviance * 100
  
  # Calculate variance explained by species
  var_exp_species <- (null_deviance - model_species$deviance) / null_deviance * 100
  
  # Store all results in the data frame
  results_df <- rbind(results_df, data.frame(
    PC = i,
    variance_explained_pc_only = var_exp_pc_only,
    variance_explained_pc_with_species = var_exp_pc_with_species,
    variance_explained_species = var_exp_species,
    p_value_pc_only = anova_1$`Pr(>Chisq)`[1],
    p_value_species = anova_2$`Pr(>Chisq)`[1],
    p_value_pc_with_species = anova_2$`Pr(>Chisq)`[2]
  ))
}



# Print results
print(results_df)

write_tsv(results_df, "data/pca_mt_results.txt")

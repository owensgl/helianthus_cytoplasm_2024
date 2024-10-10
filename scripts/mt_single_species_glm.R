library(PNWColors)
library(tidyverse)
library(broom)
library(cowplot)
library(MASS)
library(car)
library(corrr)

color_set_1 <- pnw_palette("Bay",6)
mt_calls <- read_tsv("data/mt_phylogeny_clades.v2.txt")

sample_species <- read_tsv("data/wild.cp.diagnostic40.genotypes.txt" ) %>%
  mutate(species = case_when(species == "PetCan" ~ "NivCan",
                             TRUE ~ species)) %>%
  dplyr::select(sample, species, pop)
final_diagnostic <- inner_join(mt_calls, sample_species)
env_data <- read_tsv("env_data.txt") %>% 
  #rename(sodium=X33) %>%
  rename(pop=`Population ID`) %>%
  dplyr::select(-Latitude, -Longitude) %>%
  inner_join(read_tsv("location_data.txt"))

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
  mutate(z_score = (score - mean_score)/sd_score)

freq_with_env %>%
  dplyr::select(pop,species,clade_1_freq,variable,z_score) %>%
  pivot_wider(names_from=variable,values_from=z_score) %>%
  dplyr::select(-pop) %>%
  filter(!is.na(SOL_SALTS))-> freq_with_env_wide

freq_with_env_wide <- freq_with_env_wide %>%
  filter(species != "NivCan")
# Function to run GLM and extract statistics for each species and variable
run_species_specific_glm <- function(data, variables) {
  # Get unique species
  species_list <- unique(data$species)
  
  # Initialize empty list to store results
  results_list <- list()
  
  # Loop through each species
  for (sp in species_list) {
    # Subset data for this species
    species_data <- data[data$species == sp, ]
    
    # Loop through each variable
    for (var in variables) {
      # Create formula
      formula <- as.formula(paste("clade_1_freq ~", var))
      
      # Run GLM
      model <- tryCatch({
        glm(formula, family=quasibinomial, data = species_data)
      }, error = function(e) {
        return(NULL)
      })
      
      # If model ran successfully, extract statistics
      if (!is.null(model)) {
        # Get model summary
        model_summary <- summary(model)
        
        # Extract coefficient info
        coef_data <- tidy(model, conf.int = TRUE)
        
        # Get row for the variable of interest (not intercept)
        var_row <- coef_data[coef_data$term == var, ]
        
        # Store results
        results_list[[length(results_list) + 1]] <- data.frame(
          Species = sp,
          Variable = var,
          Estimate = var_row$estimate,
          CI_Lower = var_row$conf.low,
          CI_Upper = var_row$conf.high,
          P_Value = var_row$p.value
        )
      }
    }
  }
  
  # Combine all results into a single dataframe
  results_df <- do.call(rbind, results_list)
  
  return(results_df)
}

# Example usage:
# Define variables to test
variables_to_test <- colnames(freq_with_env_wide)[3:41]  # Add your variables here
variable_types <- read_tsv("environment_variables.txt") %>%
  rename(Variable = variable)

# Run analysis
results <- run_species_specific_glm(freq_with_env_wide, variables_to_test)

results_ordered <- results %>%
  inner_join(variable_types) %>%
  mutate(type = factor(type, 
                       levels = c("Temperature", "Precipitation", 
                                  "Combined", "Soil", "Location")))

pdf("figures/mt_single_species_lm.pdf",height=4,width=12)
results_ordered %>%
  mutate(Variable = factor(Variable, 
                           levels = results %>% 
                             arrange(type) %>% 
                             pull(Variable) %>% 
                             unique())) %>% 
  ggplot(aes(x = Variable, y = Estimate, color = Species)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = CI_Lower, ymax = CI_Upper),
                 position = position_dodge(width = 0.5)) +
  
  scale_y_continuous( breaks = c(-4, -2, 0, 2, 4)) +
  coord_cartesian(ylim = c(-5, 5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "white")  # White background for facet labels
  ) +
  facet_grid(~type, scales = "free_x", space = "free_x") +  # Facet by type
  labs(
    y = "Estimate (with 95% CI)",
    title = "Mitochondrial genome"
  ) +
  scale_color_manual(values=custom_colors)
dev.off()

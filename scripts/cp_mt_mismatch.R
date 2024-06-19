library(tidyverse)
library(ggbeeswarm)
cp_calls <- read_tsv("data/wild.cp.diagnostic40.genotypes.txt" )
mt_calls <- read_tsv("data/mt_phylogeny_clades.v2.txt")
pop_species <- read_tsv("sample_info_apr_2018.tsv") %>%
  select(population, species) %>%
  unique() %>%
  rename(pop=population,taxon=species)


mismatch_proportion <- cp_calls %>%
  select(sample, pop, taxon, clade) %>%
  rename(cp_clade = clade) %>% 
  inner_join(mt_calls) %>%
  mutate(matched = case_when(cp_clade == clade ~ "Matched",
                             TRUE ~ "Mismatch")) %>%
  group_by(pop, matched) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(matched == "Mismatch") %>%
  right_join(cp_calls %>% distinct(pop), by = "pop") %>%
  replace_na(list(freq = 0)) %>%
  select(pop, mismatch_freq = freq)

pdf("figures/cp_mt.mismatch.pdf",height=4,width=6)
cp_calls %>%
  select(sample, pop, taxon, clade) %>%
  rename(cp_clade = clade) %>% 
  inner_join(mt_calls) %>%
  mutate(matched = case_when(cp_clade == clade ~ "Matched",
                             TRUE ~ "Mismatch")) %>%
  group_by(pop) %>%
  mutate(total = n()) %>%
  filter(total > 8) %>%
  group_by(pop, matched) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  filter(!is.na(pop)) %>%
  inner_join(pop_species) %>%
  filter(taxon != "anomalus") %>% 
  full_join(mismatch_proportion) %>%
  filter(!is.na(taxon)) %>%
  filter(taxon != "Ano") %>%
  ggplot(aes(x = fct_reorder(pop, mismatch_freq), y = freq, fill = matched)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~taxon, scales = "free_x") +
  labs(x = "Population", y = "Frequency", fill = "MT-CP Match Status") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  scale_fill_brewer(palette = "Set1")
dev.off()

cp_calls %>%
  select(sample, pop, taxon, clade) %>%
  rename(cp_clade = clade) %>% 
  inner_join(mt_calls) %>%
  mutate(matched = case_when(cp_clade == clade ~ "Matched",
                             TRUE ~ "Mismatch")) %>%
  group_by(pop) %>%
  summarize(n=n()) %>%
  filter(n == 10)

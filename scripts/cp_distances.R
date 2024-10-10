library(tidyverse)

cp_calls <- read_tsv("data/wild.cp.diagnostic40.genotypes.txt" ) %>%
  select(sample, clade)

outgroups <-  read_tsv("data/wild.cp.diagnostic40.genotypes.txt" ) %>% 
  filter(species == "Gro" | species == "Gig" | species == "Div" | species == "Dec" | species == "P.tenuifolius") %>%
  select(sample)
samples <- read_tsv("../fasta/wild.snps.v2.mt.noindels.samplepruned.samples.txt",col_names = "sample")
samples$sample

distances <- read_table("../fasta/wild.snps.v2.mt.noindels.samplepruned.fa.mldist", skip = 1,
                        col_names = c("sample_1", samples$sample))

pdf("figures/cp_distances.v1.pdf",height=4,width=6)
distances %>%
  pivot_longer(-sample_1, names_to = "sample_2",values_to = "distance") %>%
  inner_join(cp_calls %>% rename(sample_1 =sample, cp_1 = clade)) %>%
  inner_join(cp_calls %>% rename(sample_2 =sample, cp_2 = clade)) %>%
  filter( sample_1 > sample_2) %>%
  filter(!sample_1 %in% outgroups$sample, !sample_2 %in% outgroups$sample) %>%
  mutate(comparison = case_when(cp_1 > cp_2 ~ paste(cp_2,"-", cp_1),
                                TRUE ~ paste(cp_1,"-", cp_2))) %>%
  ggplot(.,aes(x=comparison,y=distance)) +
  geom_boxplot() +
  theme_minimal() +
  scale_x_discrete(limit=c('Clade 1 - Clade 1', 'Clade 2 - Clade 2', 'Clade 1 - Clade 2',
                     'Ancestral clade - Clade 1', 'Ancestral clade - Clade 2', 'Ancestral clade - Ancestral clade'),
                   labels=c("C1-C1","C2-C2","C1-C2","A-C1","A-C2","A-A")) +
  ylab("Genetic Distance") +
  xlab("Comparison")
dev.off()

distances %>%
  pivot_longer(-sample_1, names_to = "sample_2",values_to = "distance") %>%
  inner_join(cp_calls %>% rename(sample_1 =sample, cp_1 = clade)) %>%
  inner_join(cp_calls %>% rename(sample_2 =sample, cp_2 = clade)) %>%
  filter( sample_1 > sample_2) %>%
  filter(!sample_1 %in% outgroups$sample, !sample_2 %in% outgroups$sample) %>%
  mutate(comparison = case_when(cp_1 > cp_2 ~ paste(cp_1,cp_2),
                                TRUE ~ paste(cp_2,cp_1))) %>%
  group_by(comparison) %>%
  summarize(median=median(distance),
            sd = sd(distance))
  

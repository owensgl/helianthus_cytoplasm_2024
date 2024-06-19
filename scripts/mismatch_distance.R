library(tidyverse)
library(castor)
library(geosphere)
library(TreeTools)
library(cowplot)
info <- read_tsv("/media/drive_5_usb/speedy/working/cytoplasm/sample_info_apr_2018.tsv") %>%
  rename(pop = population) %>%
  select(name,species, pop)
final_diagnostic <- read_tsv("data/wild.cp.diagnostic40.genotypes.txt" )

locations <- read_tsv("/media/drive_5_usb/speedy/working/cytoplasm/pop_loc_allnum.txt")
tree <- read.tree("data/wild.snps.cp.fa.gz.treefile")

mismatched_samples <- final_diagnostic %>%
  filter(species == "Ann" & clade != "Clade 2" |
           species == "Arg" & clade != "Clade 2" |
           species == "PetPet" & clade != "Clade 2" |
           species == "PetFal" & clade != "Clade 1" |
           species == "Deb" & clade != "Clade 2" |
           species == "PetCan" & clade != "Clade 1" | 
           species == "Ano" & clade != "Clade 2") %>%
  filter(clade != "Ancestral clade")
unmapped_samples <- final_diagnostic %>%
  filter(is.na(lat)) %>% filter(clade != "Ancestral clade") %>% pull(sample)
distance_percentile <- tibble()
for (i in 1:nrow(mismatched_samples)){
  mismatched_sample_n <- mismatched_samples$sample[i]
  print(mismatched_sample_n)
  samples_to_remove <- mismatched_samples$sample[c(which(mismatched_samples$sample != mismatched_sample_n))]
  reduced_tree <- drop.tip(tree,c(samples_to_remove, unmapped_samples))
  ordered_tree <- Preorder(reduced_tree)
  tmp <- getSisters(ordered_tree, mismatched_sample_n)
  
  sub_tree <- Subtree(ordered_tree, tmp[1])
  closest_relatives <- sub_tree$tip.label
  n_relatives <- length(closest_relatives)
  if (n_relatives > 20){next}
  sample_lat <- final_diagnostic %>% filter(sample == mismatched_sample_n) %>% pull(lat)
  sample_long <- final_diagnostic %>% filter(sample == mismatched_sample_n) %>% pull(long)
  dist_tibble <- tibble()
  #Check distance to related sample
  for (n in 1:n_relatives){
    test_sample <- closest_relatives[n]
    test_lat <- final_diagnostic %>% filter(sample == test_sample) %>% pull(lat)
    test_long <- final_diagnostic %>% filter(sample == test_sample) %>% pull(long)
    tmp_tib <- tibble(distance=distHaversine(c(sample_long,sample_lat),c(test_long,test_lat)))
    dist_tibble <- rbind(dist_tibble, tmp_tib)
  }
  average_dist <- mean(dist_tibble$distance)
  
  #Check distance to unrelated samples of the same clade
  mismatched_clade <- final_diagnostic %>% filter(sample == mismatched_sample_n) %>% pull(clade)
  other_samples <- final_diagnostic %>% filter(clade == mismatched_clade) %>%
    filter(!sample %in% mismatched_samples$sample) %>%
    filter(!sample %in% unmapped_samples)
  other_dist_tibble <- tibble()
  for (n in 1:nrow(other_samples)){
    test_sample <- other_samples$sample[n]
    test_lat <- other_samples$lat[n]
    test_long <- other_samples$long[n]
    tmp_tib <- tibble(distance=distHaversine(c(sample_long,sample_lat),c(test_long,test_lat)))
    other_dist_tibble <- rbind(other_dist_tibble, tmp_tib)
  }
  percentile_score <- sum(other_dist_tibble$distance < average_dist)/nrow(other_dist_tibble)
  tmp_tib <- tibble(sample=mismatched_sample_n, percent_score = percentile_score)
  distance_percentile <- rbind(distance_percentile, tmp_tib)
}
distance_percentile %>%
  inner_join(final_diagnostic) %>%
  ggplot(aes(percent_score,fill=species)) +
           geom_histogram()
distance_percentile %>%
  summarize(mean_dist = mean(percent_score), std=sd(percent_score))
hist(distance_percentile$percent_score)

summarized_dist <- distance_percentile %>%
  inner_join(final_diagnostic) %>%
  group_by(species) %>%
  summarize(mean_percent = mean(percent_score)*100, n=n(), sdev = sd(percent_score)*100/sqrt(n)) %>%
  mutate(lower= mean_percent- (1.96*sdev), upper=mean_percent+(1.96*sdev))
            
pdf("figures/mismatch_distance.v1.pdf",height=6,width=6)
distance_percentile %>%
  inner_join(final_diagnostic) %>%
  ggplot() +
  geom_jitter(aes(x=species,y=percent_score*100), width=0.2,alpha=0.2) +
  geom_point(data=summarized_dist,aes(x=species,y=mean_percent),size=3,color="red") +
  geom_linerange(data=summarized_dist,aes(x=species, ymin=lower, ymax=upper)) +
  theme_cowplot() +
  geom_hline(yintercept = 50,linetype="dotted") +
  ylab("Distance percentile") +
  xlab("Species")
dev.off()

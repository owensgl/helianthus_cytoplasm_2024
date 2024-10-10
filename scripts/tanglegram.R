library(tangler)
tree_mt <- read.tree("../fasta/wild.snps.v2.mt.noindels.samplepruned.fa.treefile")
tree_cp <- read.tree("../fasta/wild.snps.v2.cp.fa.treefile")

cp_calls <- read_tsv("data/wild.cp.diagnostic40.genotypes.txt" )
mt_calls <- read_tsv("data/mt_phylogeny_clades.v2.txt")


meta <- cp_calls %>%
  dplyr::select(sample, pop, taxon, clade) %>%
  rename(cp_clade = clade) %>% 
  inner_join(mt_calls) %>%
  mutate(matched = case_when(cp_clade == clade ~ "Matched",
                             TRUE ~ "Mismatch")) 



# Load tree 1 and use ggtree to annotate features
tree1 <- ggtree(tree_cp)   %<+% meta 
  #geom_tiplab() 

# Load tree 2
tree2 <- ggtree(tree_mt) %<+% meta

# Draw Tanglegram
pdf("figures/cp_mt.compared.tree.pdf",height=6,width=6)
simple.tanglegram(tree1, tree2, matched, Mismatch)
dev.off()

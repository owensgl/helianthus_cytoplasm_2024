library(tidyverse)
library(PNWColors)
library(ggrastr)
library(cowplot)
#This is for plotting GWAS scores

color_set_1 <- pnw_palette("Bay",3)
chr_lengths <- read_tsv("Ha412HO.chrlengths.txt")

species_abbreviations <- c("annuus","argophyllus","petpet","petfal")
species_capitals <- c("Annuus","Argophyllus","Petiolaris","Petiolaris")
tags <- c("env","gwas","petpet","petfalgwas")

for (x in 1:4){
  species_abbreviation <- species_abbreviations[x]
  species_capital <- species_capitals[x]
  tag <- tags[x]
  
  
  chosen_trait <- "cp12"
  
  
  
  
  gwas <- read_tsv(paste0("~/working/cytoplasm/gwas/",species_abbreviation,"/",species_capital,".tranche90.snp.",tag,
                          ".90.bi.remappedHa412HO.beagle.5maf.",chosen_trait,".ps.gz"),
                   col_names=c("id","beta","sd","pvalue")) %>%
    separate(id,c("chr","pos"),"_") %>%
    mutate(pos=as.numeric(pos)) %>%
    mutate(chr_n = chr, chr= paste0("Ha412HOChr",chr_n))
  
  
  gwas_cum <- chr_lengths %>%
    select(chr,end) %>%
    rename(chr_len = end) %>%
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(gwas, ., by=c("chr"="chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(chr, pos) %>%
    mutate( poscum=pos+tot)
  
  
  axisdf = gwas_cum %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2 )
  
  
  #Load permutations
  perm_pvalues <- tibble(perm_p = numeric())
  for (i in sprintf("%03d", 1:100)){
    tmp <- read_tsv(paste0("~/working/cytoplasm/gwas/",species_abbreviation,"/",species_capital,".tranche90.snp.",tag,
                           ".90.bi.remappedHa412HO.beagle.5maf.",chosen_trait,".",i,".pvalues.txt"),col_names=c("perm_p"))
    perm_pvalues <- rbind(tmp, perm_pvalues)
  }
  p_value_cutoff <- perm_pvalues %>%
    dplyr::summarize(quants = quantile(perm_p, probs = c(0.05))) %>% pull(quants)
  
  pdf(paste0("figures/",species_abbreviation,".gwas.",chosen_trait,".pdf"),height=6,width=15)
  
  gwas_cum %>%
    mutate(logp = abs(log10(pvalue))) %>%
    filter(logp > 2) %>%
    mutate(sig = case_when(pvalue < p_value_cutoff ~ "Sig",
                           TRUE ~ "NonSig")) %>%
    ggplot(.) +
    geom_point( aes(x=poscum, y=logp,color=as.factor(chr)), alpha=0.8, size=1.3) +
    # Show all points
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
    theme_bw() +
    theme(legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    ggtitle(paste(species_abbreviation,"-",chosen_trait)) +
    xlab("Chr") +  ylab("-Log10(pvalue)") + 
    geom_hline(yintercept=-log10(p_value_cutoff),linetype="dotted")
  
  for (i in sprintf("%02d", 1:17)){
    chr_chosen <- paste0("Ha412HOChr",i)
    print(
      gwas_cum %>%
        mutate(logp = abs(log10(pvalue))) %>%
        filter(chr == chr_chosen) %>%
        mutate(sig = case_when(pvalue < p_value_cutoff ~ "Sig",
                               TRUE ~ "NonSig")) %>%
        # filter(pos > 120000000, pos < 140000000) %>%
        ggplot(.) +
        geom_point_rast( aes(x=pos/1000000, y=logp,color=as.factor(sig)), alpha=0.8, size=1.3)  +
        geom_hline(yintercept=-log10(p_value_cutoff),linetype="dotted") +
        scale_color_manual(values=color_set_1,labels=c("Not Significant","Significant"),name="Permutation\nsignificance") +
        ylab("-Log10(pvalue)") + xlab("Mbp") +
        ggtitle(paste0(species_abbreviation, " ",chr_chosen))
    )
  }
  
  dev.off()
  
  max_space <- 500000
  clustered_hits <- gwas_cum %>%
    mutate(sig = case_when(pvalue < p_value_cutoff ~ "Sig",
                           TRUE ~ "NonSig")) %>%
    filter(sig == "Sig") %>%
    group_by(chr) %>%
    mutate(before = abs(lag(pos) - pos), after = abs(lead(pos) - pos)) %>%
    filter(before < max_space | after < max_space)
  
  current_count <- 0
  current_min_pvalue <- 1
  current_start <- 1
  current_chr <- "NA"
  outlier_windows <- tibble(chr=character(),start=numeric(),
                            end=numeric(),min_pvalue=numeric(),
                            outlier_count=numeric())
  for (i in 1:nrow(clustered_hits)){
    current_count <- current_count +1
    current_chr <- clustered_hits$chr[i]
    if (clustered_hits$pvalue[i] < current_min_pvalue){
      current_min_pvalue = clustered_hits$pvalue[i] 
    }
    if (current_count == 1){
      current_start <- clustered_hits$pos[i]
    }
    if (is.na(clustered_hits$after[i])){
      tmp <- tibble(chr=current_chr,start=current_start,
                    end=clustered_hits$pos[i],min_pvalue=current_min_pvalue,
                    outlier_count=current_count)
      outlier_windows <- rbind(outlier_windows, tmp)
      current_count <- 0
      current_min_pvalue <- 1
      next
    }
    if (clustered_hits$after[i] > max_space){
      tmp <- tibble(chr=current_chr,start=current_start,
                    end=clustered_hits$pos[i],min_pvalue=current_min_pvalue,
                    outlier_count=current_count)
      outlier_windows <- rbind(outlier_windows, tmp)
      current_count <- 0
      current_min_pvalue <- 1
      next
    }
  }
  outlier_windows$species <- species_abbreviation
  write_tsv(outlier_windows, paste0("data/",species_abbreviation,".gwas.outlier_windows.txt"))
}
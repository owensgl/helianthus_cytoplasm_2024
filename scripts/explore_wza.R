library(tidyverse)
library(cowplot)
#Plot WZA output
wza_ann <- read_csv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/annuus/annuus.kendall.wza.result.csv.gz") %>%
  separate(gene, c("chr","pos"),"_",convert = T) %>%
  mutate(species = "Ann") %>%
  mutate(log_p = -log10(Z_pVal))

wza_arg <- read_csv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/argophyllus/argophyllus.kendall.wza.result.csv.gz") %>%
  separate(gene, c("chr","pos"),"_",convert = T) %>%
  mutate(species = "Arg") %>%
  mutate(log_p = -log10(Z_pVal))

wza_petpet <- read_csv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/petpet/petpet.kendall.wza.result.csv.gz") %>%
  separate(gene, c("chr","pos"),"_",convert = T) %>%
  mutate(species = "PetFal") %>%
  mutate(log_p = -log10(Z_pVal))

wza_petfal <- read_csv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/petfal/petfal.kendall.wza.result.csv.gz") %>%
  separate(gene, c("chr","pos"),"_",convert = T) %>%
  mutate(species = "PetPet") %>%
  mutate(log_p = -log10(Z_pVal))

wza <- rbind(wza_ann, wza_arg, wza_petpet, wza_petfal)

pdf("figures/wza_all_results.v2.pdf",height=8,width=10)
for (i in c(1:17)){
  
  chr_chosen <- paste0("Ha412HOChr", sprintf("%02d",i))
  print(
  wza %>%
    filter(chr == chr_chosen) %>%
    filter(log_p > 0.5) %>%
    ggplot(.,aes(x=pos/1000000,y=log_p,color=species)) +
    geom_point() +
    facet_wrap(~species,ncol=1) +
    theme_cowplot() +
    xlab("MBp") + ylab("Log10(p-value)") +
    ggtitle(chr_chosen)
  )
}
dev.off()

chr_lengths <- read_tsv("../wild_gwas_2018/Ha412HO.chrlengths.txt")


wza_cum <- chr_lengths %>%
  select(chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(wza, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( poscum=pos+tot)


axisdf = wza_cum %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2 )

pdf("figures/wza_all_results_together.v2.pdf",height=8,width=10)
wza_cum %>%
  filter(log_p > 1) %>%
  ggplot(.) +
  geom_point( aes(x=poscum, y=log_p,color=as.factor(chr)), alpha=0.8, size=1.3) +
  # Show all points
  scale_color_manual(values = rep(c("#dd4124", "#0C7996"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
  theme_cowplot() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.minor.y = element_line(color = "lightgrey",
                                          linewidth = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgrey",
                                          linewidth = 0.5,
                                          linetype = 1),
        ) +
  ggtitle("Windowed cytotype-nuclear associations") +
  xlab("Chr") +
  facet_wrap(~species,ncol=1) +
  ylab("Log10(p-value)") + xlab("Chromosome")
dev.off()

  


####Investigate regions
Ha412HOChr17: 119770000


wza %>%
  filter(chr == "Ha412HOChr02") %>%
  filter(species == "Ann") %>%
  filter(pos > 1.14e8, pos <1.16e8 ) %>%
ggplot(.,aes(x=pos/1000000,y=log_p,color=species)) +
  geom_point() +
  facet_wrap(~species,ncol=1) +
  theme_cowplot()
tmp <- read_tsv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/annuus/Annuus.tranche90.snp.gwas.90.bi.remappedHa412HO.cp.cor.txt")

tmp %>%
  filter(chr == "Ha412HOChr17") %>%
  filter(pos > 1.15e8, pos <1.25e8 ) %>% 
  ggplot(.,aes(x=pos,y=log_p)) +
  geom_point(alpha=0.2) +
  theme_classic() +
  geom_point(data=wza %>%
  filter(chr == "Ha412HOChr17") %>%
  filter(species == "Ann") %>%
  filter(pos > 1.15e8, pos <1.25e8 ),
  aes(x=pos,y=log_p),color="red")

gene_locations <- read_tsv("data/FINAL_INTERPRO.AED_.sorted.interpro.genes.locations.txt",
                           col_names = c("name","chr","start","stop"))

gene_locations %>%
  filter(chr == "Ha412HOChr17") %>%
  mutate(start_abs_dist_1 =  abs(start - 119770000),
         start_abs_dist_2 = abs(start - 119780000),
         end_abs_dist_1 =  abs(stop - 119770000),
         end_abs_dist_2 = abs(stop - 119780000)) %>%
  pivot_longer(-c(name,chr,start,stop),names_to = "category", values_to = "abs_distance") %>%
  group_by(name,chr,start,stop) %>%
  summarize(min_dist = min(abs_distance)) %>%
  filter(min_dist < 10000)

system("bcftools view /media/drive_5_usb/owens/bin/cytoplasm_WGS/data/Annuus.tranche90.snp.gwas.90.bi.remappedHa412HO.5maf.vcf.gz -r Ha412HOChr02:114550000-114870000 -o data/tmp.vcf.gz -O z")
system("tabix -p vcf data/tmp.vcf.gz")

library(SNPRelate)

snpgdsVCF2GDS("data/tmp.vcf.gz", "data/tmp.vcf.gds", method="biallelic.only")

genofile <- snpgdsOpen("data/tmp.vcf.gds")
pca <- snpgdsPCA(genofile, num.thread=2,autosome.only = F)

tab <- data.frame(sample = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

final_diagnostic <- read_tsv("data/wild.cp.diagnostic40.genotypes.txt" ) %>%
  mutate(species = case_when(species == "PetCan" ~ "NivCan",
                             TRUE ~ species)) 

tab %>%
  inner_join(final_diagnostic) %>%
  ggplot(.,aes(x=EV1,y=EV2)) +
  geom_point(aes(color=clade))


tab %>%
  inner_join(final_diagnostic) %>%
  mutate(haplotype = case_when(EV1 < -0.15 ~ 0,
                               EV1 < -0.5 ~ 1,
                               TRUE ~ 2)) %>%
  group_by(pop) %>%
  filter(clade != "Ancestral clade") %>%
  mutate(clade_n = case_when(clade == "Clade 1" ~ 1,
                             TRUE ~ 0)) %>%
  summarize(mean_clade = mean(clade_n),
            mean_haplotype=mean(haplotype)) %>%
  ggplot(.,aes(y=mean_clade,x=mean_haplotype)) + 
  geom_point() + geom_smooth(method="lm")

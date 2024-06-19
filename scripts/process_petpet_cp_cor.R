library(tidyverse)

#Trying to do pearson's correlation for CP and SNP frequency.

cp_freq <- read_tsv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/petpet/cp.clade1.2.popfreq.txt",
                    col_names = c("pop","freq"))

cp_row <- cp_freq$freq

snp_freq <- read_tsv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/petpet/Petiolaris.tranche90.snp.petpet.90.bi.remappedHa412HO.popfreq.txt")
cor_results <- matrix(nrow=nrow(snp_freq),ncol=4)
col_num <- ncol(snp_freq)
for (i in 1:nrow(snp_freq)){
  snp_row <- snp_freq[i, c(3:col_num)] %>%  unlist(., use.names=FALSE)
  cor_value <- cor(cp_row, snp_row,method="kendall")
  p.value <- cor.test(cp_row, snp_row, method="kendall",alternative="two.sided")$p.value
  cor_results[i,1] <- snp_freq[i,1] %>% unlist(., use.names=FALSE)
  cor_results[i,2] <- snp_freq[i,2] %>% unlist(., use.names=FALSE)
  cor_results[i,3] <- cor_value 
  cor_results[i,4] <- p.value 
  if (i %% 10000 == 0){print(paste0("Proccessing row ",i, " position ",cor_results[i,1], " ",cor_results[i,2]))}
  
}

tmp <- as.tibble(cor_results) %>%
  dplyr::rename(chr = V1, pos = V2, rho = V3, pvalue= V4) %>%
  mutate(pvalue=as.numeric(pvalue), log_p = -log10(pvalue),
         pos = as.numeric(pos), rank_p = dense_rank(pvalue),
         emp_p = rank_p/n()) 

write_tsv(tmp, "/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/petpet/Petiolaris.tranche90.snp.petpet.90.bi.remappedHa412HO.cp.cor.txt")


average_freq <- matrix(nrow=nrow(snp_freq),ncol=3)
col_num <- ncol(snp_freq)
for (i in 1:nrow(snp_freq)){
  freq <- mean(snp_freq[i, c(3:col_num)] %>%  unlist(., use.names=FALSE) )
  chr <- snp_freq[i,1] %>% pull()
  pos <- snp_freq[i,2] %>% pull()
  average_freq[i,1] <- chr 
  average_freq[i,2] <- pos 
  average_freq[i,3] <- freq 
  if (i %% 10000 == 0){print(paste0("Proccessing row ",i, " position ",snp_freq[i,1], " ",snp_freq[i,2]))}
  
}

average_freq_tibble <- as_tibble(average_freq) %>%
  dplyr::rename(chr = V1, pos = V2,mean_freq=V3) %>%
  mutate(pos = as.numeric(pos))

window_name <- tmp %>%
  dplyr::select(chr,pos) %>%
  mutate(window_pos = floor(pos/10000)*10000) %>%
  mutate(window = paste0(chr,"_",window_pos)) %>%
  dplyr::select(window)



#Write WZA input
cbind(tmp,window_name ) %>%
  inner_join(average_freq_tibble) %>%
  dplyr::select(window, emp_p, mean_freq) %>%
  dplyr::rename(window_name = window) %>%
  write_csv(.,"/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/petpet/Petiolaris.tranche90.snp.petpet.90.bi.remappedHa412HO.wzainput.txt")
# 
# #Plot WZA output
# wza <- read_csv("/media/drive_5_usb/owens/bin/cytoplasm_WGS/bin/WZA/test") %>%
#   separate(gene, c("chr","pos"),"_",convert = T) 
# 
# wza %>%
#   filter(chr == "Ha412HOChr03") %>%
#   mutate(log_p = -log10(Z_pVal)) %>%
#   ggplot(.,aes(x=pos,y=log_p)) +
#   geom_point()
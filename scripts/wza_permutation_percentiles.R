library(tidyverse)
library(cowplot)


#Plot permutation quantiles for raw correlations

cor_perms <- tibble()
for (i in 1:100){
  n_perm = str_pad(i, 3, pad = "0")
  tmp <- read_csv(paste0("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/petpet/","perm.",n_perm,".wzaout.txt")) %>%
    filter(Z_pVal > 0) %>%
    mutate(log_p = -log10(Z_pVal)) %>%
    summarize(quant99 = quantile(log_p, probs = 0.99,na.rm=T),
              quant999 = quantile(log_p, probs = 0.999,na.rm=T),
              quant9999 = quantile(log_p, probs = 0.9999,na.rm=T),
              quant99999 = quantile(log_p, probs = 0.99999,na.rm=T),
              quant999999 = quantile(log_p, probs = 0.999999,na.rm=T),
              quant9999999 = quantile(log_p, probs = 0.9999999,na.rm=T),
              perm=n_perm,species="permutation",
              dataset="petpet")
  cor_perms <- rbind(cor_perms,tmp)
}
for (i in 1:100){
  n_perm = str_pad(i, 3, pad = "0")
  tmp <- read_csv(paste0("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/petfal/","perm.",n_perm,".wzaout.txt")) %>%
    filter(Z_pVal > 0) %>%
    mutate(log_p = -log10(Z_pVal)) %>%
    summarize(quant99 = quantile(log_p, probs = 0.99,na.rm=T),
              quant999 = quantile(log_p, probs = 0.999,na.rm=T),
              quant9999 = quantile(log_p, probs = 0.9999,na.rm=T),
              quant99999 = quantile(log_p, probs = 0.99999,na.rm=T),
              quant999999 = quantile(log_p, probs = 0.999999,na.rm=T),
              quant9999999 = quantile(log_p, probs = 0.9999999,na.rm=T),
              perm=n_perm,species="permutation",
              dataset="petfal")
  cor_perms <- rbind(cor_perms,tmp)
}
for (i in 1:100){
  n_perm = str_pad(i, 3, pad = "0")
  tmp <- read_csv(paste0("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/annuus/","perm.",n_perm,".wzaout.txt")) %>%
    filter(Z_pVal > 0) %>%
    mutate(log_p = -log10(Z_pVal)) %>%
    summarize(quant99 = quantile(log_p, probs = 0.99,na.rm=T),
              quant999 = quantile(log_p, probs = 0.999,na.rm=T),
              quant9999 = quantile(log_p, probs = 0.9999,na.rm=T),
              quant99999 = quantile(log_p, probs = 0.99999,na.rm=T),
              quant999999 = quantile(log_p, probs = 0.999999,na.rm=T),
              quant9999999 = quantile(log_p, probs = 0.9999999,na.rm=T),
              perm=n_perm,species="permutation",
              dataset="ann")
  cor_perms <- rbind(cor_perms,tmp)
}
for (i in 1:100){
  n_perm = str_pad(i, 3, pad = "0")
  tmp <- read_csv(paste0("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/argophyllus/","perm.",n_perm,".wzaout.txt")) %>%
    filter(Z_pVal > 0) %>%
    mutate(log_p = -log10(Z_pVal)) %>%
    summarize(quant99 = quantile(log_p, probs = 0.99,na.rm=T),
              quant999 = quantile(log_p, probs = 0.999,na.rm=T),
              quant9999 = quantile(log_p, probs = 0.9999,na.rm=T),
              quant99999 = quantile(log_p, probs = 0.99999,na.rm=T),
              quant999999 = quantile(log_p, probs = 0.999999,na.rm=T),
              quant9999999 = quantile(log_p, probs = 0.9999999,na.rm=T),
              perm=n_perm,species="permutation",
              dataset="arg")
  cor_perms <- rbind(cor_perms,tmp)
}
#Grab empirical scores
cor_empirical <- tibble()

tmp <- read_csv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/petpet/petpet.kendall.wza.result.csv.gz") %>%
  filter(Z_pVal > 0) %>%
  mutate(log_p = -log10(Z_pVal)) %>%
  summarize(quant99 = quantile(log_p, probs = 0.99,na.rm=T),
            quant999 = quantile(log_p, probs = 0.999,na.rm=T),
            quant9999 = quantile(log_p, probs = 0.9999,na.rm=T),
            quant99999 = quantile(log_p, probs = 0.99999,na.rm=T),
            quant999999 = quantile(log_p, probs = 0.999999,na.rm=T),
            quant9999999 = quantile(log_p, probs = 0.9999999,na.rm=T),
            perm="real",species="petpet",
            dataset="petpet")
cor_empirical <- rbind(cor_empirical,tmp)

tmp <- read_csv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/petfal/petfal.kendall.wza.result.csv.gz") %>%
  filter(Z_pVal > 0) %>%
  mutate(log_p = -log10(Z_pVal)) %>%
  summarize(quant99 = quantile(log_p, probs = 0.99,na.rm=T),
            quant999 = quantile(log_p, probs = 0.999,na.rm=T),
            quant9999 = quantile(log_p, probs = 0.9999,na.rm=T),
            quant99999 = quantile(log_p, probs = 0.99999,na.rm=T),
            quant999999 = quantile(log_p, probs = 0.999999,na.rm=T),
            quant9999999 = quantile(log_p, probs = 0.9999999,na.rm=T),
            perm="real",species="petfal",
            dataset="petfal")
cor_empirical <- rbind(cor_empirical,tmp)

tmp <- read_csv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/annuus/annuus.kendall.wza.result.csv.gz") %>%
  filter(Z_pVal > 0) %>%
  mutate(log_p = -log10(Z_pVal)) %>%
  summarize(quant99 = quantile(log_p, probs = 0.99,na.rm=T),
            quant999 = quantile(log_p, probs = 0.999,na.rm=T),
            quant9999 = quantile(log_p, probs = 0.9999,na.rm=T),
            quant99999 = quantile(log_p, probs = 0.99999,na.rm=T),
            quant999999 = quantile(log_p, probs = 0.999999,na.rm=T),
            quant9999999 = quantile(log_p, probs = 0.9999999,na.rm=T),
            perm="real",species="ann",
            dataset="ann")
cor_empirical <- rbind(cor_empirical,tmp)

tmp <- read_csv("/media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/argophyllus/argophyllus.kendall.wza.result.csv.gz") %>%
  filter(Z_pVal > 0) %>%
  mutate(log_p = -log10(Z_pVal)) %>%
  summarize(quant99 = quantile(log_p, probs = 0.99,na.rm=T),
            quant999 = quantile(log_p, probs = 0.999,na.rm=T),
            quant9999 = quantile(log_p, probs = 0.9999,na.rm=T),
            quant99999 = quantile(log_p, probs = 0.99999,na.rm=T),
            quant999999 = quantile(log_p, probs = 0.999999,na.rm=T),
            quant9999999 = quantile(log_p, probs = 0.9999999,na.rm=T),
            perm="real",species="arg",
            dataset="arg")
cor_empirical <- rbind(cor_empirical,tmp)



pdf("figures/wza_permutations.v1.pdf",height=6,width=8)
cor_perms %>%
  pivot_longer(-c(perm,species,dataset), names_to = "quantile",values_to = "log_p") %>%
  group_by(quantile, dataset) %>%
  summarise(percentile = scales::percent(c(0.05, 0.5, 0.95)),
            score = quantile(log_p, c(0.05,0.5, 0.95))) %>%
  pivot_wider(names_from = percentile, values_from = score) %>%
  ggplot(.,aes(x=quantile,y=`50%`)) +
  geom_point() +
  geom_linerange(aes(x=quantile,ymin=`5%`,ymax=`95%`)) +
  facet_wrap(~dataset) +
  geom_line(data=cor_empirical %>%
              pivot_longer(-c(perm,species, dataset), names_to = "quantile",values_to = "log_p"),
            aes(y=log_p,x=quantile,group=perm),color="#0C7996") +
  theme_cowplot() +
  scale_x_discrete(breaks=c("quant99","quant999","quant9999","quant99999","quant999999","quant9999999"),
                   labels=c("99%", "99.9%", "99.99%",  "99.999%", "99.9999%", "99.99999%"),
                   guide = guide_axis(n.dodge=2)) +
  ylab("Log10(p-value)") + xlab("Percentile") 
dev.off()




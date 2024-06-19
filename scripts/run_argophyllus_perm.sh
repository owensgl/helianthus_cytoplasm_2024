
perm_n=$1
#Rscript cor_permutation_arg.R $perm_n
cd /media/drive_5_usb/speedy/working/cytoplasm/pop_gwas/argophyllus

cut -f 7 Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO.cp.cor.$perm_n.txt > tmp.$perm_n.txt

cut -f 1,3 -d "," Argophyllus.tranche90.snp.gwas.90.bi.remappedHa412HO.wzainput.txt | paste -d , tmp.$perm_n.txt - > $perm_n.txt

python /media/drive_5_usb/owens/bin/cytoplasm_WGS/bin/WZA/general_WZA_script.py --correlations $perm_n.txt --summary_stat emp_p --window window_name --MAF mean_freq --output perm.$perm_n.wzaout.txt --sep ,

rm tmp.$perm_n.txt
rm $perm_n.txt


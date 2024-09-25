awk 'NR==FNR{a[$1"\t"$2]=$0;next} ($1"\t"$2) in a{print $3"\t"a[($1"\t"$2)]; next}' sn_SLC8A1_SNPs.txt dbSNP_38.vcf > sn_SLC8A1_rsid.txt
awk '{print $1}' sn_SLC8A1_rsid.txt > sn_SLC8A1_SNPs.txt

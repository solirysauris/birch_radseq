
```
plink2 --vcf newly_merged.vcf.gz --make-bed --allow-extra-chr --out newly_merged
gemma -bfile newly_merged -gk 1 -miss 1 -maf 0 -hwe 0 -r2 0 -o newly_merged
gemma -bfile newly_merged  -lmm 4 -n 1 -miss 1 -maf 0 -hwe 0 -r2 0 -o new_renamed -k output/newly_merged.cXX.txt
```

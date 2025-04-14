Used code: 

Making bed: 

```
plink2 --vcf file.vcf.gz --make-bed --allow-extra-chr --out file
```
(not to forget to change phenotype (.fam) file according to genotype info)

Computing kinship & running gwas: 

```
gemma -bfile file -gk 1 -miss 1 -maf 0 -hwe 0 -r2 0 -o file
gemma -bfile file  -lmm 4 -n 1 -miss 1 -maf 0 -hwe 0 -r2 0 -o file -k output/file.cXX.txt
```

Parameters:

1. -gk 1 -miss 1 -maf 0 -hwe 0 -r2 0 for not applying additional filtration (can be removed, but gemma keeps like 20% SNPs) 
2. -lmm 4 for choosing model with kinship

Datasets: 

1. New data with lesser depth filtration, 0.8 max missing, 0.1 MAF
2. Combined data from new dataset and data from old run (RADSeq for 192 trees from two crosses) 

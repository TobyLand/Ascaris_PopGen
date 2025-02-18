# Analysis 

## Table of contents
0. [Project setup](#setup)
1. [Population Structure](#structure)
2. [Selection](#selection)

## 00 - Project setup <a name="setup"></a>
### Setup a working environment for the analysis.
```
WORKING_DIR=${PWD}
cd ${WORKING_DIR}/06_ANALYSIS

# make working directories
mkdir 00_SCRIPTS 01_STRUCTURE 02_SELECTION 03_ASSOCIATION
```
___
## 01 - Population structure <a name="structure"></a>
### Create plink files
```
cd 01_STRUCTURE

# Convert vcf to plink .bed format, select autosomes. 
plink2 --vcf ${WORKING_DIR}/06_ANALYSIS/FINAL/As_POPGEN.vcf --chr CM024168.1, CM024169.1, CM024170.1, CM024171.1, CM024172.1, CM024173.1, CM024174.1, CM024175.1, CM024176.1, CM024177.1, CM024178.1, CM024179.1, CM024180.1, CM024181.1, CM024182.1, CM024183.1, CM024184.1, CM024185.1, CM024186.1 --make-bed --allow-extra-chr --set-all-var-ids @_# --out autosomes_unfiltered

# Remove variants in strong linkage disequilibrium
plink2 --bfile autosomes_unfiltered --allow-extra-chr --set-all-var-ids @_# --indep-pairwise 50 10 0.15
plink2 --bfile autosomes_unfiltered --allow-extra-chr --set-all-var-ids @_# --extract autosomes_unfiltered.prune.in --out prunedData --make-bed
```
### Principal component analysis
```
plink2 --bfile prunedData --allow-extra-chr --set-all-var-ids @_# --pca
```
### Neighbour-joining phylogeny
```
# Produce a distance matrix
plink2 --bfile prunedData --allow-extra-chr --set-all-var-ids @_# --distance square 1-ibs
paste <( cut -f2 prunedData_tree.mdist.id) prunedData_tree.mdist | cat <(cut -f2 prunedData_tree.mdist.id | tr '\n' '\t' | sed -e '1s/^/\t/g') - > autosomes.mdist
```
```
### Nucleotide diversity, F<sub>ST</sub> and d<sub>XY</sub>
```
# Create list for each level of population
cut -f1,2 ${WORKING_DIR}/00_METADATA/suppdata.txt | grep ALB > host.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/suppdata.txt | grep ALB > household.list
cut -f1,3,6 ${WORKING_DIR}/00_METADATA/suppdata.txt | grep -v Treat | grep ALB > treatment.list

# Subset the allsites VCF for each chromosome
parallel --dry-run "vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FINAL/As_POPGEN.allsites.vcf --chr {} --recode-INFO-all --recode --out As_POPGEN.allsites.{}.vcf" ::: CM024168.1, CM024169.1, CM024170.1, CM024171.1, CM024172.1, CM024173.1, CM024174.1, CM024175.1, CM024176.1, CM024177.1, CM024178.1, CM024179.1, CM024180.1, CM024181.1, CM024182.1, CM024183.1, CM024184.1, CM024185.1, CM024186.1

# Run PIXY (this is an example which will write command for calculate the 3 statistics for each population and chromosome)
parallel --dry-run "pixy --stats fst dxy pi 
--vcf ${WORKING_DIR}/06_ANALYSIS/FINAL/As_POPGEN.allsites.{1}.vcf 
--zarr_path zarr/ 
--window_size {2}
#--reuse_zarr yes # Add after the first run has been completed
--populations {3}.list # Repeat for each of the 3 other lists
--variant_filter_expression 'DP>=10,GQ>=20,RGQ>=20' 
--invariant_filter_expression 'DP>=10,RGQ>=20' 
--outfile_prefix output/pixy.{1}.{2}.{3}" ::: CM024168.1, CM024169.1, CM024170.1, CM024171.1, CM024172.1, CM024173.1, CM024174.1, CM024175.1, CM024176.1, CM024177.1, CM024178.1, CM024179.1, CM024180.1, CM024181.1, CM024182.1, CM024183.1, CM024184.1, CM024185.1, CM024186.1::: 5000 25000 ::: host household treat

# Create headers for each statistic (doesn't matter which file you use, as long as it's specific to each statistic)
head -1 pixy.As_1.5000.school_fst.txt > fst.header
head -1 pixy.As_1.5000.school_dxy.txt > dxy.header
head -1 pixy.As_1.5000.school_pi.txt > pi.header

# Combine files for each chromosome
cat pixy.As_*.5000.household_pi.txt | grep -v pop | cat pi.header - | grep -v  > all.pi.pixy.household.txt
cat pixy.As_*.5000.host_pi.txt | grep -v pop | cat pi.header - | grep -v ZW > pi.host.txt 
cat pixy.As_*.5000.treatment_pi.txt | grep -v pop | cat pi.header - | grep -v > all.pi.treat.fix.txt
cat pixy.As_*.5000.household_fst.txt | grep -v pop | cat fst.header - | grep -v > autosomes.fst.5kb.household.txt
cat pixy.As_*.5000.household_dxy.txt | grep -v pop | cat dxy.header - | grep -v  > autosomes.dxy.5kb.houshold.txt 
cat pixy.As_*.5000.treatment_fst.txt | grep -v pop | cat fst.header - | grep -v  > fst.treatment.txt 
cat pixy.As_*.5000.treatment_dxy.txt | grep -v pop | cat dxy.header - | grep -v > autosomes.dxy.5kb.treatment.txt 
cat pixy.As_*.5000.host_fst.txt | grep -v pop | cat fst.header - | grep -v > autosomes.fst.5kb.host.txt
cat pixy.As_*.2000.treatment_fst.txt | grep -v pop | cat fst.header - | grep -v  > fst.treatment.2kb.txt 
cat pixy.As_*.5000.host_fst.txt | grep -v pop | cat fst.header | grep -v > fst.host.txt 
```
### Recombination
```
# Create input files for each treatment population
cut -f1,3 ${WORKING_DIR}/00_METADATA/suppdata.txt | grep ALB | grep "KorkeDoge" > korkedoge.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/suppdata.txt | grep ALB | grep -v "KorkeDoge" | shuf | head -17 > Compliant.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep ALB | grep -v "KorkeDoge" | grep -v -f Compliant.list | shuf | head -17 > SemiCompliant.list
cut -f1,3 ${WORKING_DIR}/00_METADATA/supplementary_data_10.txt | grep ALB | grep -v "KorkeDoge" | grep -v -f <(cat Compliant.list SemiCompliant.list) | shuf | head -17 > NonCompliant.list

parallel --dry-run "vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FINAL/As_POPGEN.vcf --geno-r2 --out {.} --keep {} --min-r2 0.1 --ld-window-bp 50000 --maf 0.05" ::: korkedoge.list Compliant.list SemiCompliant.list NonCompliant.list 

# Calculate median r2 for each distance
parallel "awk '{print $1,$3-$2,$5}' {}.geno.ld | grep -v CHR | tr ' ' '\t' | sort -k1,1 -k2,2 | datamash -g 1,2 median 3 > {}.median.txt" ::: korkedoge Compliant SemiCompliant NonCompliant

```
___
## 02 - Selection <a name="selection"></a>
### Create input VCFs
```
cd ${WORKING_DIR}/06_ANALYSIS/02_SELECTION
# Select only biallelic variants
vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FINAL/As_POPGEN.vcf --recode --recode-INFO-all --out As_POPGEN.biallelic.vcf --min-alleles 2 --max-alleles 2

# Split VCF into per-chromosome VCFs
parallel --dry-run "vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FINAL/As_POPGEN.vcf --recode --recode-INFO-all --out As_POPGEN.biallelic.{}.vcf --chr {}" ::: CM024168.1, CM024169.1, CM024170.1, CM024171.1, CM024172.1, CM024173.1, CM024174.1, CM024175.1, CM024176.1, CM024177.1, CM024178.1, CM024179.1, CM024180.1, CM024181.1, CM024182.1, CM024183.1, CM024184.1, CM024185.1, CM024186.1, CM024187.1, CM024188.1, CM024189.1, CM024190.1, CM024191.1 

# Example
vcftools --vcf As_POPGEN.biallelic.vcf --recode --recode-INFO-all --out As_POPGEN.biallelic.CM024168.1.vcf --chr CM024168.1
```
### Create genetic maps
```
# Produce a list of variant sites
bcftools query -f '%CHROM\t%POS\n' As_POPGEN.biallelic.vcf > As_POPGEN.biallelic.txt

# Produce an approximate genetic map for each chromosome (using a per-chromosome recombination rate estimates). First three chromosomes are described below with following chromosomes continuing concordantly  
awk '{$3=((($2-3000)*3.701)/1000000)}{print $1,".",$3,$2}' <(grep "CM024168.1" As_POPGEN.biallelic.txt) > As_1.gmap
awk '{$3=((($2-7397)*4.975)/1000000)}{print $1,".",$3,$2}' <(grep "CM024169.1" As_POPGEN.biallelic.txt) > As_2.gmap
awk '{$3=((($2-19901)*3.236)/1000000)}{print $1,".",$3,$2}' <(grep "CM024170.1" As_POPGEN.biallelic.txt) > As_3.gmap
```
### Statistically phase variants using BEAGLE
```
# Phase variants for each vcf file using BEAGLE
parallel --dry-run "java -jar beagle.28Sep18.793.jar gt=As_POPGEN.biallelic.{}.vcf out=As_POPGEN.biallelic.{}.beagle map={}.gmap nthreads=4 iterations=100 burnin=10 ne=60000" ::: CM024168.1, CM024169.1, CM024170.1, CM024171.1, CM024172.1, CM024173.1, CM024174.1, CM024175.1, CM024176.1, CM024177.1, CM024178.1, CM024179.1, CM024180.1, CM024181.1, CM024182.1, CM024183.1, CM024184.1, CM024185.1, CM024186.1, CM024187.1, CM024188.1, CM024189.1, CM024190.1, CM024191.1
```
### Produce subsets for each population
```
cat ${WORKING_DIR}/00_METADATA/suppdata.txt | grep "gz" | grep 'NonCompliant' | cut -f4 > NonCompliant.list
cat ${WORKING_DIR}/00_METADATA/supodata.txt | grep "gz" | grep 'SemiCompliant' | cut -f4 > SemiCompliant.list
cat ${WORKING_DIR}/00_METADATA/supodata.txt | grep "gz" | grep 'FullyCompliant' | cut -f4 > FullyCompliant.list


#EDIT KEEP MISSING VARIANTS
parallel --dry-run "vcftools --vcf As_POPGEN.biallelic.{}.beagle.vcf.gz --recode --recode-INFO-all --out As_POPGEN.biallelic.{1}.NonCompliant.vcf --chr {1} --keep {2}" ::: CM024168.1, CM024169.1, CM024170.1, CM024171.1, CM024172.1, CM024173.1, CM024174.1, CM024175.1, CM024176.1, CM024177.1, CM024178.1, CM024179.1, CM024180.1, CM024181.1, CM024182.1, CM024183.1, CM024184.1, CM024185.1, CM024186.1, CM024187.1, CM024188.1, CM024189.1, CM024190.1, CM024191.1:::: NonCompliant.list

parallel --dry-run "vcftools --vcf As_POPGEN.biallelic.{}.beagle.vcf.gz --recode --recode-INFO-all --out As_POPGEN.biallelic.{1}.SemiCompliant.vcf --chr {1} --keep {2}" ::: CM024168.1, CM024169.1, CM024170.1, CM024171.1, CM024172.1, CM024173.1, CM024174.1, CM024175.1, CM024176.1, CM024177.1, CM024178.1, CM024179.1, CM024180.1, CM024181.1, CM024182.1, CM024183.1, CM024184.1, CM024185.1, CM024186.1, CM024187.1, CM024188.1, CM024189.1, CM024190.1, CM024191.1:::: SemiCompliant.list

parallel --dry-run "vcftools --vcf As_POPGEN.biallelic.{}.beagle.vcf.gz --recode --recode-INFO-all --out As_POPGEN.biallelic.{1}.FullyCompliant.vcf --chr {1} --keep {2}" ::: CM024168.1, CM024169.1, CM024170.1, CM024171.1, CM024172.1, CM024173.1, CM024174.1, CM024175.1, CM024176.1, CM024177.1, CM024178.1, CM024179.1, CM024180.1, CM024181.1, CM024182.1, CM024183.1, CM024184.1, CM024185.1, CM024186.1, CM024187.1, CM024188.1, CM024189.1, CM024190.1, CM024191.1:::: FullyCompliant.list

```
### Calculate IHS for each district (per chromosome)
```
# Calculate IHS
parallel --dry-run "selscan --ihs --vcf As_POPGEN.biallelic.{1}.{2}.vcf --map {1}.gmap --threads 2 --out {1}.{2}.{3}" ::: CM024168.1, CM024169.1, CM024170.1, CM024171.1, CM024172.1, CM024173.1, CM024174.1, CM024175.1, CM024176.1, CM024177.1, CM024178.1, CM024179.1, CM024180.1, CM024181.1, CM024182.1, CM024183.1, CM024184.1, CM024185.1, CM024186.1, CM024187.1, CM024188.1, CM024189.1, CM024190.1, CM024191.1::: Fully Compliant NonCompliant SemiCompliant

# Normalize
norm --ihs --files CM024168.1.FullyCompliant.ihs.out CM024168.1.FullyCompliant.ihs.out CM024169.1.FullyCompliant.ihs.out CM024170.1.FullyCompliant.ihs.out CM024171.1.FullyCompliant.ihs.out CM024172.1.FullyCompliant.ihs.out CM024173.1.FullyCompliant.ihs.out CM024174.1.FullyCompliant.ihs.out CM024175.1.FullyCompliant.ihs.out CM024176.1.FullyCompliant.ihs.out CM024177.1.FullyCompliant.ihs.out CM024178.1.FullyCompliant.ihs.out CM024179.1.FullyCompliant.ihs.out CM024180.1.FullyCompliant.ihs.out CM024181.1.FullyCompliant.ihs.out CM024182.1.FullyCompliant.ihs.out CM024183.1.FullyCompliant.ihs.out CM024184.1.FullyCompliant.ihs.out CM024185.1.FullyCompliant.ihs.out CM024186.1.FullyCompliant.ihs.out CM024188.1.FullyCompliant.ihs.out
CM024189.1.FullyCompliant.ihs.out CM024189.1.FullyCompliant.ihs.out
CM024190.1.FullyCompliant.ihs.out CM024191.1.FullyCompliant.ihs.out

norm --ihs --files CM024168.1.FullyCompliant.ihs.out CM024168.1.NonCompliant.ihs.out CM024169.1.NonCompliant.ihs.out CM024170.1.NonCompliant.ihs.out CM024171.1.NonCompliant.ihs.out CM024172.1.NonCompliant.ihs.out CM024173.1.NonCompliant.ihs.out CM024174.1.NonCompliant.ihs.out CM024175.1.NonCompliant.ihs.out CM024176.1.NonCompliant.ihs.out CM024177.1.NonCompliant.ihs.out CM024178.1.NonCompliant.ihs.out CM024179.1.NonCompliant.ihs.out CM024180.1.NonCompliant.ihs.out CM024181.1.NonCompliant.ihs.out CM024182.1.NonCompliant.ihs.out CM024183.1.NonCompliant.ihs.out CM024184.1.NonCompliant.ihs.out CM024185.1.NonCompliant.ihs.out CM024186.1.NonCompliant.ihs.out CM024188.1.NonCompliant.ihs.out
CM024189.1.NonCompliant.ihs.out CM024189.1.NonCompliant.ihs.out
CM024190.1.NonCompliant.ihs.out CM024191.1.NonCompliant.ihs.out



# Add chromosome names to normalized IHS results
parallel --dry-run "sed -e 's/^/{1}\t/g' {1}.{2}.ihs.out.100bins.norm > {1}.{2}.ihs.out.100bins.norm.temp" ::: CM024168.1, CM024169.1, CM024170.1, CM024171.1, CM024172.1, CM024173.1, CM024174.1, CM024175.1, CM024176.1, CM024177.1, CM024178.1, CM024179.1, CM024180.1, CM024181.1, CM024182.1, CM024183.1, CM024184.1, CM024185.1, CM024186.1, CM024187.1, CM024188.1, CM024189.1, CM024190.1, CM024191.1::: NonComnpliant FullyCompliant

cat *.FullyCompliant.ihs.out.100bins.norm.temp | sed 's/As_//g' > ALL.FullyCompliant.IHS.ihs.out.100bins.norm.txt
cat *.NonCompliant.ihs.out.norm.temp | sed 's/As_//g' > ALL.NonCompliant.IHS.ihs.out.100bins.norm.txt
```
### Calculate XP-EHH between districts (per chromosome)
```
parallel --dry-run "selscan --xpehh --vcf As_POPGEN.biallelic.{}.FullyCompliant.vcf --vcf-ref As_POPGEN.biallelic.{}.NonCompliant.vcf --map {}.gmap --threads 2 --out {}" ::: CM024168.1, CM024169.1, CM024170.1, CM024171.1, CM024172.1, CM024173.1, CM024174.1, CM024175.1, CM024176.1, CM024177.1, CM024178.1, CM024179.1, CM024180.1, CM024181.1, CM024182.1, CM024183.1, CM024184.1, CM024185.1, CM024186.1, CM024187.1, CM024188.1, CM024189.1, CM024190.1, CM024191.1

#Normalize
norm --xpehh --files CM024168.1.xpehh.out CM024169.1.xpehh.out CM024170.1.xpehh.out CM024171.1.xpehh.out CM024172.1.xpehh.out CM024173.1.xpehh.out CM024174.1.xpehh.out CM024175.1.xpehh.out CM024176.1.xpehh.out CM024177.1.xpehh.out CM024178.1.xpehh.out CM024179.1.xpehh.out CM024180.1.xpehh.out CM024181.1.xpehh.out CM024182.1.xpehh.out CM024183.1.xpehh.out CM024184.1.xpehh.out CM024185.1.xpehh.out CM024186.1.xpehh.out CM024187.1.xpehh.out CM024187.1.xpehh.out CM024188.1.xpehh.out CM024189.1.xpehh.out CM024190.1.xpehh.out CM024191.1.xpehh.out

# Create header for XP-EHH results
head -1 CM024168.1.xpehh.out.norm | sed -e '1s/id/chromosome\tid/g' 
head -1 CM024169.1.xpehh.out.norm | sed -e '1s/id/chromosome\tid/g'  > header.xpehh
```
# Add chromosome names to normalised XP-EHH results
parallel --dry-run "sed -e 's/^/{1}\t/g' {}.xpehh.out.norm > {}.xpehh.out.norm.temp" ::: CM024168.1, CM024169.1, CM024170.1, CM024171.1, CM024172.1, CM024173.1, CM024174.1, CM024175.1, CM024176.1, CM024177.1, CM024178.1, CM024179.1, CM024180.1, CM024181.1, CM024182.1, CM024183.1, CM024184.1, CM024185.1, CM024186.1, CM024187.1, CM024188.1, CM024189.1, CM024190.1, CM024191.1

cat *.xpehh.out.norm.temp > ALL.FullyvsNonCompliant.xpehh.xpehh.out.norm.txt
```
### Calculate F<sub>ST</sub> between districts (per chromosome)
```
# Run PIXY (described above too, the specific command is)
parallel --dry-run "pixy --stats fst 
--vcf As_POPGEN.allsites.{1}.vcf 
--zarr_path zarr/ 
--window_size 2000 # Or 5kb where needed
--reuse_zarr yes
--populations Household.list 
--variant_filter_expression 'DP>=10,GQ>=20,RGQ>=20' 
--invariant_filter_expression 'DP>=10,RGQ>=20' 
--outfile_prefix output/pixy.{1}.25000.district" ::: CM024168.1, CM024169.1, CM024170.1, CM024171.1, CM024172.1, CM024173.1, CM024174.1, CM024175.1, CM024176.1, CM024177.1, CM024178.1, CM024179.1, CM024180.1, CM024181.1, CM024182.1, CM024183.1, CM024184.1, CM024185.1, CM024186.1, CM024187.1, CM024188.1, CM024189.1, CM024190.1, CM024191.1


cat pixy.As_*.2000.district_fst.txt | grep -v pop | cat fst.header - > fst.district.txt 

# RUN VCFtools (alternative to PIXY for FST calculations) for between household comparisons

vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FINAL/As_POPGEN.vcf --weir-fst-pop hh_cluster1 --weir-fst-pop hh_cluster2 --fst-window-size 2000 --out HH_2000.windowed.weir.txt

```
___
```
```
```
### Tajima's D
```
# Calculate Tajima's D for each treatment subpopulation in 2 kb windows. 
vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FINAL/As_POPGEN.vcf --keep FullyCompliant.list --TajimaD 2000 --out FullyCompliant_TAJIMA_D_2000
vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FINAL/As_POPGEN.vcf --keep SemiCompliant.list --TajimaD 2000 --out SemiCompliant_TAJIMA_D_2000
vcftools --vcf ${WORKING_DIR}/06_ANALYSIS/FINAL/As_POPGEN.vcf --keep NonCompliant.list --TajimaD 2000 --out NonCompliant_TAJIMA_D_2000

cat FullyCompliant_TAJIMA_D_2000.Tajima.D | grep -v Taj | awk '{print $1,$2,$3,$4,"Kocoge"}' | sed 's/ / /g' >> FullyCompliant_TAJIMA_D_2000.Tajima.D.2kb.txt
cat SemiCompliant_TAJIMA_D_2000.Tajima.D | grep -v Taj | awk '{print $1,$2,$3,$4,"Bugoto"}' | sed 's/ / /g' >> SemiCompliant_TAJIMA_D_2000.Tajima.D.2kb.txt
cat NonCompliant_TAJIMA_D_2000.Tajima.D | grep -v Taj | awk '{print $1,$2,$3,$4,"Bwondha"}' | sed 's/ / /g' >> SemiCompliant_TAJIMA_D_2000.Tajima.D 2kb.txt



library(rehh)
library(ggplot2)
library(patchwork)
library(wesanderson)
pal <- wes_palette("Darjeeling1")

unique(combined$CHR)

combined$CHR <- as.factor(combined$CHR)
###### BEGIN REHH ANALYSIS 

#1

haplohhData <- data2haplohh("Alum.biallelic.1.vcf.recode.vcf_phased.vcf.gz", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = "CM024168.1")
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs1 <- ihh2ihs(scan, freqbin = 1)
distribplot(scanihs1$ihs$IHS, xlab = "iHS")

distribplot(scanihs1$ihs$IHS, 
            xlab = "iHS", 
            qqplot = TRUE)

manhattanplot(scanihs1, pval = TRUE, threshold = 4,
              main = "iHS (A.Lumbricoides Chr 1)")

ggplot(combined, aes(x = POSITION, y = LOGPVALUE, colour = combined$CHR)) + geom_point()




#2

haplohhData <- data2haplohh("Alum.biallelic.2.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = "CM024169.1")
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs2 <- ihh2ihs(scan, freqbin = 1)
distribplot(scanihs2$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs2$ihs$IHS, xlab = "iHS")


manhattanplot(scanihs2, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 2)")
ggplot(scanihs2$ihs, aes(POSITION, IHS)) + geom_point() + scale_y_continuous(expand = c(0, 0) )

plot1 <- ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)
plot2 <- ggplot(scanihs2$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 3
haplohhData <- data2haplohh("Alum.biallelic.3.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = "CM024170.1")
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs3 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs3$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs3$ihs$IHS, xlab = "iHS")


manhattanplot(scanihs3, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 3)")
ggplot(scanihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 4
haplohhData <- data2haplohh("Alum.biallelic.4.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = "CM024171.1")
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs4 <- ihh2ihs(scan, freqbin = 1)
distribplot(scanihs4$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs4$ihs$IHS, xlab = "iHS")


manhattanplot(scanihs4, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 4)")

ggplot(scanihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 5
haplohhData <- data2haplohh("Alum.biallelic.5.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = "CM024172.1")

hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs5 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs6$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs6$ihs$IHS, xlab = "iHS")


manhattanplot(scanihs5, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 5)")

ggplot(scanihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 6
haplohhData <- data2haplohh("Alum.biallelic.6.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, "CM024173.1")
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs6 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs6$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs6$ihs$IHS, xlab = "iHS")


manhattanplot(scanihs6, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 6)")

ggplot(scanihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 7
haplohhData <- data2haplohh("Alum.biallelic.7.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = 'CM024174.1' )
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs7 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs7$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs7$ihs$IHS, xlab = "iHS")


manhattanplot(scanihs7, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 7)")

ggplot(scanihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 8
haplohhData <- data2haplohh("Alum.biallelic.8.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = 'CM024175.1')
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs8 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs8$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs8$ihs$IHS, xlab = "iHS")


manhattanplot(scanihs8, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 8)")

ggplot(scanihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 9
haplohhData <- data2haplohh("Alum.biallelic.9.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = 'CM024176.1')
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs9 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs8$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs8$ihs$IHS, xlab = "iHS")


manhattanplot(scanihs9, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 9)")
#Alt Plots
ggplot(scanihs9$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 10
haplohhData <- data2haplohh("Alum.biallelic.10.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = 'CM024177.1')
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs10 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs10$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs10$ihs$IHS, xlab = "iHS")


manhattanplot(scanihs10, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 10)")

#Alt plots
ggplot(scanihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 11
haplohhData <- data2haplohh("Alum.biallelic.11.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = 'CM024178.1')
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs11 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs11$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs11$ihs$IHS, xlab = "iHS")


manhattanplot(scanihs11, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 11)")


ggplot(scanihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 12
haplohhData <- data2haplohh("Alum.biallelic.12.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = 'CM024179.1')
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs12 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs12$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs12$ihs$IHS, xlab = "iHS")


manhattanplot(scanihs12, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 12)")


ggplot(scanihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 13
haplohhData <- data2haplohh("Alum.biallelic.13.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = 'CM024180.1')
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs13 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs13$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs13$ihs$IHS, xlab = "iHS")


manhattanplot(scanihs13, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 13)")


ggplot(scanihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 14
haplohhData <- data2haplohh("Alum.biallelic.14.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = 'CM024181.1')
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs14 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs14$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs14$ihs$IHS, xlab = "iHS")


manhattanplot(scanihs14, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 14)")


ggplot(scanihs14$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs14$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 15
haplohhData <- data2haplohh("Alum.biallelic.15.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = 'CM024182.1')
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs15 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs15$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs15$ihs$IHS, xlab = "iHS")

manhattanplot(scanihs15, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 15)")

ggplot(scanihs15$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 16
haplohhData <- data2haplohh("Alum.biallelic.16.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = 'CM024183.1')
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs16 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs16$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs16$ihs$IHS, xlab = "iHS")

manhattanplot(scanihs16, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 16)")

ggplot(scanihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 17
haplohhData <- data2haplohh("Alum.biallelic.17.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = 'CM024184.1')
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs17 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs17$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs17$ihs$IHS, xlab = "iHS")

manhattanplot(scanihs17, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 17)")

ggplot(scanihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 18
haplohhData <- data2haplohh("Alum.biallelic.18.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = 'CM024185.1')
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs18 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs18$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs18$ihs$IHS, xlab = "iHS")

manhattanplot(scanihs18, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 18)")


ggplot(scanihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

# 19
haplohhData <- data2haplohh("Alum.biallelic.19.vcf.recode.vcf_phased.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = 'CM024186.1')
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihs19 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihs19$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihs19$ihs$IHS, xlab = "iHS")

manhattanplot(scanihs19, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr 19)")

ggplot(scanihs$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

#X1

haplohhData <- data2haplohh("RmDupAlumX1.phased.haps.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = "CM024187.1")
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihsX1 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihsX1$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihsX1$ihs$IHS, xlab = "iHS")

manhattanplot(scanihsX1, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr X1)")

ggplot(scanihsX1$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihsX1$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

#X2

haplohhData <- data2haplohh("RmDupAlumX2.phased.haps.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = "CM024188.1")
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihsX2 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihsX1$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihsX1$ihs$IHS, xlab = "iHS")

manhattanplot(scanihsX2, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr X1)")

ggplot(scanihsX1$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihsX1$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

#X3

haplohhData <- data2haplohh("RmDupAlumX3.phased.haps.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = "CM024189.1")
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihsX3 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihsX3$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihsX3$ihs$IHS, xlab = "iHS")

manhattanplot(scanihsX3, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr X3)")

ggplot(scanihsX1$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihsX1$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

#X4

haplohhData <- data2haplohh("RmDupAlumX4.phased.haps.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = "CM024190.1")
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihsX4 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihsX1$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihsX1$ihs$IHS, xlab = "iHS")

manhattanplot(scanihsX4, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr X3)")

ggplot(scanihsX1$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihsX1$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)

#X5

haplohhData <- data2haplohh("RmDupAlumX5.phased.haps.vcf", polarize_vcf = FALSE, remove_multiple_markers = TRUE, chr.name = "CM024191.1")
hap_f <- subset(haplohhData, min_maf = 0.05)

scan <- scan_hh(hap_f, polarized = FALSE)
scanihsX5 <- ihh2ihs(scan, freqbin = 1)

distribplot(scanihsX1$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(scanihsX1$ihs$IHS, xlab = "iHS")

manhattanplot(scanihsX5, pval = TRUE, threshold = 4, main = "iHS (A.Lumbricoides Chr X5)")

ggplot(scanihsX1$ihs, aes(POSITION, IHS)) + geom_point()

ggplot(scanihsX1$ihs, aes(POSITION, LOGPVALUE)) + geom_point() + ylim(0,7)


####### Plot GWAS ---------------

df1 <- as.data.frame(scanihs1$ihs)
df2 <- as.data.frame(scanihs2$ihs)
df3 <- as.data.frame(scanihs3$ihs)
df4 <- as.data.frame(scanihs4$ihs)
df5 <- as.data.frame(scanihs5$ihs)
df6 <- as.data.frame(scanihs6$ihs)
df7 <- as.data.frame(scanihs7$ihs)
df8 <- as.data.frame(scanihs8$ihs)
df9 <- as.data.frame(scanihs9$ihs)
df10 <- as.data.frame(scanihs10$ihs)
df11 <- as.data.frame(scanihs11$ihs)
df12 <- as.data.frame(scanihs12$ihs)
df13 <- as.data.frame(scanihs13$ihs)
df14 <- as.data.frame(scanihs14$ihs)
df15 <- as.data.frame(scanihs15$ihs)
df16 <- as.data.frame(scanihs16$ihs)
df17 <- as.data.frame(scanihs17$ihs)
df18 <- as.data.frame(scanihs18$ihs)
df19 <- as.data.frame(scanihs19$ihs)


dfX1 <- as.data.frame(scanihsX1$ihs)
dfX2 <- as.data.frame(scanihsX2$ihs)
dfX3 <- as.data.frame(scanihsX3$ihs)
dfX4 <- as.data.frame(scanihsX4$ihs)
dfX5 <- as.data.frame(scanihsX5$ihs)

#Autosomes

df_list <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17, df18, df19)
combined <- Reduce(function(x, y) merge(x, y, all = TRUE), df_list)

#SexChromosomes

df_list <- list(dfX1, dfX2, dfX3, dfX4, dfX5)
combined <- Reduce(function(x, y) merge(x, y, all = TRUE), df_list)

unique(combined$CHR)

combined$CHR <- as.factor(combined$CHR)

manhattanplot(combined, main = "iHS (A.Lumbricoides Chr X1-X5)")



###### Calcucating canditate regions ------

cr.cgu <- calc_candidate_regions(combined,
                                 threshold = 4,
                                 pval = TRUE,
                                 window_size = 1E6,
                                 overlap = 1E5,
                                 min_n_extr_mrk = 2)
cr.cgu

###### Plotting the results ----------
dev.new(width=20, height=10, unit="in", noRStudioGD=TRUE)

plot <- manhattanplot(combined, 
              pval = TRUE,
              threshold = 4, cr = cr.cgu)




plot <- manhattanplot(combined, 
                      pval = TRUE,
                      threshold = 4, cr = cr.cgu, inset = 0)

ihs <- combined
# create new data frame
wgscan.cgu.ihs.qqman <- data.frame(
  CHR = as.integer(factor(ihs$CHR, 
                          levels = unique(ihs$CHR))),
  # chromosomes as integers
  BP = ihs$POSITION,         # base pairs
  P = 10**(-ihs$LOGPVALUE),  # transform back to p-values
  SNP = row.names(ihs)       # SNP names
)

library(qqman)

manhattan(wgscan.cgu.ihs.qqman, chr = 'CHR', bp='BP', p='P', snp='SNP', ylim=c(0,7), suggestiveline = 4, col = pal, annotatePval = 0.01)



install.packages("qqman")
library(qqman)
library(wesanderson)
library(tidyverse)

# Load necessary libraries
library(data.table)

pal <- wes_palette("Darjeeling1")
##### Plotting Fst -------
fst_data <- read.table("Pi10kTotal.windowed.pi", sep='\t' ,header = TRUE)
unique(fst_data$CHROM)

fst_data$CHROM[fst_data$CHROM == "CM024168.1"] <- "1"
fst_data$CHROM[fst_data$CHROM == "CM024169.1"] <- "2"
fst_data$CHROM[fst_data$CHROM == "CM024170.1"] <- "3"
fst_data$CHROM[fst_data$CHROM == "CM024171.1"] <- "4"
fst_data$CHROM[fst_data$CHROM == "CM024172.1"] <- "5"
fst_data$CHROM[fst_data$CHROM == "CM024173.1"] <- "6"
fst_data$CHROM[fst_data$CHROM == "CM024174.1"] <- "7"
fst_data$CHROM[fst_data$CHROM == "CM024175.1"] <- "8"
fst_data$CHROM[fst_data$CHROM == "CM024176.1"] <- "9"
fst_data$CHROM[fst_data$CHROM == "CM024177.1"] <- "10"
fst_data$CHROM[fst_data$CHROM == "CM024178.1"] <- "11"
fst_data$CHROM[fst_data$CHROM == "CM024179.1"] <- "12"
fst_data$CHROM[fst_data$CHROM == "CM024180.1"] <- "13"
fst_data$CHROM[fst_data$CHROM == "CM024181.1"] <- "14"
fst_data$CHROM[fst_data$CHROM == "CM024182.1"] <- "15"
fst_data$CHROM[fst_data$CHROM == "CM024183.1"] <- "16"
fst_data$CHROM[fst_data$CHROM == "CM024184.1"] <- "17"
fst_data$CHROM[fst_data$CHROM == "CM024185.1"] <- "18"
fst_data$CHROM[fst_data$CHROM == "CM024186.1"] <- "19"
fst_data$CHROM[fst_data$CHROM == "CM024187.1"] <- "20"
fst_data$CHROM[fst_data$CHROM == "CM024188.1"] <- "21"
fst_data$CHROM[fst_data$CHROM == "CM024189.1"] <- "22"
fst_data$CHROM[fst_data$CHROM == "CM024190.1"] <- "23"
fst_data$CHROM[fst_data$CHROM == "CM024191.1"] <- "24"


colnames(fst_data)
length_ <- dim(fst_data)[1]
length_
fst_data$SNP <- paste('SNP',1:length_)
head(fst_data)
fst_data$CHROM <- as.numeric(fst_data$CHROM)
unique(fst_data$CHROM)

dev.new(width=20, height=10, unit="in", noRStudioGD=TRUE)

manhattan(fst_data, chr = 'CHROM', bp= 'BIN_START', p = 'PI', snp = 'SNP', logp = FALSE, ylim=c(0, 0.2), col = pal)

manhattan(fst_data)
str(fst)

#HIGHLIGHT OUTLIERS

quantile(fst$WEIGHTED_FST, c(0.975, 0.995), na.rm = T)
my_threshold <- quantile(fst$WEIGHTED_FST, 0.975, na.rm = T)
fst <- fst %>% mutate(outlier = ifelse(WEIGHTED_FST > my_threshold, "outlier", "background"))
fst %>% group_by(outlier) %>% tally()

##### Plotting Tajimas D -------
Taj <- read.table("GenomScan_10kb.Tajima.D", sep='\t' ,header = TRUE)

colnames(Taj)
length_ <- dim(Taj)[1]
length_
Taj$SNP <- paste('SNP',1:length_)


manhattan(fst, chr = 'CHROM', bp= 'BIN_START', p = 'MEAN_FST', snp = 'SNP', col = pal, logp = TRUE, ylab = 'TajimaD', xlab = 'CHR', ylim=c(-0.8, 2))

manhattan(fst)
str(fst)


windowed_fst <- read.table("GenomeScanFst10k.weir.fst.edit", sep="\t", header=TRUE)
str(windowed_fst)

quantile(windowed_fst$WEIGHTED_FST, probs = c(.95, .99, .999))

pdf("fst_starlings_windowed.pdf", width=10, height=5)

ggplot(windowed_fst, aes(x=X1, y=WEIGHTED_FST)) + 
  geom_point() + 
  geom_hline(yintercept=0.35, linetype="dashed", color = "red") +
  labs(x = "Window Number") +
  theme_classic()

dev.off()

q()


####### Plotting Pi ------------

fst_data <- read.table("Pi10kTotal.windowed.pi", sep='\t' ,header = TRUE)
unique(fst_data$CHROM)
fst_data$CHROM[fst_data$CHROM == "CM024168.1"] <- "1"
fst_data$CHROM[fst_data$CHROM == "CM024169.1"] <- "2"
fst_data$CHROM[fst_data$CHROM == "CM024170.1"] <- "3"
fst_data$CHROM[fst_data$CHROM == "CM024171.1"] <- "4"
fst_data$CHROM[fst_data$CHROM == "CM024172.1"] <- "5"
fst_data$CHROM[fst_data$CHROM == "CM024173.1"] <- "6"
fst_data$CHROM[fst_data$CHROM == "CM024174.1"] <- "7"
fst_data$CHROM[fst_data$CHROM == "CM024175.1"] <- "8"
fst_data$CHROM[fst_data$CHROM == "CM024176.1"] <- "9"
fst_data$CHROM[fst_data$CHROM == "CM024177.1"] <- "10"
fst_data$CHROM[fst_data$CHROM == "CM024178.1"] <- "11"
fst_data$CHROM[fst_data$CHROM == "CM024179.1"] <- "12"
fst_data$CHROM[fst_data$CHROM == "CM024180.1"] <- "13"
fst_data$CHROM[fst_data$CHROM == "CM024181.1"] <- "14"
fst_data$CHROM[fst_data$CHROM == "CM024182.1"] <- "15"
fst_data$CHROM[fst_data$CHROM == "CM024183.1"] <- "16"
fst_data$CHROM[fst_data$CHROM == "CM024184.1"] <- "17"
fst_data$CHROM[fst_data$CHROM == "CM024185.1"] <- "18"
fst_data$CHROM[fst_data$CHROM == "CM024186.1"] <- "19"
fst_data$CHROM[fst_data$CHROM == "CM024187.1"] <- "20"
fst_data$CHROM[fst_data$CHROM == "CM024188.1"] <- "21"
fst_data$CHROM[fst_data$CHROM == "CM024189.1"] <- "22"
fst_data$CHROM[fst_data$CHROM == "CM024190.1"] <- "23"
fst_data$CHROM[fst_data$CHROM == "CM024191.1"] <- "24"

colnames(fst_data)
length_ <- dim(fst_data)[1]
length_
fst_data$SNP <- paste('SNP',1:length_)
head(fst_data)
fst_data$CHROM <- as.numeric(fst_data$CHROM)
unique(fst_data$CHROM)

manhattan(fst_data, chr = 'CHROM', bp= 'BIN_START', p = 'PI', snp = 'SNP', logp = FALSE, ylim=c(-0.02, 0.02), col = pal)


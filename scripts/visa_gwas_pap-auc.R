# Robert A Petit III
# Made using RStudio

library(ggplot2)
source('~/repos/visa-gwas/scripts/visa_gwas_functions.R')

###############################################################################
SNP <- "~/repos/visa-gwas/data/genotype.txt"
MIC <- "~/repos/visa-gwas/data/phenotype_pap-auc.txt"

GENOME_SIZE <- 2814816
MAF <- 0.05
MIC_THRESHOLD <- 0.9

OUT_DIR <- "~/gwas/pap-auc/"
RT_DIR <- paste(OUT_DIR, 'roadtrips/', sep='')
QRT_DIR <- paste(OUT_DIR, 'qroadtrips/', sep='')

###############################################################################
Sys.time()
system(paste("rm -rf ", OUT_DIR, sep=''))

# Read input & sort by strain
mic <- read.table(MIC, header=TRUE)
mic <- mic[ order(mic[,1]), ]

snps <- read.table(SNP, header=TRUE, comment.char="", nrows=75)
snps <- snps[ order(snps[,1]), ]

###############################################################################
# Basic SNP stats (names and minor allele frequency)
snp_stats <- data.frame(name=colnames(snps)[2:ncol(snps)],
                        freq=colSums(snps[2:ncol(snps)]) / nrow(snps))

bonferroni_correction <- -1 * log10(0.05 / nrow(snp_stats))
###############################################################################
# ROADTRIPS
write_roadtrips(snps, mic, MIC_THRESHOLD, RT_DIR)
rt <- run_roadtrips(RT_DIR)
rt <- filter_results(rt, snp_stats, MAF)

# Manhattan plot & QQplot
manhattan(rt, bonferroni_correction, GENOME_SIZE)
qqplot(rt[with(rt, order(logp)), ], bonferroni_correction)

###############################################################################
# QROADTRIPS
write_qtrips(snps, mic, QRT_DIR)
qt <- run_qtrips(QRT_DIR)
qt <- filter_results(qt, snp_stats, MAF)

# Manhattan plot & QQplot
manhattan(qt, bonferroni_correction, GENOME_SIZE)
qqplot(qt[with(qt, order(logp)), ], bonferroni_correction)

Sys.time()

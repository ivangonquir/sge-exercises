###############################################################################
# BSG - Statistical Genetics - Practical 1  (answers only, base R)
# Run it line by line (Ctrl+Enter) and read the results in the console.
###############################################################################

############################ SNP DATASET ######################################

cat("\014")
rm(list = ls())
setwd("G:/My Drive/University/Masters/MDS/q3/BSG/sge-exercises")

# Set the working directory to the folder with the data file, or use file.choose()
raw <- read.table("TSICHR22RAW.raw", header = TRUE)
X   <- raw[, 7:ncol(raw)]                 # only the genetic information
dim(X)                                    # individuals, variants

## SNP 1: number of variants and % missing
ncol(X)
100 * mean(is.na(X))

## SNP 2: monomorphic variants
nB   <- colSums(X, na.rm = TRUE)          # copies of allele B
nTot <- 2 * colSums(!is.na(X))            # total alleles typed
mono <- nB == 0 | nB == nTot              # only one allele present
sum(mono)
100 * mean(mono)                          # % monomorphic
X <- X[, !mono]                           # remove them
ncol(X)                                   # variants remaining

## SNP 3: rs8138488_C
g  <- X[, "rs8138488_C"]
gc <- table(factor(g, levels = 0:2, labels = c("AA", "AB", "BB")))
gc                                        # genotype counts
count_A <- 2 * gc[["AA"]] + gc[["AB"]]
count_B <- 2 * gc[["BB"]] + gc[["AB"]]        # B = the C allele
min(count_A, count_B)                     # minor allele count
min(count_A, count_B) / (count_A + count_B)   # MAF

## SNP 4: MAF of all markers
p   <- colSums(X, na.rm = TRUE) / (2 * colSums(!is.na(X)))
maf <- pmin(p, 1 - p)
hist(maf, breaks = 50, col = "steelblue",
     main = "Minor allele frequency of SNPs (chr22, TSI)",
     xlab = "Minor allele frequency", ylab = "Number of SNPs")
100 * mean(maf < 0.05)                    # % MAF < 0.05
100 * mean(maf < 0.01)                    # % MAF < 0.01

## SNP 5: observed heterozygosity
Ho <- colMeans(X == 1, na.rm = TRUE)
hist(Ho, breaks = 50, col = "darkorange",
     main = "Observed heterozygosity of SNPs",
     xlab = "Observed heterozygosity (Ho)", ylab = "Number of SNPs")
range(Ho)

## SNP 6: expected heterozygosity
He <- 1 - (p^2 + (1 - p)^2)
hist(He, breaks = 50, col = "seagreen",
     main = "Expected heterozygosity of SNPs",
     xlab = "Expected heterozygosity (He)", ylab = "Number of SNPs")
range(He)
mean(He)


############################ STR DATASET ######################################

library(HardyWeinberg)
data(NistSTRs)
# If HardyWeinberg can't be installed: load("NistSTRs.rda")  (see notes)

## STR 1: individuals and STRs
nrow(NistSTRs)                            # individuals
ncol(NistSTRs) / 2                        # STRs (2 columns per STR)

K <- ncol(NistSTRs) / 2
str_names <- colnames(NistSTRs)[seq(1, ncol(NistSTRs), by = 2)]

# alleles of STR k as text, so 14.3 is different from 14 and 15
alleles_of <- function(k) {
  a <- c(as.character(NistSTRs[, 2*k - 1]), as.character(NistSTRs[, 2*k]))
  a[!is.na(a)]
}

## STR 2: number of alleles per STR
n_alleles <- function(k) length(unique(alleles_of(k)))
nA <- sapply(1:K, n_alleles)
names(nA) <- str_names
nA
mean(nA); sd(nA); median(nA); min(nA); max(nA)

## STR 3: table and barplot
tab <- table(nA)
tab
barplot(tab, col = "slateblue",
        main = "Number of STRs by number of alleles",
        xlab = "Number of alleles", ylab = "Number of STRs")
names(tab)[tab == max(tab)]               # most common number of alleles

## STR 4: expected heterozygosity
He_str <- sapply(1:K, function(k) {
  f <- table(alleles_of(k)) / length(alleles_of(k))
  1 - sum(f^2)
})
names(He_str) <- str_names
hist(He_str, col = "seagreen",
     main = "Expected heterozygosity of STRs",
     xlab = "Expected heterozygosity (He)", ylab = "Number of STRs")
mean(He_str)

## STR 5: observed heterozygosity and plot
Ho_str <- sapply(1:K, function(k) {
  a1 <- as.character(NistSTRs[, 2*k - 1])
  a2 <- as.character(NistSTRs[, 2*k])
  mean(a1 != a2, na.rm = TRUE)
})
names(Ho_str) <- str_names
plot(He_str, Ho_str, pch = 19,
     xlim = range(c(He_str, Ho_str)), ylim = range(c(He_str, Ho_str)),
     main = "Observed vs expected heterozygosity (STRs)",
     xlab = "Expected heterozygosity (He)", ylab = "Observed heterozygosity (Ho)")
abline(0, 1, lty = 2, col = "red")        # line Ho = He
cor(Ho_str, He_str)

## STR 6: comparison SNPs vs STRs
c(SNP_mean_He = mean(He), STR_mean_He = mean(He_str))
boxplot(list(SNPs = He, STRs = He_str), col = c("seagreen", "slateblue"),
        main = "Expected heterozygosity: SNPs vs STRs",
        xlab = "Marker type", ylab = "Expected heterozygosity (He)")
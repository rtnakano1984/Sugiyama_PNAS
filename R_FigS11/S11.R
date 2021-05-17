# Sugiyama et al., PNAS, 2021
# Fig S11
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

# cleanup
rm(list=ls())

# initialization
library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/path_to_data_dir/R_FigS11/"
source("/path_to_data_dir/function.R")

# import
dat_all <- read.table(paste(data_dir, "R_Fig S11.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat_all$batch <- factor(dat_all$batch)

genotype <- c("wt", "bglu28bglu30")
dat_all$genotype <- factor(dat_all$genotype, levels=genotype)

sulfur <- c("S1500", "S150", "S15", "S0")
dat_all$sulfur <- factor(dat_all$sulfur, levels=sulfur)

# empty df for estimates
estim_df <- data.frame(group=NULL, estim=NULL, cpd=NULL)


# individual glucosinolate species
for(cpd in unique(dat_all$cpd)){
  idx <- dat_all$cpd == cpd
  dat <- dat_all[idx,]

  pdf(paste(data_dir, cpd, "-data_nomral.pdf", sep=""))
    hist(dat$value, breaks=20)
    qqnorm(dat$value)
    qqline(dat$value, col='red')

    hist(log10(dat$value), breaks=20)
    qqnorm(log10(dat$value))
    qqline(log10(dat$value), col='red')

    hist(sqrt(dat$value), breaks=20)
    qqnorm(sqrt(dat$value))
    qqline(sqrt(dat$value), col='red')
  dev.off()

  lmer_fit <- lmer(sqrt(value) ~ sulfur:genotype - 1 + (1|batch), dat)
  dat_lmer <- summary(lmer_fit)

  sink(paste(data_dir, cpd, "-lmer_summary.txt", sep=""))
    print(dat_lmer)
  sink()

  pdf(paste(data_dir, cpd, "-diagnosis.pdf", sep=""))
    plot(fitted(lmer_fit), resid(lmer_fit))
    abline(0, 0, col="red")

    qqnorm(dat_lmer$resid)
    qqline(dat_lmer$resid, col="red")
  dev.off()

  adj_p.val <- pairwise_ttest(dat_lmer)
  p.letters <- multcompLetters(adj_p.val)

  names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "sulfur", "")
  names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "genotype", "")

  sort <- apply(expand.grid(genotype, sulfur)[, 2:1], 1, function(x) paste(x, collapse=":"))
  p.letters_df <- as.data.frame(p.letters$Letters[sort])

  write.table(p.letters_df, file=paste(data_dir, cpd, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}



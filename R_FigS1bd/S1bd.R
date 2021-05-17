# Sugiyama et al., PNAS, 2021
# Fig S1A,C
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

# cleanup
rm(list=ls())

# initialization
library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/path_to_data_dir/R_FigS1bd/"
source("/path_to_data_dir/function.R")

# import
dat <- read.table(paste(data_dir, "R_Fig S1bd.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$genotype <- str_replace(dat$genotype, " ", "")
dat$batch <- factor(dat$batch)

sulfur <- c("S150", "S0", "4MTB", "allyl", "PhE", "I3G")
dat$sulfur <- factor(dat$sulfur, levels=sulfur)

# for each compound
for(cpd in unique(dat$cpd)){
  idx <- dat$cpd == cpd
  temp <- dat[idx,]

  # distribution
  pdf(paste(data_dir, cpd, "-data_nomral.pdf", sep=""))
    hist(temp$value, breaks=20)
    qqnorm(temp$value)
    qqline(temp$value, col='red')

    hist(log10(temp$value), breaks=20)
    qqnorm(log10(temp$value))
    qqline(log10(temp$value), col='red')
  dev.off()

  # fit
  lmer_fit <- lmer(value ~ sulfur - 1 + (1|batch), temp)
  dat_lmer <- summary(lmer_fit)

  sink(paste(data_dir, cpd, "-lmer_summary.txt", sep=""))
  print(dat_lmer)
  sink()

  # diagnosis
  pdf(paste(data_dir, cpd, "-diagnosis.pdf", sep=""))
    plot(fitted(lmer_fit), resid(lmer_fit))
    abline(0, 0, col="red")

    qqnorm(dat_lmer$resid)
    qqline(dat_lmer$resid, col="red")
  dev.off()

  # pairwise test
  adj_p.val <- pairwise_ttest(dat_lmer)
  p.letters <- multcompLetters(adj_p.val)

  names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "sulfur", "")
  sort <- sulfur

  write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, cpd, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}



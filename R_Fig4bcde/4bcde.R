# Sugiyama et al., PNAS, 2021
# Fig 4B-E
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

# cleanup
rm(list=ls())

# initialization
library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/path_to_data_dir/R_Fig3b/"
source("/path_to_data_dir/function.R")

# import
dat_all <- read.table(paste(data_dir, "R_Fig 4bcde.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat_all$batch <- factor(dat_all$batch)

# for each treatment
for(treatment in unique(dat_all$treatment)){
  idx <- dat_all$treatment == treatment
  dat <- dat_all[idx,]

  # distribution
  pdf(paste(data_dir, treatment, "-data_nomral.pdf", sep=""))
    hist(dat$value, breaks=20)
    qqnorm(dat$value)
    qqline(dat$value, col='red')

    hist(log10(dat$value), breaks=20)
    qqnorm(log10(dat$value))
    qqline(log10(dat$value), col='red')
  dev.off()

  # fit
  lmer_fit <- lmer(log10(value) ~ time:genotype - 1 + (1|batch), dat)
  dat_lmer <- summary(lmer_fit)

  sink(paste(data_dir, treatment, "-lmer_summary.txt", sep=""))
  print(dat_lmer)
  sink()

  # diagnosis
  pdf(paste(data_dir, treatment, "-diagnosis.pdf", sep=""))
    plot(fitted(lmer_fit), resid(lmer_fit))
    abline(0, 0, col="red")

    qqnorm(dat_lmer$resid)
    qqline(dat_lmer$resid, col="red")
  dev.off()

  # pairwise test
  adj_p.val <- pairwise_ttest(dat_lmer)
  p.letters <- multcompLetters(adj_p.val)

  names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "time", "")
  names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "genotype", "")

  sort <- apply(expand.grid(c("wt", "oxp1"), c("00h", "09h", "24h", "48h"))[, 2:1], 1, function(x) paste(x, collapse=":"))
  write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, treatment, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}


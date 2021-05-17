# Sugiyama et al., PNAS, 2021
# Fig 5E and S12
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

# cleanup
rm(list=ls())

# initialization
library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/path_to_data_dir/R_FigS12/"
source("/path_to_data_dir/function.R")

# import
dat_all <- read.table(paste(data_dir, "R_Fig S12.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat_all$batch <- factor(dat_all$batch)

genotype <- c("Col", "bglu2830", "tgg12", "tgg45", "pen2", "pyk10", "bglu18")
dat_all$genotype <- factor(dat_all$genotype, levels=genotype)

sulfur <- c("S1500", "S15")
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

  # logFC
  estim <- dat_lmer$coefficients[,1]
  names(estim) <- str_replace_all(names(estim), "sulfur", "")
  names(estim) <- str_replace_all(names(estim), "genotype", "")
  idx <- match(sort, names(estim))
  estim <- estim[idx]
  estim_df_temp <- data.frame(group=names(estim), estim=estim^2, cpd=cpd, row.names=NULL, letters=p.letters_df[,1], stringsAsFactors=F)
  estim_df_temp$logFC <- log2(estim_df_temp$estim / estim_df_temp$estim[estim_df_temp$group == "S1500:Col"])
  estim_df <- rbind(estim_df, estim_df_temp)

}

write.table(estim_df, file=paste(data_dir, "logFC.txt", sep=""), sep="\t", row.names=F, col.names=T, quote=F)



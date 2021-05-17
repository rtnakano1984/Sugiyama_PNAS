# Sugiyama et al., PNAS, 2021
# Fig 5D
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

# cleanup
rm(list=ls())

# initialization
library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/path_to_data_dir/R_Fig5d/"
source("/path_to_data_dir/function.R")

# import
dat <- read.table(paste(data_dir, "R_Fig 5d.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$genotype <- str_replace(dat$genotype, " ", "")
dat$batch <- factor(dat$batch)

genotype <- c("wt", "bglu28bglu30")
dat$genotype <- factor(dat$genotype, levels=genotype)

sulfur <- c("S1500", "S150", "S15", "S0")
dat$sulfur <- factor(dat$sulfur, levels=sulfur)

# distribution
pdf(paste(data_dir, "data_nomral.pdf", sep=""))
	hist(dat$value, breaks=20)
	qqnorm(dat$value)
	qqline(dat$value, col='red')

	hist(log10(dat$value), breaks=20)
	qqnorm(log10(dat$value))
	qqline(log10(dat$value), col='red')
dev.off()

# fit
lmer_fit <- lmer(log10(value) ~ sulfur:genotype - 1 + (1|batch), dat)
dat_lmer <- summary(lmer_fit)

# diagnosis
pdf(paste(data_dir, "diagnosis.pdf", sep=""))
	plot(fitted(lmer_fit), resid(lmer_fit))
	abline(0, 0, col="red")

	qqnorm(dat_lmer$resid)
	qqline(dat_lmer$resid, col="red")
dev.off()

# pairwise test
adj_p.val <- pairwise_ttest(dat_lmer)
p.letters <- multcompLetters(adj_p.val)

names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "sulfur", "")
names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "genotype", "")

sort <- apply(expand.grid(genotype, sulfur)[, 2:1], 1, function(x) paste(x, collapse=":"))
write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, "FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)



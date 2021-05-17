# Sugiyama et al., PNAS, 2021
# Fig 5F
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

# cleanup
rm(list=ls())

# initialization
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/path_to_data_dir/R_Fig5f/"

# import
dat <- read.table(paste(data_dir, "R_Fig 5f.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$batch <- as.factor(dat$batch)

genotype <- c("wt", "bglu28 bglu30")
sulfur <- c("S200", "S100", "4MTB", "4MSB", "PhE")

idx <- dat$sulfur %in% sulfur
dat <- dat[idx,]

dat$genotype <- factor(dat$genotype, genotype)
dat$sulfur   <- factor(dat$sulfur, sulfur)

# distribution
pdf(paste(data_dir, "data_nomral.pdf", sep=""))
	hist(dat$value, breaks=20)
	qqnorm(dat$value)
	qqline(dat$value, col='red')

	hist(sqrt(dat$value), breaks=20)
	qqnorm(sqrt(dat$value))
	qqline(sqrt(dat$value), col='red')

	hist(log10(dat$value), breaks=20)
	qqnorm(log10(dat$value))
	qqline(log10(dat$value), col='red')
dev.off()

# fit
fit <- lm(value ~ 0 + genotype:sulfur + batch, dat)
anova <- anova(aov(fit))

sink(paste(data_dir, "lmer_summary.txt", sep=""))
	print(anova)
sink()

# diagnosis
pdf(paste(data_dir, "diagnosis.pdf", sep=""))
	plot(fitted(fit), resid(fit))
	abline(0, 0, col="red")

	qqnorm(summary(fit)$resid)
	qqline(summary(fit)$resid, col="red")
dev.off()

# pairwise test
tukey <- TukeyHSD(aov(fit))
p.letters <- multcompLetters2(value ~ genotype:sulfur, tukey$'genotype:sulfur'[, "p adj"], dat, reversed=TRUE)

sort <- expand.grid(genotype, sulfur) %>% apply(1, paste, collapse=":")
write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, "FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)


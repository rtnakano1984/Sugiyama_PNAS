# Sugiyama et al., PNAS, 2021
# Fig 5G,H,I
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

# cleanup
rm(list=ls())

# initialization
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/path_to_data_dir/R_Fig5ghi/"

# import
dat <- read.table(paste(data_dir, "R_Fig 5ghi.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$batch <- as.factor(dat$batch)

genotype <- c("wt", "bglu28bglu30")
sulfur <- c("S200", "S100", "S0", "4MTB", "4MSB", "PhE")
GL <- c("4MSB", "4MTB", "PhE")

dat$genotype <- factor(str_replace(dat$genotype, " ", ""), genotype)

colnames(dat)[5] <- "value"

# total, Cys
for(x in unique(dat$GL)){

	# total
	idx <- dat$GL == x
	temp <- dat[idx,]

	idx <- sulfur %in% temp$sulfur
	sulfur_temp <- sulfur[idx]

	temp$sulfur <- factor(temp$sulfur, levels=sulfur_temp)

	pdf(paste(data_dir, x, "-data_nomral.pdf", sep=""))
		hist(temp$value, breaks=20)
		qqnorm(temp$value)
		qqline(temp$value, col='red')

		hist(sqrt(temp$value), breaks=20)
		qqnorm(sqrt(temp$value))
		qqline(sqrt(temp$value), col='red')

		hist(log10(temp$value), breaks=20)
		qqnorm(log10(temp$value))
		qqline(log10(temp$value), col='red')
	dev.off()

	# fit
	if(x == "4MTB") {
		fit <- lm(sqrt(value) ~ 0 + genotype:sulfur + batch, temp)
	} else {
		fit <- lm(log10(value) ~ 0 + genotype:sulfur + batch, temp)
	}
	anova <- anova(aov(fit))

	sink(paste(data_dir, x, "-anova_summary.txt", sep=""))
		print(anova)
	sink()

	# diagnosis
	pdf(paste(data_dir, x, "-diagnosis.pdf", sep=""))
		plot(fitted(fit), resid(fit))
		abline(0, 0, col="red")

		qqnorm(summary(fit)$resid)
		qqline(summary(fit)$resid, col="red")
	dev.off()

	# pairwise test
	tukey <- TukeyHSD(aov(fit))
	p.letters <- multcompLetters2(value ~ genotype:sulfur, tukey$'genotype:sulfur'[, "p adj"], temp, reversed=TRUE)

	sort <- expand.grid(unique(temp$genotype),unique(temp$sulfur)) %>% apply(1, paste, collapse=":")
	write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, x, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}



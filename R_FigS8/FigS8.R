# Sugiyama et al., PNAS, 2021
# Fig S8
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

# cleanup
rm(list=ls())

# initialization
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/path_to_data_dir/R_FigS8/"

# import
dat <- read.table(paste(data_dir, "R_Fig S8.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$batch <- as.factor(dat$batch)

sulfur <- c("S150", "S0", "4MTB", "allyl", "PhE", "I3G")
compound <- c("4MTB-NH2", "allyl-NH2", "PhE-NH2", "RA")

dat$sulfur   <- factor(dat$sulfur, sulfur)

# total, Cys
for(x in compound){

	# total
	idx <- dat$compound == x
	temp <- dat[idx,]

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

	# anova fitting
	if(x == "RA") {
		fit <- lm(sqrt(value) ~ 0 + sulfur + batch, temp)
	} else {
		fit <- lm(log10(value) ~ 0 + sulfur + batch, temp)
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
	p.letters <- multcompLetters2(value ~ sulfur, tukey$'sulfur'[, "p adj"], temp, reversed=TRUE)

	sort <- sulfur
	write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, x, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}



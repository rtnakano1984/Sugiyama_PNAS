# Ryosuke Sugiyama et al., PNAS
# Fig S7a
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

rm(list=ls())

# library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig4de_new/"
# data_dir <- "~/Desktop/ryosuke_temp/R_Fig4de_new/"

dat <- read.table(paste(data_dir, "R_Fig 4de.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$batch <- as.factor(dat$batch)

genotype <- c("wt", "oxp1")
time <- c("00h", "09h", "24h", "48h")
treatment <- c("4MSB", "RA")

idx <- dat$treatment %in% treatment
dat <- dat[idx,]

dat$genotype <- factor(dat$genotype, genotype)
dat$time     <- factor(dat$time, time)


# total, Cys
for(x in unique(dat$treatment)){

	# total
	idx <- dat$treatment == x
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

	fit <- lm(sqrt(value) ~ 0 + genotype:time + batch, temp)
	anova <- anova(aov(fit))

	sink(paste(data_dir, x, "-lmer_summary.txt", sep=""))
	print(anova)
	sink()

	pdf(paste(data_dir, x, "-diagnosis.pdf", sep=""))
		plot(fitted(fit), resid(fit))
		abline(0, 0, col="red")

		qqnorm(summary(fit)$resid)
		qqline(summary(fit)$resid, col="red")
	dev.off()

	tukey <- TukeyHSD(aov(fit))
	idx <- order(-tukey$'genotype:time'[, 'diff'])
	p.letters <- multcompLetters(tukey$'genotype:time'[idx, 'p adj'])

	sort <- expand.grid(genotype, time) %>% apply(1, paste, collapse=":")
	write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, x, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}



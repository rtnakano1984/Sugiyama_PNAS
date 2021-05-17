# Sugiyama et al., PNAS, 2021
# Fig S2B
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

# cleanup
rm(list=ls())

# initialization
library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/path_to_data_dir/R_FigS2b_ratio/"
source("/path_to_data_dir/function.R")

# import
dat <- read.table(paste(data_dir, "R_Fig S2b_ratio.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$sulfur <- str_replace(dat$sulfur, "-", "_")
dat$batch <- as.factor(dat$batch)

sulfur <- c("S1500", "S150", "4MSB", "4MSB_34S")
dat$sulfur <- factor(dat$sulfur, levels=sulfur)

# for each compound
for(cpd in unique(dat$cpd)){

	# total
	idx <- dat$cpd == cpd
	temp <- dat[idx,]
	
	pdf(paste(data_dir, cpd, "-data_nomral.pdf", sep=""))
		hist(temp$value, breaks=20)
		qqnorm(temp$value)
		qqline(temp$value, col='red')

		hist(log10(temp$value), breaks=20)
		qqnorm(log10(temp$value))
		qqline(log10(temp$value), col='red')
	dev.off()

	lmer_fit <- lmer(log10(value) ~ sulfur - 1 + (1|batch), temp)
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
	write.table(as.data.frame(p.letters$Letters[sulfur]), file=paste(data_dir, cpd, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}



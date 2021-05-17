# Sugiyama et al., PNAS, 2021
# Fig S9B
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

# cleanup
rm(list=ls())

# initialization
library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/path_to_data_dir/R_FigS9b/"
source("/path_to_data_dir/function.R")

# import
dat_all <- read.table(paste(data_dir, "R_Fig S9b.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat_all$batch <- factor(dat_all$batch)
dat_all$raw <- dat_all$value

# for each gene
for(x in unique(dat_all$gene)){

	idx <- dat_all$gene == x
	dat <- dat_all[idx,]

	# distribution
	pdf(paste(data_dir, x, "-data_nomral.pdf", sep=""))
		hist(dat$value, breaks=20)
		qqnorm(dat$value)
		qqline(dat$value, col='red')

		hist(log10(dat$value), breaks=20)
		qqnorm(log10(dat$value[dat$value != 0]))
		qqline(log10(dat$value[dat$value != 0]), col='red')

		hist(sqrt(dat$value), breaks=20)
		qqnorm(sqrt(dat$value))
		qqline(sqrt(dat$value), col='red')
	dev.off()

	# fit
	lmer_fit <- lmer(sqrt(value) ~ cpd - 1 + (1|batch), dat)
	dat_lmer <- summary(lmer_fit)

	sink(paste(data_dir, x, "-lmer_summary.txt", sep=""))
		print(dat_lmer)
	sink()

	# diagnosis
	pdf(paste(data_dir, x, "-diagnosis.pdf", sep=""))
		plot(fitted(lmer_fit), resid(lmer_fit))
		abline(0, 0, col="red")

		qqnorm(dat_lmer$resid)
		qqline(dat_lmer$resid, col="red")
	dev.off()

	# pairwise test
	adj_p.val <- pairwise_ttest(dat_lmer)
	p.letters <- multcompLetters(adj_p.val)

	names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "cpd", "")

	sort <- c("Allyl", "3MSP", "4MSB", "4MSB(en)", "4MTB", "Bn", "PhE", "I3G", "1MI3G")
	write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, x, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}



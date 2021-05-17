# Sugiyama et al., PNAS, 2021
# Fig 4F,G
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

# cleanup
rm(list=ls())

# initialization
library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/path_to_data_dir/R_Fig4fg/"
source("/path_to_data_dir/function.R")

# import
dat <- read.table(paste(data_dir, "R_Fig 4fg.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$batch <- as.factor(dat$batch)

genotype <- c("wt", "oxp1")
sulfur <- c("S1500", "S15")
substrate <- c("pCys", "RA", "5OP")

dat$genotype <- factor(dat$genotype, genotype)
dat$sulfur   <- factor(dat$sulfur, sulfur)

# for each substrate
for(x in unique(dat$substrate)){

	# total
	idx <- dat$substrate == x
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

	if(x == "pCys"){
		lmer_fit <- lmer(log10(value) ~ genotype:sulfur - 1 + (1|batch), temp)		
	} else {
		lmer_fit <- lmer(sqrt(value) ~ genotype:sulfur - 1 + (1|batch), temp)
	}
	dat_lmer <- summary(lmer_fit)

	sink(paste(data_dir, x, "-lmer_summary.txt", sep=""))
	print(dat_lmer)
	sink()

	pdf(paste(data_dir, x, "-diagnosis.pdf", sep=""))
		plot(fitted(lmer_fit), resid(lmer_fit))
		abline(0, 0, col="red")

		qqnorm(dat_lmer$resid)
		qqline(dat_lmer$resid, col="red")
	dev.off()

	adj_p.val <- pairwise_ttest(dat_lmer)
	p.letters <- multcompLetters(adj_p.val)

	sort <- expand.grid(paste("genotype", genotype, sep=""), paste("sulfur", sulfur, sep="")) %>% apply(1, paste, collapse=":")
	write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, x, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}


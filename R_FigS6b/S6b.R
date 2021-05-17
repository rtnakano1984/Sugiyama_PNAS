# Sugiyama et al., PNAS, 2021
# Fig S6B
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

# cleanup
rm(list=ls())

# initialization
library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/path_to_data_dir/R_FigS6b/"
source("/path_to_data_dir/function.R")

# import
dat_all <- read.table(paste(data_dir, "R_Fig S6b.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat_all$treatment <- str_replace(dat_all$treatment, "-", "_")
dat_all$batch <- as.factor(dat_all$batch)


idx <- dat_all$treatment == "BSO_4MSB"
dat_all$genotype[idx] <- "wt_BSO"

genotype <- c("wt", "wt_BSO", "pad2")
dat_all$genotype <- factor(dat_all$genotype, levels=genotype)

time <- c("00h", "09h", "24h", "48h")
dat_all$time <- factor(dat_all$time, levels=time)

for(cpd in unique(dat_all$cpd)){
	idx <- dat_all$cpd == cpd
	temp <- dat_all[idx,]

	pdf(paste(data_dir, cpd, "-data_nomral.pdf", sep=""))
		hist(temp$value, breaks=20)
		qqnorm(temp$value)
		qqline(temp$value, col='red')

		hist(log10(temp$value), breaks=20)
		qqnorm(log10(temp$value))
		qqline(log10(temp$value), col='red')

		hist(sqrt(temp$value), breaks=20)
		qqnorm(sqrt(temp$value))
		qqline(sqrt(temp$value), col='red')
	dev.off()

	if(cpd %in% c("4MSB", "GSH", "SFN-GSH", "SFN-rEC")){
		lmer_fit <- lmer(sqrt(value) ~ time:genotype - 1 + (1|batch), temp)
	} else if(cpd %in% c("SFN-Cys", "4MSB")) {
		lmer_fit <- lmer(log10(value) ~ time:genotype - 1 + (1|batch), temp)		
	} else {
		lmer_fit <- lmer(value ~ time:genotype - 1 + (1|batch), temp)		
	}
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

	names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "time", "")
	names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "genotype", "")

	sort <- apply(expand.grid(genotype, time)[2:1], 1, function(x) paste(x, collapse=":"))
	write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, cpd, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}



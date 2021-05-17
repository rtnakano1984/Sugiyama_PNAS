# Sugiyama et al., PNAS, 2021
# Fig 5A,B
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

# cleanup
rm(list=ls())

# initialization
library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/path_to_data_dir/R_Fig5ab/"
source("/path_to_data_dir/function.R")

# import
dat_all <- read.table(paste(data_dir, "R_Fig 5ab.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat_all$batch <- factor(dat_all$batch)
dat_all$raw <- dat_all$value

# relative activity, taking the maximum mean as 1
norm_fac <- dat_all %>%
	group_by(gene, cpd) %>%
	summarise(mean=mean(value)) %>%
	group_by(gene) %>%
	summarise(norm_fac=max(mean)) %>%
	data.frame(., stringsAsFactors=F)
for(x in unique(dat_all$gene)) {
	dat_all$rel[dat_all$gene==x] <- dat_all$value[dat_all$gene==x] / norm_fac$norm_fac[norm_fac$gene==x]
}

# for each gene
for(x in unique(dat_all$gene)){

	idx <- dat_all$gene == x
	dat <- dat_all[idx,]

	# distribution
	pdf(paste(data_dir, x, "-data_nomral-rel.pdf", sep=""))
		hist(dat$rel, breaks=20)
		qqnorm(dat$rel)
		qqline(dat$rel, col='red')

		hist(log10(dat$rel), breaks=20)
		qqnorm(log10(dat$rel[dat$rel != 0]))
		qqline(log10(dat$rel[dat$rel != 0]), col='red')

		hist(sqrt(dat$rel), breaks=20)
		qqnorm(sqrt(dat$rel))
		qqline(sqrt(dat$rel), col='red')
	dev.off()

	# fit
	lmer_fit <- lmer(sqrt(rel) ~ cpd - 1 + (1|batch), dat)
	dat_lmer <- summary(lmer_fit)

	sink(paste(data_dir, x, "-lmer_summary-rel.txt", sep=""))
		print(dat_lmer)
	sink()

	# diagnosis
	pdf(paste(data_dir, x, "-diagnosis-rel.pdf", sep=""))
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
	write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, x, "-FDR_letters-rel.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}



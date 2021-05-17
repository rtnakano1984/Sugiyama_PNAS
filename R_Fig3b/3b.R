# Sugiyama et al., PNAS, 2021
# Fig 3B
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

# cleanup
rm(list=ls())

# initialization
library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/path_to_data_dir/R_Fig3b/"
source("/path_to_data_dir/function.R")

# import
dat <- read.table(paste(data_dir, "R_Fig 3b.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$treatment <- str_replace(dat$treatment, "-", "_")
dat$batch <- as.factor(dat$batch)

treatment <- c("water", "4MSB", "4MSB_d5")
dat$treatment <- factor(dat$treatment, levels=treatment)

time <- c("03h", "09h", "24h", "48h")
dat$time <- factor(dat$time, levels=time)

# for each compound
for(cpd in unique(dat$cpd)){

	# total
	idx <- dat$cpd == cpd
	temp <- dat[idx,] %>% group_by(time, treatment, batch) %>% summarize(value=sum(value)) %>% data.frame(., stringsAsFactors=F)

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

	if(cpd == "SFN-rEC"){
		lmer_fit <- lmer(log10(value) ~ treatment:time - 1 + (1|batch), temp)
	} else if(cpd %in% c("SFN-Cys", "GSH", "RA")){
		lmer_fit <- lmer(sqrt(value) ~ treatment:time - 1 + (1|batch), temp)
	} else {
		lmer_fit <- lmer(value ~ treatment:time - 1 + (1|batch), temp)
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
	names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "treatment", "")
	sort <- apply(expand.grid(time, treatment)[, 2:1], 1, function(x) paste(x, collapse=":"))

	p.letters <- p.letters$Letters[sort]

	write.table(as.data.frame(p.letters), file=paste(data_dir, cpd, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)


	if(!cpd %in% c("GSH", "RA")){
		# isotope 0
		idx <- dat$cpd == cpd & dat$isotope == 0
		temp <- dat[idx,]

		pdf(paste(data_dir, cpd, "-isotope_0-data_nomral.pdf", sep=""))
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

		if(cpd %in% c("4MSB", "4MSB-NH2", "SFN-GSH")){
			lmer_fit <- lmer(log10(value) ~ treatment:time - 1 + (1|batch), temp)
		} else {
			lmer_fit <- lmer(value ~ treatment:time - 1 + (1|batch), temp)
		}
		dat_lmer <- summary(lmer_fit)

		sink(paste(data_dir, cpd, "-isotope_0-lmer_summary.txt", sep=""))
		print(dat_lmer)
		sink()
		
		pdf(paste(data_dir, cpd, "-isotope_0-diagnosis.pdf", sep=""))
			plot(fitted(lmer_fit), resid(lmer_fit))
			abline(0, 0, col="red")

			qqnorm(dat_lmer$resid)
			qqline(dat_lmer$resid, col="red")
		dev.off()


		adj_p.val <- pairwise_ttest(dat_lmer)
		p.letters <- multcompLetters(adj_p.val)
		
		names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "time", "")
		names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "treatment", "")
		sort <- apply(expand.grid(time, treatment)[, 2:1], 1, function(x) paste(x, collapse=":"))

		p.letters <- p.letters$Letters[sort]
		write.table(as.data.frame(p.letters), file=paste(data_dir, cpd, "-isotope_0-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

		# isotope 5
		idx <- dat$cpd == cpd & dat$isotope == 5
		temp <- dat[idx,]

		pdf(paste(data_dir, cpd, "-isotope_5-data_nomral.pdf", sep=""))
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

		if(cpd %in% c("4MSB-NH2", "SFN-GSH")){
			lmer_fit <- lmer(log10(value) ~ treatment:time - 1 + (1|batch), temp)
		} else {
			lmer_fit <- lmer(value ~ treatment:time - 1 + (1|batch), temp)
		}
		dat_lmer <- summary(lmer_fit)

		sink(paste(data_dir, cpd, "-isotope_5-lmer_summary.txt", sep=""))
		print(dat_lmer)
		sink()

		pdf(paste(data_dir, cpd, "-isotope_5-diagnosis.pdf", sep=""))
			plot(fitted(lmer_fit), resid(lmer_fit))
			abline(0, 0, col="red")

			qqnorm(dat_lmer$resid)
			qqline(dat_lmer$resid, col="red")
		dev.off()


		adj_p.val <- pairwise_ttest(dat_lmer)
		p.letters <- multcompLetters(adj_p.val)
		
		names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "time", "")
		names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "treatment", "")
		sort <- apply(expand.grid(time, treatment)[, 2:1], 1, function(x) paste(x, collapse=":"))

		p.letters <- p.letters$Letters[sort]
		write.table(as.data.frame(p.letters), file=paste(data_dir, cpd, "-isotope_5-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)
	}
}



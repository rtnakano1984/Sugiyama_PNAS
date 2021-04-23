# Ryosuke Sugiyama et al., PNAS
# Fig 3b
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

rm(list=ls())

library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig 3b/"

dat <- read.table(paste(data_dir, "R_Fig 3b.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$treatment <- str_replace(dat$treatment, "-", "_")
dat$batch <- as.factor(dat$batch)

treatment <- c("water", "4MSB", "4MSB_d5")
dat$treatment <- factor(dat$treatment, levels=treatment)

time <- c("03h", "09h", "24h", "48h")
dat$time <- factor(dat$time, levels=time)


pairwise_test <- function(dat_lmer){

	estim <- dat_lmer$coefficients[,1]
	df <- length(dat_lmer$residuals) - length(estim) -  (sum(as.vector(dat_lmer$ngrps)) - 1 )
	vcov <- as.matrix(dat_lmer$vcov)

	idx <- order(estim)
	estim <- estim[idx]
	vcov <- vcov[idx,idx]


	group <- names(estim)

	group_comp <- c()    						
	for (i in 1:(length(group)-1)){							
	  for (j in (i+1):length(group)){
	    group_comp <- c(group_comp, paste(group[i], group[j], sep="-"))
	  }
	}

	id.mat <- matrix(0, ncol=1, nrow = length(group) )
	p.val <- c()
	for ( i in 1:(length(estim)-1)) {
	  for (j in (i+1):length(estim)){
	    id.mat.x <- id.mat
	    id.mat.x[ i, 1 ] <- 1
	    id.mat.x[ j, 1 ] <- -1
	    stder <- sqrt( t(id.mat.x) %*% vcov %*% id.mat.x )
	    t.val <- abs( estim[i]-estim[j]) / stder
	    p.val <- c( p.val, 2 * pt( t.val, df, lower.tail=F ) )
	  }
	}

	names(p.val) <- group_comp

	adj_p.val <- p.adjust(p.val, "fdr")
	p.letters <- multcompLetters(adj_p.val)
	
	names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "time", "")
	names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "treatment", "")
	sort <- apply(expand.grid(time, treatment)[, 2:1], 1, function(x) paste(x, collapse=":"))

	return(p.letters$Letters[sort])
}

# 
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

	p.letters <- pairwise_test(dat_lmer)

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

		p.letters <- pairwise_test(dat_lmer)

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

		p.letters <- pairwise_test(dat_lmer)

		write.table(as.data.frame(p.letters), file=paste(data_dir, cpd, "-isotope_5-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)
	}
}



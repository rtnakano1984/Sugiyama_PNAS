# Ryosuke Sugiyama et al., PNAS
# Fig S6b
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

rm(list=ls())

library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig S6b/"

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
	names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "genotype", "")

	sort <- apply(expand.grid(genotype, time)[2:1], 1, function(x) paste(x, collapse=":"))
	write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, cpd, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}



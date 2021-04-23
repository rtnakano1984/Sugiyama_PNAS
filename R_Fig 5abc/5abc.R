# Ryosuke Sugiyama et al., PNAS
# Fig 5abc
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

rm(list=ls())

library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig 5abc/"

dat_all <- read.table(paste(data_dir, "R_Fig 5abc.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
# dat$genotype <- str_replace(dat$genotype, " ", "")
dat_all$batch <- factor(dat_all$batch)
dat_all$raw <- dat_all$value

# rel
norm_fac <- dat_all %>% group_by(gene, cpd) %>% summarise(mean=mean(value)) %>% group_by(gene) %>% summarise(norm_fac=max(mean)) %>% data.frame(., stringsAsFactors=F)
for(x in unique(dat_all$gene)) dat_all$rel[dat_all$gene==x] <- dat_all$value[dat_all$gene==x] / norm_fac$norm_fac[norm_fac$gene==x]

# 
pdf(paste(data_dir, "data_nomral.pdf", sep=""))
	hist(dat_all$value, breaks=20)
	qqnorm(dat_all$value)
	qqline(dat_all$value, col='red')

	hist(log10(dat_all$value), breaks=20)
	qqnorm(log10(dat_all$value[dat_all$value != 0]))
	qqline(log10(dat_all$value[dat_all$value != 0]), col='red')
dev.off()


for(x in unique(dat_all$gene)){

	idx <- dat_all$gene == x
	dat <- dat_all[idx,]

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

	lmer_fit <- lmer(sqrt(rel) ~ cpd - 1 + (1|batch), dat)
	dat_lmer <- summary(lmer_fit)

	sink(paste(data_dir, x, "-lmer_summary-rel.txt", sep=""))
	print(dat_lmer)
	sink()

	pdf(paste(data_dir, x, "-diagnosis-rel.pdf", sep=""))
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

	names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "cpd", "")

	sort <- c("SIN", "GIB", "GRA", "GRE", "GEC", "Bn", "PhE", "I3G", "1MI3G")
	write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, x, "-FDR_letters-rel.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)



	idx <- dat_all$gene == x
	dat <- dat_all[idx,]

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

	# lmer_fit <- lmer(log10(value) ~ sulfur - 1 + (1|batch), dat)
	lmer_fit <- lmer(sqrt(value) ~ cpd - 1 + (1|batch), dat)
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

	names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "cpd", "")

	write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, x, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}



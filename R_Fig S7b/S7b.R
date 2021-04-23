# Ryosuke Sugiyama et al., PNAS
# Fig S7b
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

rm(list=ls())

library(lme4)
library(dplyr)
library(stringr)
library(multcompView)

data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig S7b/"

dat <- read.table(paste(data_dir, "R_Fig S7b.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$sulfur <- str_replace(dat$sulfur, "-", "_")
dat$batch <- as.factor(dat$batch)

genotype <- c("wt", "qko")
dat$genotype <- factor(dat$genotype, levels=genotype)

sulfur <- c("S1500", "S150")
dat$sulfur <- factor(dat$sulfur, levels=sulfur)


# total
pdf(paste(data_dir, "data_nomral.pdf", sep=""))
	hist(dat$value, breaks=20)
	qqnorm(dat$value)
	qqline(dat$value, col='red')

	hist(log10(dat$value), breaks=20)
	qqnorm(log10(dat$value))
	qqline(log10(dat$value), col='red')
dev.off()

total <- dat %>% group_by(sulfur, genotype, batch) %>% summarise(value=sum(value))
lmer_fit <- lmer(log10(value) ~ sulfur:genotype - 1 + (1|batch), total)
dat_lmer <- summary(lmer_fit)

sink(paste(data_dir, "lmer_summary.txt", sep=""))
print(dat_lmer)
sink()

pdf(paste(data_dir, "diagnosis.pdf", sep=""))
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

names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "sulfur", "")
names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "genotype", "")
sort <- apply(expand.grid(genotype, sulfur)[2:1], 1, function(x) paste(x, collapse=":"))

write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, "FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)




# isotopes
for(iso in unique(dat$isotope)){
	idx <- dat$isotope == iso
	temp <- dat[idx,]

	pdf(paste(data_dir, iso, "-data_nomral.pdf", sep=""))
		hist(temp$value, breaks=20)
		qqnorm(temp$value)
		qqline(temp$value, col='red')

		hist(log10(temp$value), breaks=20)
		qqnorm(log10(temp$value))
		qqline(log10(temp$value), col='red')
	dev.off()


	lmer_fit <- lmer(log10(value) ~ sulfur:genotype - 1 + (1|batch), temp)
	dat_lmer <- summary(lmer_fit)

	sink(paste(data_dir, iso, "-lmer_summary.txt", sep=""))
	print(dat_lmer)
	sink()

	pdf(paste(data_dir, iso, "-diagnosis.pdf", sep=""))
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

	names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "sulfur", "")
	names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "genotype", "")
	sort <- apply(expand.grid(genotype, sulfur)[2:1], 1, function(x) paste(x, collapse=":"))

	write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, iso, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)
}








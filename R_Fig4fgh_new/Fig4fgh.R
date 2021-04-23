# Ryosuke Sugiyama et al., PNAS
# Fig S7a
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

rm(list=ls())

library(lme4)
library(dplyr)
library(stringr)
library(multcompView)
library(ggplot2)

data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig4fgh_new/"
# data_dir <- "~/Desktop/ryosuke_temp/R_Fig4fgh_new/"

dat <- read.table(paste(data_dir, "R_Fig 4fgh.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat$batch <- as.factor(dat$batch)

genotype <- c("wt", "oxp1")
sulfur <- c("S1500", "S15")
substrate <- c("pCys", "RA", "5OP")


dat$genotype <- factor(dat$genotype, genotype)
dat$sulfur   <- factor(dat$sulfur, sulfur)

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

	return(p.letters)
}


# total, Cys
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

	p.letters <- pairwise_test(dat_lmer)

	sort <- expand.grid(paste("genotype", genotype, sep=""), paste("sulfur", sulfur, sep="")) %>% apply(1, paste, collapse=":")
	write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, x, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}

p <- ggplot(dat, aes(x=sulfur, y=value, fill=genotype)) +
	geom_boxplot(outlier.shape=NA) +
	geom_point(size=.5, alpha=.75, position=position_jitterdodge(jitter.width=.25)) +
	scale_y_sqrt() +
	facet_grid(. ~ substrate) +
	theme_bw()
ggsave(p, file=paste(data_dir, "plot.pdf", sep=""), width=5, height=4, bg="transparent")


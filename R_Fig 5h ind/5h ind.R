# Ryosuke Sugiyama et al., PNAS
# Fig S10b
# By Ryohei Thomas Nakano, PhD; nakano@mpipz.mpg.de

rm(list=ls())

library(lme4)
library(dplyr)
library(stringr)
library(multcompView)
library(ggplot2)
library(patchwork)
source("/biodata/dep_psl/grp_psl/ThomasN/scripts/ggplot-themes_RTN.R")


data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig 5h ind/"

dat_all <- read.table(paste(data_dir, "R_Fig 5h.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat_all$batch <- factor(dat_all$batch)

genotype <- c("wt", "bglu28bglu30")
dat_all$genotype <- factor(dat_all$genotype, levels=genotype)

sulfur <- c("S1500", "S30", "S3")
dat_all$sulfur <- factor(dat_all$sulfur, levels=sulfur)

# empty df for estimates
estim_df <- data.frame(group=NULL, estim=NULL, cpd=NULL, batch=NULL)

# individual glucosinolate species
for(cpd in unique(dat_all$cpd)){
  for(batch in unique(dat_all$batch)){
    idx <- dat_all$cpd == cpd & dat_all$batch == batch
    dat <- dat_all[idx,]

    pdf(paste(data_dir, cpd, "-", batch, "-data_nomral.pdf", sep=""))
      hist(dat$value, breaks=20)
      qqnorm(dat$value)
      qqline(dat$value, col='red')

      hist(log10(dat$value), breaks=20)
      qqnorm(log10(dat$value))
      qqline(log10(dat$value), col='red')

      hist(sqrt(dat$value), breaks=20)
      qqnorm(sqrt(dat$value))
      qqline(sqrt(dat$value), col='red')
    dev.off()

    lmer_fit <- lmer(sqrt(value) ~ sulfur:genotype - 1 + sample + (1|pot), dat)

    dat_lmer <- summary(lmer_fit)

    sink(paste(data_dir, cpd, "-", batch, "-lmer_summary.txt", sep=""))
    print(dat_lmer)
    sink()

    pdf(paste(data_dir, cpd, "-", batch, "-diagnosis.pdf", sep=""))
      plot(fitted(lmer_fit), resid(lmer_fit))
      abline(0, 0, col="red")

      qqnorm(dat_lmer$resid)
      qqline(dat_lmer$resid, col="red")
    dev.off()


    estim <- dat_lmer$coefficients[,1]
    df <- length(dat_lmer$residuals) - length(estim) -  (sum(as.vector(dat_lmer$ngrps)) - 1 )
    vcov <- as.matrix(dat_lmer$vcov)

    idx <- names(estim) != "sample"
    estim <- estim[idx]

    idx <- order(estim)
    estim <- estim[idx]

    idx <- match(names(estim), rownames(vcov))
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

    write.table(data.frame(group_comp=group_comp, p.val=p.val, adj_p=adj_p.val), paste(data_dir, cpd, "-", batch, "-P_values.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")

    p.letters <- multcompLetters(adj_p.val)

    names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "sulfur", "")
    names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "genotype", "")

    sort <- apply(expand.grid(genotype, sulfur)[, 2:1], 1, function(x) paste(x, collapse=":"))
    p.letters_df <- as.data.frame(p.letters$Letters[sort])
    write.table(p.letters_df, file=paste(data_dir, cpd, "-", batch, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)


    names(estim) <- str_replace_all(names(estim), "sulfur", "")
    names(estim) <- str_replace_all(names(estim), "genotype", "")
    idx <- match(sort, names(estim))
    estim <- estim[idx]
    estim_df_temp <- data.frame(group=names(estim), estim=estim^2, cpd=cpd, batch=batch, row.names=NULL, letters=p.letters_df[,1], stringsAsFactors=F)
    estim_df_temp$logFC <- log2(estim_df_temp$estim / estim_df_temp$estim[estim_df_temp$group == "S1500:wt"])

    estim_df <- rbind(estim_df, estim_df_temp)

  }
}

write.table(estim_df, file=paste(data_dir, "logFC.txt", sep=""), sep="\t", row.names=F, col.names=T, quote=F)


# boxplot
genotype <- c("wt", "bglu28bglu30")
dat_all$genotype <- factor(dat_all$genotype, levels=genotype)

sulfur <- c("S1500", "S30", "S3")
dat_all$sulfur <- factor(dat_all$sulfur, levels=sulfur)

p <- ggplot(dat_all, aes(x=genotype, y=value, group=genotype, shape=factor(pot), color=factor(batch))) +
  geom_boxplot(outlier.shape=NA, position=position_dodge()) +
  geom_point(, position=position_jitterdodge()) +
  facet_grid(cpd ~ batch + sulfur, scales="free") +
  theme_RTN +
  theme(axis.text.x=element_text(angle=75, hjust=1, vjust=1))
ggsave(p, file=paste(data_dir, "boxplot.pdf", sep=""), width=5, height=25, bg="transparent")



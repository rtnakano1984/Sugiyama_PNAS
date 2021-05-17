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


data_dir <- "/netscratch/dep_psl/grp_psl/ThomasN/ryosuke_temp/R_Fig 5h test/"

dat_all <- read.table(paste(data_dir, "R_Fig 5h.csv", sep=""), sep=",", header=T, stringsAsFactors=F)
dat_all$batch <- factor(dat_all$batch)

genotype <- c("wt", "bglu28bglu30")
dat_all$genotype <- factor(dat_all$genotype, levels=genotype)

sulfur <- c("S1500", "S30", "S3")
dat_all$sulfur <- factor(dat_all$sulfur, levels=sulfur)

# empty df
estim_df <- data.frame(group=NULL, estim=NULL, batch=NULL)

# individual glucosinolate species
idx <- dat_all$cpd == "G01"
dat_all <- dat_all[idx,]

pdf(paste(data_dir, "data_nomral.pdf", sep=""))
  hist(dat_all$value, breaks=20)
  qqnorm(dat_all$value)
  qqline(dat_all$value, col='red')

  hist(log10(dat_all$value), breaks=20)
  qqnorm(log10(dat_all$value))
  qqline(log10(dat_all$value), col='red')

  hist(sqrt(dat_all$value), breaks=20)
  qqnorm(sqrt(dat_all$value))
  qqline(sqrt(dat_all$value), col='red')
dev.off()

lmer_fit <- lmer(sqrt(value) ~ sulfur:genotype - 1 + sample + (1|pot) + (1|batch), dat_all)
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

idx <- names(estim) != "sample"
estim <- estim[idx]

estim_df_temp <- data.frame(group=names(estim), estim=estim^2, batch="ALL", row.names=NULL)
estim_df <- rbind(estim_df, estim_df_temp)

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

write.table(data.frame(group_comp=group_comp, p.val=p.val, adj_p=adj_p.val), paste(data_dir, "P_values.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")

p.letters <- multcompLetters(adj_p.val)

names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "sulfur", "")
names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "genotype", "")


sort <- apply(expand.grid(genotype, sulfur)[, 2:1], 1, function(x) paste(x, collapse=":"))

write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, "FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)





# batch by batch

  # individual glucosinolate species
for(x in unique(dat_all$batch)){
  idx <- dat_all$batch == x
  dat <- dat_all[idx,]

  pdf(paste(data_dir, x, "-data_nomral.pdf", sep=""))
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

  idx <- names(estim) != "sample"
  estim <- estim[idx]

  estim_df_temp <- data.frame(group=names(estim), estim=estim^2, batch=x, row.names=NULL)
  estim_df <- rbind(estim_df, estim_df_temp)

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

  write.table(data.frame(group_comp=group_comp, p.val=p.val, adj_p=adj_p.val), paste(data_dir, x, "-P_values.txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")

  p.letters <- multcompLetters(adj_p.val)

  names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "sulfur", "")
  names(p.letters$Letters) <- str_replace(names(p.letters$Letters), "genotype", "")


  sort <- apply(expand.grid(genotype, sulfur)[, 2:1], 1, function(x) paste(x, collapse=":"))

  write.table(as.data.frame(p.letters$Letters[sort]), file=paste(data_dir, x, "-FDR_letters.txt", sep=""), sep="\t", row.names=T,col.names=NA, quote=F)

}



# fold changes
estim_df$group <- str_replace_all(estim_df$group, "sulfur", "")
estim_df$group <- str_replace_all(estim_df$group, "genotype", "")
for(x in unique(estim_df$batch)){
  idx <- estim_df$batch == x
  estim_df$logFC[idx] <- log2(estim_df$estim[idx] / estim_df$estim[idx & estim_df$group == "S1500:wt"])
}

idx <- estim_df$batch == "ALL"
mean_FC_estim_all <- estim_df[idx,c("group", "logFC")]
mean_FC_estim_all <- mean_FC_estim_all[order(mean_FC_estim_all$group),]
mean_FC_estim_ind <- estim_df[!idx,] %>% group_by(group) %>% summarise(logFC=mean(logFC)) %>% data.frame(., stringsAsFactors=F)

mean_FC_val_all   <- dat_all %>% mutate(group=paste(sulfur, genotype, sep=":")) %>% group_by(group) %>% summarise(logFC=log2(mean(relative)/100)) %>% arrange(group) %>% data.frame(., stringsAsFactors=F)
mean_FC_val_ind   <- dat_all %>% mutate(group=paste(sulfur, genotype, sep=":")) %>% group_by(group, batch) %>% summarise(mean=mean(relative)) %>% group_by(group) %>% summarise(logFC=log2(mean(mean)/100)) %>% data.frame(., stringsAsFactors=F)

pdf(paste(data_dir, "logFC_compare.pdf", sep=""))
  plot(mean_FC_val_ind$logFC ~ mean_FC_estim_ind$logFC); abline(0, 1, col="red")
  plot(mean_FC_val_all$logFC ~ mean_FC_estim_all$logFC); abline(0, 1, col="red")
  plot(mean_FC_val_ind$logFC ~ mean_FC_val_all$logFC); abline(0, 1, col="red")
  plot(mean_FC_estim_ind$logFC ~ mean_FC_estim_all$logFC); abline(0, 1, col="red")
dev.off()




###############################
### Figures in UKB
###############################

library(data.table)
setwd("C:/Users/Admin/Documents/GeneticNurtureNonCog")

siblings <- read.table("UKB/Siblings/summary_mean_CI_siblings_UKB_pop_lm_20210519_BOTH.csv")
rownames(siblings) <- siblings$Estimates
siblings$Estimates <- NULL
siblings <- t(siblings)
siblings <- as.data.frame(siblings)
siblings$Measure <- row.names(siblings)
siblings$Methods <- "Siblings_UKB"
siblings$pheno <- "EA"
siblings$Zscore <- siblings$original/siblings$se
siblings$pvalue <- 2*pnorm(-abs(siblings$Zscore))

adoption <- read.table("UKB/Adoptees/summary_mean_CI_adoption_UKB_20200529.csv")
rownames(adoption) <- adoption$Estimates
adoption$Estimates <- NULL
adoption <- t(adoption)
adoption <- as.data.frame(adoption)
adoption$Measure <- row.names(adoption)
adoption$Methods <- "Adoption_UKB"
adoption$pheno <- "EA"
adoption$Zscore <- adoption$original/adoption$se
adoption$pvalue <- 2*pnorm(-abs(adoption$Zscore))

ntrsib <- read.table("NTR/summary_mean_CI_siblings_NTR_EA_pop_lm_20210519_BOTH.csv")
rownames(ntrsib) <- ntrsib$Estimates
ntrsib$Estimates <- NULL
ntrsib <- t(ntrsib)
ntrsib <- as.data.frame(ntrsib)
ntrsib$Measure <- row.names(ntrsib)
ntrsib$Methods <- "Siblings_NTR"
ntrsib$pheno <- "EA"
ntrsib$Zscore <- ntrsib$original/ntrsib$se
ntrsib$pvalue <- 2*pnorm(-abs(ntrsib$Zscore))

ntrsibcito <- read.table("NTR/summary_mean_CI_siblings_NTR_CITO_both_lm_20210519.csv")
rownames(ntrsibcito) <- ntrsibcito$Estimates
ntrsibcito$Estimates <- NULL
ntrsibcito <- t(ntrsibcito)
ntrsibcito <- as.data.frame(ntrsibcito)
ntrsibcito$Measure <- row.names(ntrsibcito)
ntrsibcito$Methods <- "Siblings_NTR"
ntrsibcito$pheno <- "CITO"
ntrsibcito$Zscore <- ntrsibcito$original/ntrsibcito$se
ntrsibcito$pvalue <- 2*pnorm(-abs(ntrsibcito$Zscore))

ntrtrioea <- read.table("NTR/summary_mean_CI_trios_NTR_EA_lm_20210518.csv")
rownames(ntrtrioea) <- ntrtrioea$Estimates
ntrtrioea$Estimates <- NULL
ntrtrioea <- t(ntrtrioea)
ntrtrioea <- as.data.frame(ntrtrioea)
ntrtrioea$Measure <- row.names(ntrtrioea)
ntrtrioea$Methods <- "Trios_NTR"
ntrtrioea$pheno <- "EA"
ntrtrioea$Zscore <- ntrtrioea$original/ntrtrioea$se
ntrtrioea$pvalue <- 2*pnorm(-abs(ntrtrioea$Zscore))

ntrtriocito <- read.table("NTR/summary_mean_CI_trios_NTR_CITO_lm_20210518.csv")
rownames(ntrtriocito) <- ntrtriocito$Estimates
ntrtriocito$Estimates <- NULL
ntrtriocito <- t(ntrtriocito)
ntrtriocito <- as.data.frame(ntrtriocito)
ntrtriocito$Measure <- row.names(ntrtriocito)
ntrtriocito$Methods <- "Trios_NTR"
ntrtriocito$pheno <- "CITO"
ntrtriocito$Zscore <- ntrtriocito$original/ntrtriocito$se
ntrtriocito$pvalue <- 2*pnorm(-abs(ntrtriocito$Zscore))

tedstwin12 <- read.table("TEDS/summary_mean_CI_siblings_TEDS_12_pop_lm_20210519_BOTH.csv")
rownames(tedstwin12) <- tedstwin12$Estimates
tedstwin12$Estimates <- NULL
tedstwin12 <- t(tedstwin12)
tedstwin12 <- as.data.frame(tedstwin12)
tedstwin12$Measure <- row.names(tedstwin12)
tedstwin12$Methods <- "Siblings_TEDS"
tedstwin12$pheno <- "12yo"
tedstwin12$Zscore <- tedstwin12$original/tedstwin12$se
tedstwin12$pvalue <- 2*pnorm(-abs(tedstwin12$Zscore))

tedstwinGCSE <- read.table("TEDS/summary_mean_CI_siblings_TEDS_16_pop_lm_20210519_BOTH.csv")
rownames(tedstwinGCSE) <- tedstwinGCSE$Estimates
tedstwinGCSE$Estimates <- NULL
tedstwinGCSE <- t(tedstwinGCSE)
tedstwinGCSE <- as.data.frame(tedstwinGCSE)
tedstwinGCSE$Measure <- row.names(tedstwinGCSE)
tedstwinGCSE$Methods <- "Siblings_TEDS"
tedstwinGCSE$pheno <- "GCSE"
tedstwinGCSE$Zscore <- tedstwinGCSE$original/tedstwinGCSE$se
tedstwinGCSE$pvalue <- 2*pnorm(-abs(tedstwinGCSE$Zscore))



# Rename sibling design measure to fit other names
siblings$Measure[siblings$Measure == "pop_NonCog"] <- "total_NonCog"
siblings$Measure[siblings$Measure == "pop_Cog"] <- "total_Cog"
siblings$Measure[siblings$Measure == "ratio_pop_NonCog"] <- "ratio_tot_NonCog"
siblings$Measure[siblings$Measure == "ratio_pop_Cog"] <- "ratio_tot_Cog"

ntrsib$Measure[ntrsib$Measure == "pop_NonCog"] <- "total_NonCog"
ntrsib$Measure[ntrsib$Measure == "pop_Cog"] <- "total_Cog"
ntrsib$Measure[ntrsib$Measure == "ratio_pop_NonCog"] <- "ratio_tot_NonCog"
ntrsib$Measure[ntrsib$Measure == "ratio_pop_Cog"] <- "ratio_tot_Cog"

ntrsibcito$Measure[ntrsibcito$Measure == "pop_NonCog"] <- "total_NonCog"
ntrsibcito$Measure[ntrsibcito$Measure == "pop_Cog"] <- "total_Cog"
ntrsibcito$Measure[ntrsibcito$Measure == "ratio_pop_NonCog"] <- "ratio_tot_NonCog"
ntrsibcito$Measure[ntrsibcito$Measure == "ratio_pop_Cog"] <- "ratio_tot_Cog"

tedstwin12$Measure[tedstwin12$Measure == "pop_NonCog"] <- "total_NonCog"
tedstwin12$Measure[tedstwin12$Measure == "pop_Cog"] <- "total_Cog"
tedstwin12$Measure[tedstwin12$Measure == "ratio_pop_NonCog"] <- "ratio_tot_NonCog"
tedstwin12$Measure[tedstwin12$Measure == "ratio_pop_Cog"] <- "ratio_tot_Cog"

tedstwinGCSE$Measure[tedstwinGCSE$Measure == "pop_NonCog"] <- "total_NonCog"
tedstwinGCSE$Measure[tedstwinGCSE$Measure == "pop_Cog"] <- "total_Cog"
tedstwinGCSE$Measure[tedstwinGCSE$Measure == "ratio_pop_NonCog"] <- "ratio_tot_NonCog"
tedstwinGCSE$Measure[tedstwinGCSE$Measure == "ratio_pop_Cog"] <- "ratio_tot_Cog"

data <- rbind(siblings, adoption, ntrsib, ntrsibcito, ntrtrioea, ntrtriocito, tedstwin12, tedstwinGCSE)
head(data)
library(stringr)
# Limit to the Measures we want 
datafig <- data[which(data$Measure == "direct_Cog" |
                      data$Measure == "direct_NonCog" |
                      data$Measure == "indirect_Cog" |
                      data$Measure == "indirect_NonCog" |
                      data$Measure == "total_Cog" |
                      data$Measure == "total_NonCog" |
                      data$Measure == "ratio_tot_Cog" |
                      data$Measure == "ratio_tot_NonCog"),]
data <- datafig

#Split Measure into Type and PRS and change format back to dataframe
split <- unlist(strsplit(data$Measure, "_(?=[^_]+$)", perl=TRUE))
cols <- c("Type", "PRS")
nC <- length(cols)
ind <- seq(from=1, by=nC, length=nrow(data))
for(i in 1:nC) {
  data[, cols[i]] <- split[ind + i - 1]
}
head(data)


#write.table(data, "dataforfigures_20210519.csv", quote=F, row.names=F)
#data <- read.table("dataforfigures_20200607.csv", header=T)

cogonly <- data[data$PRS == "Cog", ]
cogonly$pvalue_bonf <- p.adjust(cogonly$pvalue, "bonferroni")
noncogonly <- data[data$PRS == "NonCog", ]
noncogonly$pvalue_bonf <- p.adjust(noncogonly$pvalue, "bonferroni")
reshaped <- cbind(cogonly, noncogonly)
head(reshaped)

write.table(reshaped, "datafortable_20210519.csv", quote=F, row.names=F)




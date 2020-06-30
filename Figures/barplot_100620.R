#rosa
#bar plot with full results for indirect direct tot and ratio
rm(list = ls())

setwd("~/Dropbox/ROSA/Perline_cognoncog/results")

library(data.table)
library(psych)
library(stringr)
library(ggplot2)
library(wesanderson)
library(readr)
Trios_NTR_CITO <-read.table('summary_mean_CI_trios_NTR_CITO_20200531.csv')
rownames(Trios_NTR_CITO) <- Trios_NTR_CITO$Estimates
Trios_NTR_CITO <- Trios_NTR_CITO[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog","total_Cog","total_NonCog","ratio_tot_Cog","ratio_tot_NonCog") ] #Same variables of interest
Trios_NTR_CITO <- t(Trios_NTR_CITO)
Trios_NTR_CITO <- as.data.frame(Trios_NTR_CITO)
Trios_NTR_CITO$Measure <- row.names(Trios_NTR_CITO)
Trios_NTR_CITO$Methods <- "Trios_NTR_CITO"

Trios_NTR_EA <- read.table("summary_mean_CI_trios_NTR_EA_20200531.csv")
rownames(Trios_NTR_EA) <- Trios_NTR_EA$Estimates
Trios_NTR_EA <- Trios_NTR_EA[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog","total_Cog","total_NonCog","ratio_tot_Cog","ratio_tot_NonCog") ] #Same variables of interest
Trios_NTR_EA <- t(Trios_NTR_EA)
Trios_NTR_EA <- as.data.frame(Trios_NTR_EA)
Trios_NTR_EA$Measure <- row.names(Trios_NTR_EA)
Trios_NTR_EA$Methods <- "Trios_NTR_EA"

Sibs_NTR_CITO <- read.table("summary_mean_CI_siblings_NTR_CITO_20200511.csv")
rownames(Sibs_NTR_CITO) <- Sibs_NTR_CITO$Estimates
Sibs_NTR_CITO <- Sibs_NTR_CITO[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog","total_Cog","total_NonCog","ratio_tot_Cog","ratio_tot_NonCog") ] #Same variables of interest
Sibs_NTR_CITO <- t(Sibs_NTR_CITO)
Sibs_NTR_CITO <- as.data.frame(Sibs_NTR_CITO)
Sibs_NTR_CITO$Measure <- row.names(Sibs_NTR_CITO)
Sibs_NTR_CITO$Methods <- "Siblings_NTR_CITO"

Sibs_NTR_EA<- read.table("summary_mean_CI_siblings_NTR_EA_20200511.csv")
rownames(Sibs_NTR_EA) <- Sibs_NTR_EA$Estimates
Sibs_NTR_EA <- Sibs_NTR_EA[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog","total_Cog","total_NonCog","ratio_tot_Cog","ratio_tot_NonCog") ] #Same variables of interest
Sibs_NTR_EA <- t(Sibs_NTR_EA)
Sibs_NTR_EA <- as.data.frame(Sibs_NTR_EA)
Sibs_NTR_EA$Measure <- row.names(Sibs_NTR_EA)
Sibs_NTR_EA$Methods <- "Siblings_NTR_EA"


DZs_TEDS_12 <- read.table("summary_mean_CI_siblings_TEDS_12_270520.csv")
head(DZs_TEDS_12)
names(DZs_TEDS_12)
rownames(DZs_TEDS_12) <- DZs_TEDS_12$Estimates
DZs_TEDS_12 <- DZs_TEDS_12[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog","total_Cog","total_NonCog","ratio_tot_Cog","ratio_tot_NonCog") ] #Same variables of interest
DZs_TEDS_12 <- t(DZs_TEDS_12)
DZs_TEDS_12 <- as.data.frame(DZs_TEDS_12)
DZs_TEDS_12$Measure <- row.names(DZs_TEDS_12)
DZs_TEDS_12$Methods <- "DZs_TEDS_age_12"


DZs_TEDS_16 <- read.table("summary_mean_CI_siblings_TEDS_GCSE_270520.csv")
head(DZs_TEDS_16)
rownames(DZs_TEDS_16) <- DZs_TEDS_16$Estimates
DZs_TEDS_16 <- DZs_TEDS_16[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog","total_Cog","total_NonCog","ratio_tot_Cog","ratio_tot_NonCog") ] #Same variables of interest
DZs_TEDS_16 <- t(DZs_TEDS_16)
DZs_TEDS_16 <- as.data.frame(DZs_TEDS_16)
DZs_TEDS_16$Measure <- row.names(DZs_TEDS_16)
DZs_TEDS_16$Methods <- "DZs_TEDS_GCSE"

siblings <- read.table("summary_mean_CI_siblings_UKB_20200529.csv")
rownames(siblings) <- siblings$Estimates
siblings <- siblings[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog","total_Cog","total_NonCog","ratio_tot_Cog","ratio_tot_NonCog") ] #Same variables of interest
siblings <- t(siblings)
siblings <- as.data.frame(siblings)
siblings$Measure <- row.names(siblings)
siblings$Methods <- "Siblings_UKB_EA"

adoption <- read.table("summary_mean_CI_adoption_UKB_20200529.csv")
rownames(adoption) <- adoption$Estimates
adoption <- adoption[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog","total_Cog","total_NonCog","ratio_tot_Cog","ratio_tot_NonCog") ] #Same variables of interest
adoption <- t(adoption)
adoption <- as.data.frame(adoption)
adoption$Measure <- row.names(adoption)
adoption$Methods <- "Adoption_UKB_EA"
head(adoption)

data <- rbind(Trios_NTR_CITO, Trios_NTR_EA, Sibs_NTR_CITO, Sibs_NTR_EA, DZs_TEDS_12,DZs_TEDS_16,siblings, adoption)

#Split Measure into Type and PRS and change format back to dataframe
split <- unlist(strsplit(data$Measure, "_(?=[^_]+$)", perl=TRUE))
cols <- c("Type", "PRS")
nC <- length(cols)
ind <- seq(from=1, by=nC, length=nrow(data))
for(i in 1:nC) {
  data[, cols[i]] <- split[ind + i - 1]
}
head(data)

write.table(data, "dataforfigures_20200531.csv", quote=F, row.names=T)
data <- read.table("dataforfigures_20200531.csv", header=T)

#add pvalue in data
data$Zscore <- data$original/data$se
data$pvalue <- 2*pnorm(-abs(data$Zscore))

data$pvaluebonf <- p.adjust(data$pvalue, "bonferroni")
data$pvaluebonf

data$Type <- factor(data$Type, 
                    levels=c("total", "direct", "indirect", "ratio_tot"))
data$Methods <- factor(data$Methods, 
                       levels=c("Siblings_UKB_EA","Adoption_UKB_EA","DZs_TEDS_age_12","DZs_TEDS_GCSE","Siblings_NTR_CITO","Siblings_NTR_EA","Trios_NTR_CITO","Trios_NTR_EA"))
ggplot(data, aes(x=Type, y=original, fill=PRS)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=leftCI, ymax=rightCI),
                width=.4,                    # Width of the error bars
                position=position_dodge(.9)) + 
  scale_fill_manual(values = c("#1E90FF","#ff9933")) + 
  #coord_flip()+
  scale_x_discrete(limits = levels(data$Type)) +
  guides(fill = guide_legend(reverse = F))+
  theme_light()+
  xlab("Genetic Effects")+
  facet_grid(. ~ Methods)+
  #scale_y_continuous(expand = c(0,0)) +
  ylim(-0.6, 1.1) +
  theme(axis.title.y=element_blank(),
        legend.title=element_text("PRS"),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(colour = "black"))+
  theme(axis.text.x = element_text(size = 5)) 
  
  
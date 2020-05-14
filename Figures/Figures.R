###############################
### Figures in UKB
###############################

library(data.table)
setwd("C:/Users/Admin/Documents/GeneticNurtureNonCog")

siblings <- read.table("UKB/Siblings/summary_mean_CI_siblings_UKB_20200430.csv")
head(siblings)
siblings <- t(siblings)
siblings <- as.data.frame(siblings)
siblings$Measure <- row.names(siblings)
siblings$Methods <- "Siblings_UKB"
siblings$pheno <- "EA"

adoption <- read.table("UKB/Adoptees/summary_mean_CI_adoption_UKB_20200430.csv")
head(adoption)
adoption <- t(adoption)
adoption <- as.data.frame(adoption)
adoption$Measure <- row.names(adoption)
adoption$Methods <- "Adoption_UKB"
adoption$pheno <- "EA"
head(adoption)

ntrsib <- read.table("NTR/summary_mean_CI_siblings_NTR_EA_20200511.csv")
rownames(ntrsib) <- ntrsib$Estimates
ntrsib$Estimates <- NULL
ntrsib <- t(ntrsib)
ntrsib <- as.data.frame(ntrsib)
ntrsib$Measure <- row.names(ntrsib)
ntrsib$Methods <- "Siblings_NTR"
ntrsib$pheno <- "EA"

ntrsibcito <- read.table("NTR/summary_mean_CI_siblings_NTR_CITO_20200511.csv")
rownames(ntrsibcito) <- ntrsibcito$Estimates
ntrsibcito$Estimates <- NULL
ntrsibcito <- t(ntrsibcito)
ntrsibcito <- as.data.frame(ntrsibcito)
ntrsibcito$Measure <- row.names(ntrsibcito)
ntrsibcito$Methods <- "Siblings_NTR"
ntrsibcito$pheno <- "CITO"

ntrtrioea <- read.table("NTR/summary_mean_CI_trios_NTR_EA_20200511.csv")
rownames(ntrtrioea) <- ntrtrioea$Estimates
ntrtrioea$Estimates <- NULL
ntrtrioea <- t(ntrtrioea)
ntrtrioea <- as.data.frame(ntrtrioea)
ntrtrioea$Measure <- row.names(ntrtrioea)
ntrtrioea$Methods <- "Trios_NTR"
ntrtrioea$pheno <- "EA"

ntrtriocito <- read.table("NTR/summary_mean_CI_trios_NTR_CITO_20200511.csv")
rownames(ntrtriocito) <- ntrtriocito$Estimates
ntrtriocito$Estimates <- NULL
ntrtriocito <- t(ntrtriocito)
ntrtriocito <- as.data.frame(ntrtriocito)
ntrtriocito$Measure <- row.names(ntrtriocito)
ntrtriocito$Methods <- "Trios_NTR"
ntrtriocito$pheno <- "CITO"


data <- rbind(ntrsib, ntrsibcito, ntrtrioea, ntrtriocito)
head(data)
library(stringr)
# Limit to the Measures we want 
datafig <- data[which(data$Measure == "direct_Cog" | data$Measure == "direct_NonCog" | 
                      data$Measure == "indirect_Cog" | data$Measure == "indirect_NonCog" |
                        data$Measure == "total_Cog" | data$Measure == "total_NonCog"| 
                        data$Measure == "ratio_tot_Cog" | data$Measure == "ratio_tot_NonCog"),]
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

#write.table(data, "dataforfigures_20200501.csv", quote=F, row.names=T)


library(ggplot2)
data$Type <- factor(data$Type, levels=c("total", "direct", "indirect", "ratio_tot"))
data$Methods <- factor(data$Methods, levels=c("Siblings_NTR", "Trios_NTR"))
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
  facet_grid(. ~ Methods*pheno)+
  #scale_y_continuous(expand = c(0,0)) +
  ylim(-0.6, 1) +
  theme(axis.title.y=element_blank(),
        legend.title=element_text("PRS"),
         panel.grid.major = element_blank(),
        axis.line.x = element_line(colour = "black"))





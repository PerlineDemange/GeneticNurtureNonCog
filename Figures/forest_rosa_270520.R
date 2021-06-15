##############################################################
### Forest plot for cognoncog genetic nurture results
# may 2020
# rosa
##############################################################

rm(list = ls())
#setwd("~/Dropbox/ROSA/Perline_cognoncog/results")
setwd("C:/Users/Admin/Documents/GeneticNurtureNonCog")

library(data.table)
library(psych)
library(stringr)
library(ggplot2)
library(wesanderson)
library(readr)

Trios_NTR_CITO <-read.table("NTR/summary_mean_CI_trios_NTR_CITO_lm_20210518.csv")
rownames(Trios_NTR_CITO) <- Trios_NTR_CITO$Estimates
Trios_NTR_CITO <- Trios_NTR_CITO[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog") ] #Same variables of interest
Trios_NTR_CITO <- t(Trios_NTR_CITO)
Trios_NTR_CITO <- as.data.frame(Trios_NTR_CITO)
Trios_NTR_CITO$Measure <- row.names(Trios_NTR_CITO)
Trios_NTR_CITO$Methods <- "Trios_NTR_CITO"

Trios_NTR_EA <- read.table("NTR/summary_mean_CI_trios_NTR_EA_lm_20210518.csv")
rownames(Trios_NTR_EA) <- Trios_NTR_EA$Estimates
Trios_NTR_EA <- Trios_NTR_EA[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog") ] #Same variables of interest
Trios_NTR_EA <- t(Trios_NTR_EA)
Trios_NTR_EA <- as.data.frame(Trios_NTR_EA)
Trios_NTR_EA$Measure <- row.names(Trios_NTR_EA)
Trios_NTR_EA$Methods <- "Trios_NTR_EA"

Sibs_NTR_CITO <- read.table("NTR/summary_mean_CI_siblings_NTR_CITO_both_lm_20210519.csv")
rownames(Sibs_NTR_CITO) <- Sibs_NTR_CITO$Estimates
Sibs_NTR_CITO <- Sibs_NTR_CITO[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog") ] #Same variables of interest
Sibs_NTR_CITO <- t(Sibs_NTR_CITO)
Sibs_NTR_CITO <- as.data.frame(Sibs_NTR_CITO)
Sibs_NTR_CITO$Measure <- row.names(Sibs_NTR_CITO)
Sibs_NTR_CITO$Methods <- "Siblings_NTR_CITO"

Sibs_NTR_EA<- read.table("NTR/summary_mean_CI_siblings_NTR_EA_pop_lm_20210519_BOTH.csv")
rownames(Sibs_NTR_EA) <- Sibs_NTR_EA$Estimates
Sibs_NTR_EA <- Sibs_NTR_EA[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog") ] #Same variables of interest
Sibs_NTR_EA <- t(Sibs_NTR_EA)
Sibs_NTR_EA <- as.data.frame(Sibs_NTR_EA)
Sibs_NTR_EA$Measure <- row.names(Sibs_NTR_EA)
Sibs_NTR_EA$Methods <- "Siblings_NTR_EA"


DZs_TEDS_12 <- read.table("TEDS/summary_mean_CI_siblings_TEDS_12_pop_lm_20210519_BOTH.csv")
head(DZs_TEDS_12)
names(DZs_TEDS_12)
rownames(DZs_TEDS_12) <- DZs_TEDS_12$Estimates
DZs_TEDS_12 <- DZs_TEDS_12[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog") ] #Same variables of interest
DZs_TEDS_12 <- t(DZs_TEDS_12)
DZs_TEDS_12 <- as.data.frame(DZs_TEDS_12)
DZs_TEDS_12$Measure <- row.names(DZs_TEDS_12)
DZs_TEDS_12$Methods <- "DZs_TEDS_age_12"


DZs_TEDS_16 <- read.table("TEDS/summary_mean_CI_siblings_TEDS_16_pop_lm_20210519_BOTH.csv")
head(DZs_TEDS_16)
rownames(DZs_TEDS_16) <- DZs_TEDS_16$Estimates
DZs_TEDS_16 <- DZs_TEDS_16[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog") ] #Same variables of interest
DZs_TEDS_16 <- t(DZs_TEDS_16)
DZs_TEDS_16 <- as.data.frame(DZs_TEDS_16)
DZs_TEDS_16$Measure <- row.names(DZs_TEDS_16)
DZs_TEDS_16$Methods <- "DZs_TEDS_GCSE"

siblings <- read.table("UKB/Siblings/summary_mean_CI_siblings_UKB_pop_lm_20210519_BOTH.csv")
rownames(siblings) <- siblings$Estimates
siblings <- siblings[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog") ] #Same variables of interest
siblings <- t(siblings)
siblings <- as.data.frame(siblings)
siblings$Measure <- row.names(siblings)
siblings$Methods <- "Siblings_UKB_EA"
head(siblings)

adoption <- read.table("UKB/Adoptees/summary_mean_CI_adoption_UKB_20200529.csv")
rownames(adoption) <- adoption$Estimates
adoption <- adoption[, c("direct_NonCog", "direct_Cog","indirect_NonCog", "indirect_Cog") ] #Same variables of interest
adoption <- t(adoption)
adoption <- as.data.frame(adoption)
adoption$Measure <- row.names(adoption)
adoption$Methods <- "Adoption_UKB_EA"
head(adoption)



all <- rbind(Trios_NTR_CITO, Trios_NTR_EA, Sibs_NTR_CITO, Sibs_NTR_EA, DZs_TEDS_12,DZs_TEDS_16,siblings, adoption)


head(all)
split <- str_split_fixed(all$Measure, "_", 2)
colnames(split) <- c("Type", "PRS")
head(split)
split <- print.data.frame(as.data.frame(split), quote=FALSE)
all <- cbind(all, split)
head(all)

split2 <- str_split_fixed(all$Methods, "_", 2)
colnames(split2) <- c("method", "Sample")
head(split2)
split2 <- print.data.frame(as.data.frame(split2), quote=FALSE)

all <- cbind(all, split2)
head(all)


all$Method=ifelse(all$method=="DZs"|all$method=="Siblings", "Sibling",
                  ifelse(all$method=="Adoption", "Adoption",
                         ifelse(all$method=="Trios", "Trio", NA)))

plotdata=all[,c("original", "mean", "bias","leftCI", "rightCI", "error","Measure", "Method","Sample")]

str(plotdata)
str(plotdata$Measure)

plotdata$Measure <- factor(plotdata$Measure, levels=c("direct_NonCog", "direct_Cog", "indirect_NonCog","indirect_Cog","total_NonCog", "total_Cog","ratio_NonCog", "ratio_Cog"))
plotdata$Method <- factor(plotdata$Method, levels=c("Sibling", "Adoption", "Trio"))


data<-plotdata
data$Facet <- ifelse(data$Measure== "direct_NonCog", "NonCog.Direct",
                     ifelse(data$Measure== "direct_Cog", "Cog.Direct",
                            ifelse(data$Measure=="indirect_NonCog","NonCog.Indirect",
                                   ifelse(data$Measure=="indirect_Cog","Cog.Indirect", NA))))
data$Factor <- factor(data$Facet, levels=c("Cog.Direct", "NonCog.Direct","Cog.Indirect","NonCog.Indirect"))

# plot with original and corrected CI 
# Supp Figure X. 

plot <-
  data %>%
  dplyr::mutate(var.plot = paste0(Sample, Method)) %>%
  ggplot(.,
         aes(
           x = var.plot,
           y = original,
           colour = Sample
         ))  +
  geom_point(aes(shape= Method, color = Sample), size=4)+#,position = position_dodge(4))+
  geom_errorbar(aes(x=var.plot,ymin=leftCI, ymax=rightCI),  width=.5) +#, position=position_dodge(4)) + 
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(y = "Estimated effect of polygenic score on educational outcome", x = " ") +
  facet_grid(~Facet,  switch="both") #puts facet names at bottom
plot

plot + 
  scale_color_manual(values=c("blue", "green","yellow", "orange","red"))+
  scale_shape_manual(values=c( 16,15, 17))

# # Plot with mean and uncorrected SE 
# plot <-
#   data %>%
#   dplyr::mutate(var.plot = paste0(Sample, Method)) %>%
#   ggplot(.,
#          aes(
#            x = var.plot,
#            y = mean,
#            colour = Sample
#          ))  +
#   geom_point(aes(shape= Method, color = Sample), size=4)+#,position = position_dodge(4))+
#   geom_errorbar(aes(x=var.plot,ymin=mean-error, ymax=mean+error),  width=.5) +#, position=position_dodge(4)) + 
#   theme_bw()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())+
#   labs(y = "Estimated effect of polygenic score on educational outcome", x = " ") +
#   facet_grid(~Facet,  switch="both") #puts facet names at bottom
# 
# plot
# 
# #pick colours
# #order of colours ntr cito, ntr ea, teds 12, teds 16, uk b sib, ukb adopt
# plot + 
#   scale_color_manual(values=c("blue", "green","yellow", "orange","red"))+
#   scale_shape_manual(values=c( 16,15, 17))
# 
# 
# ggplot(data=data, aes(x= Sample, y=bias,  color=Method))+ 
#   geom_point()+ 
#   facet_grid(~Factor)



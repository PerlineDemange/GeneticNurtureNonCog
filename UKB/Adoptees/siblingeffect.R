#### get data on number of siblings 

module load pre2019
module load python/3.4.2

echo "Field,Description" > NBsib.in.csv
echo "1873,FullBro" >> NBsib.in.csv
echo "1883,FullSis" >> NBsib.in.csv
echo "3972,AdoptBro" >> NBsib.in.csv
echo "3982,AdoptSis" >> NBsib.in.csv

python3 /project/ukbaumc/UKBGWAS/phenotypes/extract_ukbb_variables.py \
/project/ukbaumc/UKBGWAS/phenotypes/ukb30545.csv \
NBsib.in.csv \
NBsib.out.csv 

head NBsib.out.csv

module unload python 

module load 2019
module load R/3.5.1-foss-2019b

R
library(data.table)
library(psych)
set.seed(42)


pheno <- read.csv("NBsib.out.csv", header=T, sep=",")
pheno$Fullbro <- pheno$X1873.0.0
pheno$Fullsis <- pheno$X1883.0.0
pheno$Adoptbro <- pheno$X3972.0.0
pheno$Adoptsis <- pheno$X3982.0.0
pheno$Fullbro[pheno$Fullbro == -3] <- NA
pheno$Fullbro[pheno$Fullbro == -1] <- NA
pheno$Fullsis[pheno$Fullsis == -3] <- NA
pheno$Fullsis[pheno$Fullsis == -1] <- NA
pheno$Adoptbro[pheno$Adoptbro == -3] <- NA
pheno$Adoptbro[pheno$Adoptbro == -1] <- NA
pheno$Adoptsis[pheno$Adoptsis == -3] <- NA
pheno$Adoptsis[pheno$Adoptsis == -1] <- NA
pheno <- pheno[-c(2:13)]

pheno <- as.data.table(pheno)
pheno[, Fulltot := ifelse(is.na(pheno$Fullbro) & is.na(pheno$Fullsis), NA_real_, rowSums(pheno[, c("Fullbro","Fullsis")], na.rm = T))]
pheno[, Adopttot := ifelse(is.na(pheno$Adoptbro) & is.na(pheno$Adoptsis), NA_real_, rowSums(pheno[, c("Adoptbro","Adoptsis")], na.rm = T))]

# Merge with adoption data 
finaladop <- read.table("../Analysis/Data_scores_adop_UKB_20200529.csv", header=T)
finalcontrol <- read.table("../Analysis/Data_scores_nonadop_UKB_20200529.csv", header=T) #6500

finaladop <- merge(finaladop, pheno, by.x="ID1", by.y="eid") #6407
finalcontrol <- merge(finalcontrol, pheno, by.x="ID1", by.y="eid")

# Get tables 
summary(finaladop) #some adoptees did answer the full sib question.. 
table(finaladop$Adopttot)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
# 2357 1725  709  336  175   99   60   35   26    8    6    8    3    3    1    1
# 16   18   20
# 2    1    2

table(finaladop$Fulltot)
# 0  1  2  4  5  6 11
# 35 12  7  5  1  1  1

finaladopwithfull <- finaladop[!is.na(finaladop$Fulltot),] #62
summary(finaladopwithfull) #all have Na in adop data 


summary(finalcontrol)#nobody answered the adopted sib questions
table(finalcontrol$Fulltot)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
# 830 2111 1599  926  436  235  142   83   46   45   15   15    6    2    2    1


#Look up effect of number of sib 

control <- lm(EA_sc~scoreNonCog_sc + scoreCog_sc + Fulltot +
               sex + array + Age + sex*Age + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=finalcontrol)
summary(control)

# Estimate Std. Error t value Pr(>|t|)
# (Intercept)     1.398334   0.122656  11.400  < 2e-16 ***
#   scoreNonCog_sc  0.201900   0.012282  16.439  < 2e-16 ***
#   scoreCog_sc     0.254080   0.012139  20.931  < 2e-16 ***
#   Fulltot        -0.060969   0.006528  -9.340  < 2e-16 ***
#   sex            -0.728400   0.170861  -4.263 2.04e-05 ***
#   array           0.078491   0.037932   2.069  0.03856 *
#   Age            -0.024357   0.002033 -11.978  < 2e-16 ***
#   PC1            -1.635250   0.748217  -2.186  0.02889 *
#   PC2            -1.683564   0.890855  -1.890  0.05883 .
# PC3            -2.928589   1.129350  -2.593  0.00953 **
#   PC4            -1.984089   1.382298  -1.435  0.15123
# PC5            -4.451401   1.512889  -2.942  0.00327 **
#   PC6            -2.416539   1.480312  -1.632  0.10263
# PC7             1.573326   1.559695   1.009  0.31314
# PC8            -1.020912   1.653070  -0.618  0.53687
# PC9             0.765414   1.621818   0.472  0.63698
# PC10           -2.044254   1.619156  -1.263  0.20680
# sex:Age         0.014650   0.002969   4.934 8.27e-07 ***
  
adop <- lm(EA_sc~scoreNonCog_sc + scoreCog_sc + Adopttot +
                  sex + array + Age + sex*Age + 
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=finaladop)
summary(adop)
    
# Estimate Std. Error t value Pr(>|t|)
# (Intercept)     1.467952   0.124665  11.775  < 2e-16 ***
#   scoreNonCog_sc  0.187510   0.013708  13.679  < 2e-16 ***
#   scoreCog_sc     0.188827   0.013321  14.175  < 2e-16 ***
#   Adopttot       -0.032709   0.007585  -4.312 1.64e-05 ***
#   sex            -0.241508   0.170963  -1.413 0.157819
# array           0.143121   0.039695   3.606 0.000314 ***
#   Age            -0.028537   0.002095 -13.619  < 2e-16 ***
#   PC1            -1.009640   0.659556  -1.531 0.125879
# PC2            -0.462486   0.896052  -0.516 0.605779



# Subset per number of sib 
res <- NULL
for (nbsib in 0:4){ 
  data <- finaladop[which(finaladop$Adopttot == nbsib),]
  model <- lm(EA_sc~scoreNonCog_sc + scoreCog_sc +
                  sex + array + Age + sex*Age +
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=data)

  est.NonCog <- summary(model)$coefficients[2,1]
  est.Cog <- summary(model)$coefficients[3,1]
  SE.NonCog <- summary(model)$coefficients[2,2]
  SE.Cog <- summary(model)$coefficients[3,2]
  t.NonCog <- summary(model)$coefficients[2,3]
  t.Cog <- summary(model)$coefficients[3,3]
  p.NonCog <- summary(model)$coefficients[2,4]
  p.Cog <- summary(model)$coefficients[3,4]
  group <- "adop"
  dat <- cbind(group, nbsib, est.NonCog, est.Cog, SE.NonCog, SE.Cog, t.NonCog, t.Cog, p.NonCog, p.Cog)
  res <- rbind(res, dat)
} 
adopres <- as.data.frame(res)


res <- NULL
for (nbsib in 0:4){ 
  data <- finalcontrol[which(finalcontrol$Fulltot == nbsib),]
  model <- lm(EA_sc~scoreNonCog_sc + scoreCog_sc +
                sex + array + Age + sex*Age +
                PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=data)
  
  est.NonCog <- summary(model)$coefficients[2,1]
  est.Cog <- summary(model)$coefficients[3,1]
  SE.NonCog <- summary(model)$coefficients[2,2]
  SE.Cog <- summary(model)$coefficients[3,2]
  t.NonCog <- summary(model)$coefficients[2,3]
  t.Cog <- summary(model)$coefficients[3,3]
  p.NonCog <- summary(model)$coefficients[2,4]
  p.Cog <- summary(model)$coefficients[3,4]
  group <- "control"
  dat <- cbind(group, nbsib, est.NonCog, est.Cog, SE.NonCog, SE.Cog, t.NonCog, t.Cog, p.NonCog, p.Cog)
  res <- rbind(res, dat)
} 
controlres <- as.data.frame(res)

result <- rbind(controlres, adopres)

write.csv(result, "siblingeffect_adop_04.csv", row.names = F)

# group nbsib            noncog               cog
# 1  control     0 0.206289022062309 0.252664356915483
# 2  control     1 0.167896384434339 0.235978399940887
# 3  control     2 0.249905111993391 0.271105118164205
# 4  control     3 0.184975297638503 0.260884113290899
# 5  control     4 0.210233578066139 0.290107185019964
# 6     adop     0 0.171649618150879 0.172412387923481
# 7     adop     1 0.207545277527636 0.207545509422723
# 8     adop     2 0.157506668032175 0.173816493327903
# 9     adop     3 0.280020692479752 0.272818100565479
# 10    adop     4 0.146842304262247 0.200611788331715



###### ON PC ############
version #R version 3.5.2
setwd("D:/Backup home computer/GeneticNurtureNonCog/UKB/Adoptees")
result <- read.csv("siblingeffect_adop_04.csv")
result

noncog <- result[,c(1:3, 5,7,9)]
names(noncog) <- c("group", "nbsib", "est", "SE", "t", "p")
noncog$PGS <- "NonCog"
cog <- result[,c(1,2,4,6,8,10)]
names(cog) <- c("group", "nbsib", "est", "SE", "t", "p")
cog$PGS <- "Cog"
result <- rbind(noncog,cog)
result
library(ggplot2)
result$group <- factor(result$group, levels = c("adop", "control"),
                  labels = c("Adopted", "Non-adopted"))

# Figure

ggplot(result, aes(x=nbsib, y=est, color=PGS)) + 
  geom_point(stat="identity", position=position_dodge(0.3), size=5) +
  geom_smooth(method='lm', formula= y~x, se=F) +
  geom_errorbar(aes(ymin=est-1.96*SE, ymax=est+1.96*SE),
                width=.4,position=position_dodge(0.3), size=1) + 
  scale_color_manual(values = c("#1E90FF","#ff9933")) + 
  #coord_flip()+
  #scale_x_discrete(limits = levels(data$Type)) +
  guides(fill = guide_legend(reverse = F))+
  theme_light()+
  xlab("Number of siblings (Adopted/Full-Siblings)") +
  ylab("Estimated effect of PGS on EA") +
  facet_grid(. ~ group, scales = "free", space = "free")+
  #scale_y_continuous(expand = c(0,0)) +
  #ylim(-0.6, 1) +
  theme(strip.text = element_text(size=10, colour = 'black', face = "bold"),legend.title=element_text("PGS"),
        strip.background = element_rect(colour="grey", fill="white"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())



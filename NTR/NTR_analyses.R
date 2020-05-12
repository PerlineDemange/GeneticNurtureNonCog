
## Project: Genetic Nuture of NonCog 2020 ####
## Starting Date : 20/03/20 ####
## Script: All analyses in NTR ####
## Author : Perline Demange ####

version # At Home the R version is R version 3.6.3 (2020-02-29)
library(data.table)
library(sjlabelled)
library(tidyverse)
library(dplyr)
library("gee")
library(nlme)
library(boot)
set.seed(42)
setwd("C:/Users/Admin/Documents/GeneticNurtureNonCog/NTR")

# 1. Load Data ##################################

# Get phenotypic data from NTR
pheno <- fread('Phenotypic_data/NTR-DAC-1927_Demange_parental_env_EA_20200429.csv', header=T, colClasses=c("fisnumber"="character")) #data was changed from spss to csv in spss #FISNumber is read as integer64 so fix to to be read as chr 
head(pheno)
str(pheno)
summary(pheno) #24984

names(pheno)[names(pheno) == 'fisnumber'] <- 'FISNumber'

# do with old pheno to test
#pheno <- read.csv2("D:/Documents/CogNonCog/NTRphenotypes/educational_attainment_NTR.csv", sep=";")
#colnames(pheno)[colnames(pheno)=="?..FISNumber"] <- "FISNumber" # all first columns are not weel named when exporting csv from spss in csv 

# get variable and value labels
label <- read_spss('Phenotypic_data/NTR-DAC-1927_Demange_parental_env_EA_20200429.sav')
label.var <- get_label(label)
label.val <- get_labels(label,  values = "as.name")
label.val

#Load covariate file 
cov <- fread("Scores_NTR/MRG10_Pedigree_V2_NTR_only.sav/MRG10_Pedigree_V2_NTR_only.csv", colClasses=c("ID"="character"))
head(cov)
str(cov)
cov <- as.data.frame(cov)



# * 1.1 Clean and adjust phenotypic data: EA + CITO ======================
# * * 1.1.1 Clean education data, recode to EduYears -----------------------------

summary(is.na(pheno$educat_a)) #17891
summary(is.na(pheno$educat_c)) # same number of missing data for both 
#Change values to NA
pheno$educat_a[pheno$educat_a == -9] <- NA
pheno$educat_a[pheno$educat_a == -8] <- NA
pheno$educat_a[pheno$educat_a == -2] <- NA
pheno$educat_a[pheno$educat_a == -1] <- NA
pheno$educat_c[pheno$educat_c == -9] <- NA
pheno$educat_c[pheno$educat_c == -8] <- NA
pheno$educat_c[pheno$educat_c == -2] <- NA
pheno$educat_c[pheno$educat_c == -1] <- NA

summary(is.na(pheno$educat_a)) #13971
summary(is.na(pheno$educat_c)) #12293 for educat_c 

# recode education in ISCED following Okbay et al 2016 Supplementary Table 1.3 (only use category c)
pheno$ISCED <- pheno$educat_c
pheno$ISCED[pheno$ISCED == 0] <- "ISCED0"
pheno$ISCED[pheno$ISCED == 1] <- "ISCED1"
pheno$ISCED[pheno$ISCED == 2] <- "ISCED2"
pheno$ISCED[pheno$ISCED == 3] <- "ISCED2"
pheno$ISCED[pheno$ISCED == 4] <- "ISCED3"
pheno$ISCED[pheno$ISCED == 5] <- "ISCED3"
pheno$ISCED[pheno$ISCED == 6] <- "ISCED5"
pheno$ISCED[pheno$ISCED == 7] <- "ISCED5"
pheno$ISCED[pheno$ISCED == 8] <- "ISCED6"
# recode in EduYears following Okbay et al 2016 Supplementary Table 1.6
pheno$Eduyears <-  pheno$ISCED
pheno$Eduyears[pheno$Eduyears == "ISCED0"] <- 1
pheno$Eduyears[pheno$Eduyears == "ISCED1"] <- 7
pheno$Eduyears[pheno$Eduyears == "ISCED2"] <- 10
pheno$Eduyears[pheno$Eduyears == "ISCED3"] <- 13
pheno$Eduyears[pheno$Eduyears== "ISCED5"] <- 19
pheno$Eduyears[pheno$Eduyears == "ISCED6"] <- 22
pheno$Eduyears <- as.numeric(pheno$Eduyears)

hist(pheno$Eduyears)

# * * 1.1.2 Clean CITO ------------------------------------
# Cito score is one variable cito_final
hist(pheno$cito_final)

# 2. Sibling Analyses ##################################
# * 2.1 Clean Data for Sibling analysis =======================
# Order by family number 
phenoorder <- pheno[order(pheno$FamilyNumber), ]
summary(phenoorder)

# Get Only unique occurence of an individual
phenounique <- unique(phenoorder)
#View(phenounique)

# Get number of individuals in each family
number <- as.data.frame(table(phenounique$FamilyNumber))
number
onemember <- number[which(number$Freq == 1),] #identify family with only one member in this data #2003
head(onemember)
list_onemember <- onemember$Var1 

# Get data for which several members of family are present 
datafamily <- phenounique[!which(phenounique$FamilyNumber %in% list_onemember), ] #22981
head(datafamily)


# Select only siblings
label.val$Extension
# 1 to 6 is multiple
# 10 to 24 is sibling #no value between 6 and 10
allsib <- datafamily[which(datafamily$Extension > 0 & datafamily$Extension <24),] #14803



# * * 2.1.1 Remove MZ: we can not use pairs that are only MZ because their PGS is identical ------------------------------------
# get all multiple 
multiple <- datafamily[which(datafamily$Extension > 0 & datafamily$Extension <6),] #12002
label.val$twzyg
summary(multiple$twzyg) # 55 NAs
MZ <- multiple[which(multiple$twzyg == 1 | multiple$twzyg == 3),] #6654
DZ <- multiple[which(multiple$twzyg == 2 | multiple$twzyg > 3),] #5293
summary(DZ$multiple_type)

# investigate NAs in zygosity
head(multiple[is.na(multiple$twzyg),])
summary(multiple[is.na(multiple$twzyg),]$famzyg) # all are MZ or DZ or NA
multiple[is.na(multiple$twzyg),][is.na(multiple[is.na(multiple$twzyg),]$famzyg)] 
# multiple with both zygosity and fam zygosity missing are in different family, except 2 triplets in fam 21783
# Get multiple with family zygosity being MZ (while zygosity is na)
MZ_multipleNa <- multiple[is.na(multiple$twzyg),][which(multiple[is.na(multiple$twzyg),]$famzyg == 1 | multiple[is.na(multiple$twzyg),]$famzyg == 3),]
MZ_multipleNa # all of them are part of triplets # 19

MZ <- rbind(MZ, MZ_multipleNa) #6673

## Several options can be investigated
## 1. Remove all MZ 
## 2. Take one MZ per pair randomly 
## If the sample size is not so much lower when removing all MZ we will remove all MZ 

# 1. Remove all MZ 
# mzlist <- MZ$FISNumber
# datawoMZ <- allsib[!which(allsib$FISNumber %in% mzlist),] #12704
# 
# number <- as.data.frame(table(datawoMZ$FamilyNumber))
# number
# onemember <- number[which(number$Freq == 1),] #identify family with only one member in this data #397
# head(onemember)
# list_onemember <- onemember$Var1 
# datawoMZ_sibonly<- datawoMZ[!which(datawoMZ$FamilyNumber %in% list_onemember), ] #11677
# head(datawoMZ_sibonly)

# By removing the MZ twins we had to remove an additional 397 families: 
# 397 sibships where made of a MZ pair + one sibling
# So if we take one MZ per pair, we would increase the total sample size by 397*2= 794 to 12471
#1. This is not substancial so I just remove all MZ pairs 
#allsib <- datawoMZ_sibonly

#2. Keep one MZ per pair
#This might be substancial so I keep one MZ per pair

oneMZ <- MZ %>% group_by(FamilyNumber) %>% sample_n(1) # take at random one of MZ per family 
oneMZ <- as.data.frame(oneMZ)
a <- unique(oneMZ$FamilyNumber) #3411, same as oneMZ so only one per family 
pickoneMZ <- oneMZ$FISNumber
datawithrestMZ <- MZ[!which(MZ$FISNumber %in% pickoneMZ),] #3262
datawithrestMZvec <- datawithrestMZ$FISNumber
dataoneMZ <- allsib[!which(allsib$FISNumber %in% datawithrestMZvec),] # 11541

number <- as.data.frame(table(dataoneMZ$FamilyNumber))
number
onemember <- number[which(number$Freq == 1),] #identify family with only one member in this data #2754
head(onemember)
list_onemember <- onemember$Var1 
dataoneMZ_sibonly<- dataoneMZ[!which(dataoneMZ$FamilyNumber %in% list_onemember), ] #8787
head(dataoneMZ_sibonly)

allsib <- dataoneMZ_sibonly

# * * 2.1.2 Add scores for all individuals  --------------------------------

noncog <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_ALL_NONCOG_CEULDpred_inf.profile", header=T, stringsAsFactors =F, colClasses=c("IID"="character"))
head(noncog)
cog <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_ALL_COG_CEULDpred_inf.profile", header=T, stringsAsFactors =F, colClasses=c("IID"="character"))
head(cog)
head(allsib)
scores <- merge(noncog, cog, by= c("FID", "IID", "PHENO", "CNT", "CNT2"), suffixes = c( ".NonCog", ".Cog"))
str(scores)

datasib <- merge(allsib, scores, by.x = "FISNumber", by.y= "IID")
head(datasib)  #8787

datasib$SCORE.Cog <- -1 * datasib$SCORE.Cog
datasib$SCORE.NonCog <- -1 * datasib$SCORE.NonCog

# Add covariates (PC + chips) 
#cov_label <- read_spss("Scores_NTR/MRG10_Pedigree_V2_NTR_only.sav/MRG10_Pedigree_V2_NTR_only.sav") 
#label.val <- get_labels(cov_label,  values = "as.name")
#label.val

data <- merge(datasib, cov, by.x= "FISNumber", by.y= "ID")
head(data)
summary(data$Platform) 

datasib <- data

# * 2.2 Sibling analysis with EA ==========================
# * * 2.2.1 Clean data and scale variables ------
# Keep only data with EA
datasibEA <- datasib[!is.na(datasib$Eduyears),] #3129
head(datasibEA)

number <- as.data.frame(table(datasibEA$FamilyNumber))
number
onemember <- number[which(number$Freq == 1),] #identify family with only one member in this data #695
list_onemember <- onemember$Var1 
datasibEA_sibonly<- datasibEA[!which(datasibEA$FamilyNumber %in% list_onemember), ] #2438
head(datasibEA_sibonly)

finalsib <- as.data.frame(datasibEA_sibonly)
# sample size = 2438

# scale variables
finalsib[,c("EA_sc","scoreNonCog_sc", "scoreCog_sc")]<-apply(finalsib[,c("Eduyears","SCORE.NonCog", "SCORE.Cog")],
                                                             2,
                                                             scale)
hist(finalsib$EA_sc)

# save data
# write.table(finalsib, "Data_siblings_NTR_EA_20200511.csv", row.names=F, quote=F)
# finalsib <- fread("Data_siblings_NTR_EA_20200511.csv")
# head(finalsib)



# * * 2.2.2 Simple linear model -----

global <- lm(EA_sc ~ scoreNonCog_sc + scoreCog_sc + sex + yob + sex*yob + 
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=finalsib)
summary(global)

# * * 2.2.3 Create between-family and within-family estimates of PGS -----
# Between-family estimate = average per family
meanNC<-group_by(finalsib,FamilyNumber) %>% summarize(m=mean(scoreNonCog_sc))
colnames(meanNC) <- c("FamilyNumber", "GPS_B_NonCog")
meanC<-group_by(finalsib,FamilyNumber) %>% summarize(m=mean(scoreCog_sc)) 
colnames(meanC) <- c("FamilyNumber", "GPS_B_Cog")
finalsib<-merge(finalsib,meanNC,by="FamilyNumber")
finalsib<-merge(finalsib,meanC,by="FamilyNumber")

# Within-family estimates
finalsib$GPS_W_NonCog <- finalsib$scoreNonCog_sc  - finalsib$GPS_B_NonCog  
finalsib$GPS_W_Cog <- finalsib$scoreCog_sc  - finalsib$GPS_B_Cog  


# * * 2.2.4 Run mixed model between-within regression -----
final <- lme(EA_sc~GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog + 
               sex + yob + sex*yob +
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FamilyNumber, method="ML", na.action=na.omit,data=finalsib)
summary(final)


# Extract estimates 
names(summary(final))
summary(final)$tTable
# Value    Std.Error   DF    t-value      p-value
# (Intercept)  -1.294493e+01  9.172725960 1384 -1.4112414 1.583982e-01
# GPS_B_NonCog  2.700489e-01  0.026677347 1035 10.1227806 4.959911e-23
# GPS_B_Cog     2.587895e-01  0.026813999 1035  9.6512816 3.662399e-21
# GPS_W_NonCog  9.891056e-02  0.028762993 1384  3.4388131 6.016225e-04
# GPS_W_Cog     1.339745e-01  0.028051887 1384  4.7759518 1.978682e-06

# write.table(summary(final)$tTable, "Estimates_Siblings_NTR_EA_20200511.csv", quote=F )

direct_NonCog <- summary(final)$tTable[4,1]# direct effect is beta within 
direct_Cog <- summary(final)$tTable[5,1]
total_NonCog <- summary(final)$tTable[2,1] # total is beta between 
total_Cog <- summary(final)$tTable[3,1]
indirect_NonCog <- total_NonCog - direct_NonCog 
indirect_Cog <- total_Cog - direct_Cog 
ratio_NonCog <- indirect_NonCog/direct_NonCog  
ratio_Cog <- indirect_Cog/direct_Cog  
ratio_tot_NonCog <- indirect_NonCog/total_NonCog 
ratio_tot_Cog <- indirect_Cog/total_Cog 

# * * 2.2.5 Bootstrapping ------
nboot <- 10000
bootcoef<-function(data,index){
  datx<-data[index,]
  mod<-lme(EA_sc~GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog + 
               sex + yob + sex*yob +
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FamilyNumber, method="ML", na.action=na.omit, data=datx, control=lmeControl(opt = "optim"))
  fixef(mod) #get fixed effects
}

# carry out bootstrap
boot.out<-boot(finalsib,bootcoef,nboot, parallel = "multicore", ncpus=20) 

#saveRDS(boot.out, "bootstrapped_output_siblings_NTR_EA_20200511.Rda")

#plot to check bootstrapping
png("NTR.sib.EA.bootstrap_20200511.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out)
dev.off()


#Save t output of boot
bootoutput <- as.data.frame(boot.out$t)
colnames(bootoutput) <- rownames(as.data.frame(boot.out$t0))
head(bootoutput)
#write.table(bootoutput, "Data_scores_siblings_NTR_EA_bootstrapped_20200511.csv", row.names=F, quote=F)

# get values out of boot.out for all estimates + create indirect and ratio estimates
original <- as.data.frame(t(boot.out$t0)) # estimates of the original sample #best estimates of the effects

original$direct_NonCog <- original$GPS_W_NonCog
original$direct_Cog <- original$GPS_W_Cog
original$total_NonCog <- original$GPS_B_NonCog
original$total_Cog <- original$GPS_B_Cog
original$indirect_NonCog <- original$total_NonCog - original$direct_NonCog
original$indirect_Cog <- original$total_Cog - original$direct_Cog
original$ratio_NonCog <- original$indirect_NonCog / original$direct_NonCog
original$ratio_Cog <- original$indirect_Cog / original$direct_Cog
original$ratio_tot_NonCog <- original$indirect_NonCog / original$total_NonCog
original$ratio_tot_Cog <- original$indirect_Cog / original$total_Cog


bootoutput$direct_NonCog <- bootoutput$GPS_W_NonCog
bootoutput$direct_Cog <- bootoutput$GPS_W_Cog
bootoutput$total_NonCog <- bootoutput$GPS_B_NonCog
bootoutput$total_Cog <- bootoutput$GPS_B_Cog
bootoutput$indirect_NonCog <- bootoutput$total_NonCog - bootoutput$direct_NonCog
bootoutput$indirect_Cog <- bootoutput$total_Cog - bootoutput$direct_Cog
bootoutput$ratio_NonCog <- bootoutput$indirect_NonCog / bootoutput$direct_NonCog
bootoutput$ratio_Cog <- bootoutput$indirect_Cog / bootoutput$direct_Cog
bootoutput$ratio_tot_NonCog <- bootoutput$indirect_NonCog / bootoutput$total_NonCog
bootoutput$ratio_tot_Cog <- bootoutput$indirect_Cog / bootoutput$total_Cog


mean <- apply(bootoutput, 2, mean) # mean of the estimates of the bootstrap resamples
bias <- mean - original
se <- apply(bootoutput, 2, sd) #the standard deviation of the bootstrap estimates is the standard error of the sample estimates

error <- qnorm(0.975)*se
leftCI <- original - bias - error # normal Ci from boot.ci 
rightCI <- original - bias + error
# Other kind of CI given by boot.ci, not saved
# leftCI3 <- quantile(bootoutput$X, 0.025) # percentile CI from boot.ci 
# rightCI3 <- quantile(bootoutput$X, 0.975)
# leftCI4 <- 2*original$SCORE.Nontrans.Cog_sc - quantile(bootoutput$SCORE.Nontrans.Cog_sc, 0.975) #basic ci from boot.ci
# rightCI4 <- 2*original$SCORE.Nontrans.Cog_sc - quantile(bootoutput$SCORE.Nontrans.Cog_sc, 0.025)

statsoutput <- rbind(original, mean, bias, se, error, leftCI, rightCI)
statsoutput$Estimates <- c('original', 'mean', 'bias', 'se', 'error', 'leftCI', 'rightCI')
statsoutput
tot <- statsoutput[,c(ncol(statsoutput), 1:(ncol(statsoutput)-1))]
tot

#write.table(tot, "summary_mean_CI_siblings_NTR_EA_20200511.csv", row.names=T, quote=F)

bmain <- bootoutput
#Perform  t-tests
t.test(bmain$indirect_NonCog, bmain$direct_NonCog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_NonCog and bmain$direct_NonCog
# t = 158.78, df = 19628, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.07606803 0.07796959
# sample estimates:
#   mean of x  mean of y 
# 0.17362816 0.09660936 

t.test(bmain$indirect_Cog, bmain$direct_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_Cog and bmain$direct_Cog
# t = -7.8764, df = 19617, p-value = 3.543e-15
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.004598602 -0.002765914
# sample estimates:
#   mean of x mean of y 
# 0.1264471 0.1301294 

t.test(bmain$indirect_NonCog, bmain$indirect_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_NonCog and bmain$indirect_Cog
# t = 92.831, df = 19972, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.04618485 0.04817727
# sample estimates:
#   mean of x mean of y 
# 0.1736282 0.1264471 

t.test(bmain$ratio_NonCog, bmain$ratio_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$ratio_NonCog and bmain$ratio_Cog
# t = 22.071, df = 10319, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   1.061719 1.268691
# sample estimates:
#   mean of x mean of y 
# 2.272622  1.107417 

t.test(bmain$ratio_tot_NonCog, bmain$ratio_tot_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$ratio_tot_NonCog and bmain$ratio_tot_Cog
# t = 86.367, df = 19979, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.1470326 0.1538613
# sample estimates:
#   mean of x mean of y 
# 0.6408777 0.4904307 
# 

# * 2.3 Siblings analyses with CITO ==============================
#* * 2.3.1 Clean data and scale variables -----
datasibcito <- datasib[!is.na(datasib$cito_final),] #1850

# Process again to exclude individual without sib 
number <- as.data.frame(table(datasibcito$FamilyNumber))
onemember <- number[which(number$Freq == 1),] #identify family with only one member in this data #219
head(onemember)
list <- onemember$Var1 
datasibcito <- datasibcito[!which(datasibcito$FamilyNumber %in% list), ] #1631
head(datasibcito)

finalsib <- as.data.frame(datasibcito)

summary(is.na(finalsib$cito_final))  #1631: final sample size


# scale variables
finalsib[,c("CITO_sc","scoreNonCog_sc", "scoreCog_sc")]<-apply(finalsib[,c("cito_final","SCORE.NonCog", "SCORE.Cog")],
                                                             2,
                                                             scale)
# check how data looks 
hist(finalsib$cito_final)
hist(finalsib$CITO_sc)
hist(finalsib$SCORE.Cog)
hist(finalsib$scoreCog_sc)

# save data
#write.table(finalsib, "Data_siblings_NTR_CITO_20200511.csv", row.names=F, quote=F)
#finalsib <- fread("Data_siblings_NTR_CITO_20200511.csv")

# * * 2.3.2 Simple linear model -----

global <- lm(CITO_sc ~ scoreNonCog_sc + scoreCog_sc + sex + yob + sex*yob +
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=finalsib)
summary(global)

# * * 2.3.3 Create between-family and within-family estimates of PGS -------------
# Between-family estimate = average per family
meanNC<-group_by(finalsib,FamilyNumber) %>% summarize(m=mean(scoreNonCog_sc))
colnames(meanNC) <- c("FamilyNumber", "GPS_B_NonCog")
meanC<-group_by(finalsib,FamilyNumber) %>% summarize(m=mean(scoreCog_sc)) 
colnames(meanC) <- c("FamilyNumber", "GPS_B_Cog")
finalsib<-merge(finalsib,meanNC,by="FamilyNumber")
finalsib<-merge(finalsib,meanC,by="FamilyNumber")

# Within-family estimates
finalsib$GPS_W_NonCog <- finalsib$scoreNonCog_sc  - finalsib$GPS_B_NonCog  
finalsib$GPS_W_Cog <- finalsib$scoreCog_sc  - finalsib$GPS_B_Cog  

# * * 2.3.4 Run mixed model between-within regression ------------------
final <- lme(CITO_sc~GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog + 
               sex + yob + sex*yob + 
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FamilyNumber, method="ML", na.action=na.omit,data=finalsib)
summary(final)

# Extract estimates 
names(summary(final))
summary(final)$tTable
# Value    Std.Error  DF     t-value      p-value
# (Intercept)   1.416429e+01  42.93422408 858  0.32990677 7.415509e-01
# GPS_B_NonCog  1.759328e-01   0.03406690 754  5.16433315 3.090748e-07
# GPS_B_Cog     3.212078e-01   0.03382041 754  9.49745298 2.772288e-20
# GPS_W_NonCog  1.845122e-01   0.03782887 858  4.87754932 1.280036e-06
# GPS_W_Cog     2.326259e-01   0.03866158 858  6.01697913 2.629687e-09

#write.table(summary(final)$tTable, "Estimates_Siblings_NTR_CITO_20200511.csv", quote=F )


direct_NonCog <- summary(final)$tTable[4,1]# direct effect is beta within 
direct_Cog <- summary(final)$tTable[5,1]
total_NonCog <- summary(final)$tTable[2,1] # total is beta between 
total_Cog <- summary(final)$tTable[3,1]
indirect_NonCog <- total_NonCog - direct_NonCog 
indirect_Cog <- total_Cog - direct_Cog 
ratio_NonCog <- indirect_NonCog/direct_NonCog  
ratio_Cog <- indirect_Cog/direct_Cog  
ratio_tot_NonCog <- indirect_NonCog/total_NonCog
ratio_tot_Cog <- indirect_Cog/total_Cog 

# * * 2.3.5 Bootstrapping ------------
nboot <- 10000
bootcoef<-function(data,index){
  datx<-data[index,]
  mod<-lme(CITO_sc~GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog + 
               sex + yob + sex*yob +
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FamilyNumber, method="ML", na.action=na.omit, data=datx, control=lmeControl(opt = "optim"))
  fixef(mod) #get fixed effects
}

# carry out bootstrap
boot.out<-boot(finalsib,bootcoef,nboot, parallel = "multicore", ncpus=20) 

#saveRDS(boot.out, "bootstrapped_output_siblings_NTR_CITO_20200511.Rda")

#plot to check bootstrapping
png("NTR.sib.CITO.bootstrap_20200511.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out)
dev.off()

#Save t output of boot
bootoutput <- as.data.frame(boot.out$t)
colnames(bootoutput) <- rownames(as.data.frame(boot.out$t0))
head(bootoutput)
#write.table(bootoutput, "Data_scores_siblings_NTR_CITO_bootstrapped_20200511.csv", row.names=F, quote=F)

# get values out of boot.out for all estimates + create indirect and ratio estimates
original <- as.data.frame(t(boot.out$t0)) # estimates of the original sample #best estimates of the effects

original$direct_NonCog <- original$GPS_W_NonCog
original$direct_Cog <- original$GPS_W_Cog
original$total_NonCog <- original$GPS_B_NonCog
original$total_Cog <- original$GPS_B_Cog
original$indirect_NonCog <- original$total_NonCog - original$direct_NonCog
original$indirect_Cog <- original$total_Cog - original$direct_Cog
original$ratio_NonCog <- original$indirect_NonCog / original$direct_NonCog
original$ratio_Cog <- original$indirect_Cog / original$direct_Cog
original$ratio_tot_NonCog <- original$indirect_NonCog / original$total_NonCog
original$ratio_tot_Cog <- original$indirect_Cog / original$total_Cog

bootoutput$direct_NonCog <- bootoutput$GPS_W_NonCog
bootoutput$direct_Cog <- bootoutput$GPS_W_Cog
bootoutput$total_NonCog <- bootoutput$GPS_B_NonCog
bootoutput$total_Cog <- bootoutput$GPS_B_Cog
bootoutput$indirect_NonCog <- bootoutput$total_NonCog - bootoutput$direct_NonCog
bootoutput$indirect_Cog <- bootoutput$total_Cog - bootoutput$direct_Cog
bootoutput$ratio_NonCog <- bootoutput$indirect_NonCog / bootoutput$direct_NonCog
bootoutput$ratio_Cog <- bootoutput$indirect_Cog / bootoutput$direct_Cog
bootoutput$ratio_tot_NonCog <- bootoutput$indirect_NonCog / bootoutput$total_NonCog
bootoutput$ratio_tot_Cog <- bootoutput$indirect_Cog / bootoutput$total_Cog

mean <- apply(bootoutput, 2, mean) # mean of the estimates of the bootstrap resamples
bias <- mean - original
se <- apply(bootoutput, 2, sd) #the standard deviation of the bootstrap estimates is the standard error of the sample estimates

error <- qnorm(0.975)*se
leftCI <- original - bias - error # normal Ci from boot.ci 
rightCI <- original - bias + error

statsoutput <- rbind(original, mean, bias, se, error, leftCI, rightCI)
statsoutput$Estimates <- c('original', 'mean', 'bias', 'se', 'error', 'leftCI', 'rightCI')
statsoutput
tot <- statsoutput[,c(ncol(statsoutput), 1:(ncol(statsoutput)-1))]
tot

#write.table(tot, "summary_mean_CI_siblings_NTR_CITO_20200511.csv", row.names=T, quote=F)

bmain <- bootoutput

#Perform  t-tests
t.test(bmain$indirect_NonCog, bmain$direct_NonCog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_NonCog and bmain$direct_NonCog
# t = -317.14, df = 19761, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.1990520 -0.1966066
# sample estimates:
#   mean of x   mean of y 
# -0.01000833  0.18782097 

t.test(bmain$indirect_Cog, bmain$direct_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_Cog and bmain$direct_Cog
# t = -216.25, df = 19731, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.1427061 -0.1401424
# sample estimates:
#   mean of x  mean of y 
# 0.09224501 0.23366928 

t.test(bmain$indirect_NonCog, bmain$indirect_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_NonCog and bmain$indirect_Cog
# t = -151.66, df = 19947, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.1035749 -0.1009318
# sample estimates:
#   mean of x   mean of y 
# -0.01000833  0.09224501 

t.test(bmain$ratio_NonCog, bmain$ratio_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$ratio_NonCog and bmain$ratio_Cog
# t = -102.86, df = 19844, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.4567869 -0.4397035
# sample estimates:
#   mean of x    mean of y 
# 0.0008019466 0.4490471871 

t.test(bmain$ratio_tot_NonCog, bmain$ratio_tot_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$ratio_tot_NonCog and bmain$ratio_tot_Cog
# t = -113.37, df = 15063, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.3586261 -0.3464354
# sample estimates:
#   mean of x   mean of y 
# -0.07288037  0.27965037 


# 3.  Trios analyses #####################################
# * 3.1 Clean data for trios analyses ===============================
# * * 3.1.1 Get transmitted and non-transmitted scores ----------
# only upload IID and score 
noncog_dad_trans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_DAD_Trans_NONCOG_CEULDpred_inf.profile", header=T, stringsAsFactors =F)[,c(2,6)]
noncog_dad_nontrans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_DAD_NonTrans_NONCOG_CEULDpred_inf.profile", header=T, stringsAsFactors =F)[,c(2,6)]
noncog_mom_trans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_MOM_Trans_NONCOG_CEULDpred_inf.profile", header=T, stringsAsFactors =F)[,c(2,6)]
noncog_mom_nontrans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_MOM_NonTrans_NONCOG_CEULDpred_inf.profile", header=T, stringsAsFactors =F)[,c(2,6)]

cog_dad_trans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_DAD_Trans_COG_CEULDpred_inf.profile", header=T, stringsAsFactors =F)[,c(2,6)]
cog_dad_nontrans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_DAD_NonTrans_COG_CEULDpred_inf.profile", header=T, stringsAsFactors =F)[,c(2,6)]
cog_mom_trans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_MOM_Trans_COG_CEULDpred_inf.profile", header=T, stringsAsFactors =F)[,c(2,6)]
cog_mom_nontrans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_MOM_NonTrans_COG_CEULDpred_inf.profile", header=T, stringsAsFactors =F)[,c(2,6)]

head(cog_dad_nontrans)

# need to remove suffixes in IID
cog_dad_trans$IID <- str_replace(cog_dad_trans$IID, "_T", "")
cog_mom_trans$IID <- str_replace(cog_mom_trans$IID, "_T", "")
noncog_dad_trans$IID <- str_replace(noncog_dad_trans$IID, "_T", "")
noncog_mom_trans$IID <- str_replace(noncog_mom_trans$IID, "_T", "")
cog_dad_nontrans$IID <- str_replace(cog_dad_trans$IID, "_U", "")
cog_mom_nontrans$IID <- str_replace(cog_mom_trans$IID, "_U", "")
noncog_dad_nontrans$IID <- str_replace(noncog_dad_trans$IID, "_U", "")
noncog_mom_nontrans$IID <- str_replace(noncog_mom_trans$IID, "_U", "")


noncog_dad <-  merge(noncog_dad_trans, noncog_dad_nontrans, by=  "IID", suffixes = c( ".Trans", ".Nontrans"))
noncog_mom <-  merge(noncog_mom_trans, noncog_mom_nontrans, by= "IID", suffixes = c( ".Trans", ".Nontrans"))
noncog <- merge(noncog_dad, noncog_mom, by= "IID", suffixes = c( ".Dad", ".Mom"))
head(noncog)


cog_dad <-  merge(cog_dad_trans, cog_dad_nontrans, by=  "IID", suffixes = c( ".Trans", ".Nontrans"))
cog_mom <-  merge(cog_mom_trans, cog_mom_nontrans, by= "IID", suffixes = c( ".Trans", ".Nontrans"))
cog <- merge(cog_dad, cog_mom, by= "IID", suffixes = c( ".Dad", ".Mom"))
head(cog)

scores<- merge(noncog, cog, by="IID", suffixes = c(".NonCog", ".Cog"))
str(scores)               

datatrios <- merge(pheno, scores, by.x = "FISNumber", by.y= "IID") #7448
head(datatrios)
datatrios <- as.data.frame(datatrios)

# add covariates
data <- merge(datatrios, cov, by.x= "FISNumber", by.y= "ID")
head(data)
str(data)
datatrios <- data

# Reverse scores
f1 <- function(x) (x*-1)
datatrios[,c("SCORE.Nontrans.Dad.Cog", "SCORE.Nontrans.Mom.Cog", "SCORE.Trans.Dad.Cog", 
             "SCORE.Trans.Mom.Cog", "SCORE.Nontrans.Dad.NonCog", "SCORE.Nontrans.Mom.NonCog", "SCORE.Trans.Dad.NonCog", 
             "SCORE.Trans.Mom.NonCog")]<-apply(datatrios[,c("SCORE.Nontrans.Dad.Cog", "SCORE.Nontrans.Mom.Cog", 
                                                               "SCORE.Trans.Dad.Cog", "SCORE.Trans.Mom.Cog", "SCORE.Nontrans.Dad.NonCog", "SCORE.Nontrans.Mom.NonCog", 
                                                               "SCORE.Trans.Dad.NonCog", "SCORE.Trans.Mom.NonCog")],
                                                  2,
                                                  f1)


# scale scores 
datatrios[,c("CITO_sc", "EA_sc","SCORE.Nontrans.Dad.Cog_sc", "SCORE.Nontrans.Mom.Cog_sc", "SCORE.Trans.Dad.Cog_sc", 
             "SCORE.Trans.Mom.Cog_sc", "SCORE.Nontrans.Dad.NonCog_sc", "SCORE.Nontrans.Mom.NonCog_sc", "SCORE.Trans.Dad.NonCog_sc", 
             "SCORE.Trans.Mom.NonCog_sc")]<-apply(datatrios[,c("cito_final", "Eduyears","SCORE.Nontrans.Dad.Cog", "SCORE.Nontrans.Mom.Cog", 
                                                               "SCORE.Trans.Dad.Cog", "SCORE.Trans.Mom.Cog", "SCORE.Nontrans.Dad.NonCog", "SCORE.Nontrans.Mom.NonCog", 
                                                               "SCORE.Trans.Dad.NonCog", "SCORE.Trans.Mom.NonCog")],
                                                  2,
                                                  scale)

hist(datatrios$SCORE.Nontrans.Dad.Cog_sc)

# Get PGS tranmsitted and non-transmitted for both parents combined
datatrios$SCORE.Nontrans.Cog_sc <- datatrios$SCORE.Nontrans.Dad.Cog_sc + datatrios$SCORE.Nontrans.Mom.Cog_sc
datatrios$SCORE.Trans.Cog_sc <- datatrios$SCORE.Trans.Dad.Cog_sc + datatrios$SCORE.Trans.Mom.Cog_sc
datatrios$SCORE.Nontrans.NonCog_sc <- datatrios$SCORE.Nontrans.Dad.NonCog_sc + datatrios$SCORE.Nontrans.Mom.NonCog_sc
datatrios$SCORE.Trans.NonCog_sc <- datatrios$SCORE.Trans.Dad.NonCog_sc + datatrios$SCORE.Trans.Mom.NonCog_sc
hist(datatrios$SCORE.Nontrans.Cog_sc)

# Order by FamilyNumber, necessary step for using GEE 
datatrios <- datatrios[order(datatrios$FamilyNumber),] 

# Save data trios
#write.table(datatrios, "Data_trios_NTR_20200511.csv", row.names=F, quote=F)
#datatrios <- fread("Data_trios_NTR_20200511.csv", header=T, colClasses=c("FISNumber"="character"))
#str(datatrios)

# * 3.2 Trios analyses with EA =========
summary(is.na(datatrios$EA_sc)) #data for 2207
hist(datatrios$EA_sc)

datatriosEA <- datatrios[!is.na(datatrios$EA_sc),]

# * * 3.2.1 Analyses EA with PGS from both parents pulled together -------

#With gee: might lead to issues when bootstrapping
EA_bothparents <- gee(EA_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
             sex + yob + sex*yob +
             Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, id = FamilyNumber, data=datatrios, corstr = "exchangeable")
summary(EA_bothparents)
EA_bothparents_coef <- as.data.frame(summary(EA_bothparents)$coef)
head(EA_bothparents_coef)
colnames(EA_bothparents_coef) <- c("Estimate", "naive_SE", "naive_Z", "robust_SE", "robust_Z")
EA_bothparents_coef$Pval <- 2*pnorm(-abs(EA_bothparents_coef$robust_Z))
#                             Estimate    naive_SE    naive_Z   robust_SE   robust_Z         Pval
# (Intercept)              -12.21069663 18.57714156 -0.6572969 20.44943695 -0.5971165 5.504296e-01
# SCORE.Nontrans.Cog_sc      0.08353486  0.01570860  5.3177779  0.01579574  5.2884430 1.233620e-07
# SCORE.Nontrans.NonCog_sc   0.08391856  0.01616183  5.1923932  0.01672194  5.0184689 5.208491e-07
# SCORE.Trans.Cog_sc         0.15818493  0.01581062 10.0049768  0.01552264 10.1905973 2.184180e-24
# SCORE.Trans.NonCog_sc      0.14257480  0.01571810  9.0707429  0.01606072  8.8772364 6.854205e-19


# With lme 
EA_bothparents_lme <- lme(EA_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
                            sex + yob + sex*yob + 
                            Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FamilyNumber, method="ML", na.action=na.omit, data=datatriosEA)


EA_bothparents_coef <- summary(EA_bothparents_lme)$tTable
# Value    Std.Error   DF    t-value      p-value
# (Intercept)              -11.772613610 1.854943e+01 1253 -0.6346617 5.257649e-01
# SCORE.Nontrans.Cog_sc      0.083633521 1.576822e-02  935  5.3039285 1.415182e-07
# SCORE.Nontrans.NonCog_sc   0.083983778 1.622382e-02  935  5.1765731 2.766792e-07
# SCORE.Trans.Cog_sc         0.158383632 1.587046e-02  935  9.9797772 2.341669e-22
# SCORE.Trans.NonCog_sc      0.142455775 1.577912e-02  935  9.0281204 9.762234e-19

#write.table(EA_bothparents_coef, "Estimates_Trios_NTR_EA_20200511.csv", quote=F )

total_NonCog <- EA_bothparents_coef[5,1] # total is transmitted
total_Cog <- EA_bothparents_coef[4,1]
indirect_NonCog <- EA_bothparents_coef[3,1] #indirect is nontransmitted
indirect_Cog <- EA_bothparents_coef[2,1]
direct_NonCog <- total_NonCog - indirect_NonCog #0.058472
direct_Cog <- total_Cog - indirect_Cog #0.07475011
ratio_NonCog <- indirect_NonCog/direct_NonCog  #1.436308
ratio_Cog <- indirect_Cog/direct_Cog  #1.118841
ratio_tot_NonCog <- indirect_NonCog/total_NonCog #0.5895428
ratio_tot_Cog <- indirect_Cog/total_Cog #0.528044


# * * * 3.2.1.1 Bootstrapping ------
nboot <- 10000
# gee doesn't always converge when bootstrapping
# bootcoef<-function(data,index){
#   datx<-data[index,]
#   mod<-gee(EA_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
#               sex + yob + sex*yob + 
#               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, id = FamilyNumber, data=datx, corstr = "exchangeable")  
#   mod$coefficients
#     }

bootcoef<-function(data,index){
  datx<-data[index,]
  mod <- lme(EA_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
               sex + yob + sex*yob + 
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FamilyNumber, method="ML", na.action=na.omit, data=datx, control=lmeControl(opt = "optim"))
  fixef(mod) #get fixed effects
}

boot.out<-boot(datatriosEA, bootcoef, nboot, parallel = "multicore", ncpus=20) 
boot.out

saveRDS(boot.out, "bootstrapped_output_trios_NTR_EA_20200511.Rda")

#plot to check bootstrapping
png("NTR.trios.EA.bootstrap_lme.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out)
dev.off()

#Save t output of boot
bootoutput <- as.data.frame(boot.out$t)
colnames(bootoutput) <- rownames(as.data.frame(boot.out$t0))
head(bootoutput)
#write.table(bootoutput, "Data_scores_trios_NTR_EA_bootstrapped_20200511.csv", row.names=F, quote=F)

# get values out of boot.out for all estimates + create direct and ratio estimates
original <- as.data.frame(t(boot.out$t0)) # estimates of the original sample #best estimates of the effects

original$total_NonCog <- original$SCORE.Trans.NonCog_sc
original$total_Cog <- original$SCORE.Trans.Cog_sc
original$indirect_NonCog <- original$SCORE.Nontrans.NonCog_sc
original$indirect_Cog <- original$SCORE.Nontrans.Cog_sc
original$direct_NonCog <- original$total_NonCog - original$indirect_NonCog
original$direct_Cog <- original$total_Cog - original$indirect_Cog
original$ratio_NonCog <- original$indirect_NonCog / original$direct_NonCog
original$ratio_Cog <- original$indirect_Cog / original$direct_Cog
original$ratio_tot_NonCog <- original$indirect_NonCog / original$total_NonCog
original$ratio_tot_Cog <- original$indirect_Cog / original$total_Cog

bootoutput$total_NonCog <- bootoutput$SCORE.Trans.NonCog_sc
bootoutput$total_Cog <- bootoutput$SCORE.Trans.Cog_sc
bootoutput$indirect_NonCog <- bootoutput$SCORE.Nontrans.NonCog_sc
bootoutput$indirect_Cog <- bootoutput$SCORE.Nontrans.Cog_sc
bootoutput$direct_NonCog <- bootoutput$total_NonCog - bootoutput$indirect_NonCog
bootoutput$direct_Cog <- bootoutput$total_Cog - bootoutput$indirect_Cog
bootoutput$ratio_NonCog <- bootoutput$indirect_NonCog / bootoutput$direct_NonCog
bootoutput$ratio_Cog <- bootoutput$indirect_Cog / bootoutput$direct_Cog
bootoutput$ratio_tot_NonCog <- bootoutput$indirect_NonCog / bootoutput$total_NonCog
bootoutput$ratio_tot_Cog <- bootoutput$indirect_Cog / bootoutput$total_Cog

mean <- apply(bootoutput, 2, mean) # mean of the estimates of the bootstrap resamples
bias <- mean - original
se <- apply(bootoutput, 2, sd) #the standard deviation of the bootstrap estimates is the standard error of the sample estimates

error <- qnorm(0.975)*se
leftCI <- original - bias - error # normal Ci from boot.ci 
rightCI <- original - bias + error

statsoutput <- rbind(original, mean, bias, se, error, leftCI, rightCI)
statsoutput$Estimates <- c('original', 'mean', 'bias', 'se', 'error', 'leftCI', 'rightCI')
statsoutput
tot <- statsoutput[,c(ncol(statsoutput), 1:(ncol(statsoutput)-1))]
tot

write.table(tot, "summary_mean_CI_trios_NTR_EA_20200511.csv", row.names=T, quote=F)
bmain <- bootoutput

#Perform  t-tests
t.test(bmain$indirect_NonCog, bmain$direct_NonCog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_NonCog and bmain$direct_NonCog
# t = 88.786, df = 16408, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.02622384 0.02740786
# sample estimates:
#   mean of x  mean of y 
# 0.08471725 0.05790140 

t.test(bmain$indirect_Cog, bmain$direct_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_Cog and bmain$direct_Cog
# t = 37.294, df = 16189, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.01004722 0.01116195
# sample estimates:
#   mean of x  mean of y 
# 0.08448657 0.07388199 

t.test(bmain$indirect_NonCog, bmain$indirect_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_NonCog and bmain$indirect_Cog
# t = 1.0863, df = 19881, p-value = 0.2774
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.0001855639  0.0006469130
# sample estimates:
#   mean of x  mean of y 
# 0.08471725 0.08448657 

t.test(bmain$ratio_NonCog, bmain$ratio_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$ratio_NonCog and bmain$ratio_Cog
# t = 2.4037, df = 10116, p-value = 0.01625
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.1414117 1.3913488
# sample estimates:
#   mean of x mean of y 
# 2.207711  1.441331  

t.test(bmain$ratio_tot_NonCog, bmain$ratio_tot_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$ratio_tot_NonCog and bmain$ratio_tot_Cog
# t = 32.768, df = 19155, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.06032857 0.06800511
# sample estimates:
#   mean of x mean of y 
# 0.6059586 0.5417918



# * * 3.2.2 Difference between parents -----

summary(lm(EA_sc ~ SCORE.Nontrans.Dad.Cog_sc + SCORE.Nontrans.Mom.Cog_sc + SCORE.Trans.Dad.Cog_sc + SCORE.Trans.Mom.Cog_sc + 
             SCORE.Nontrans.Dad.NonCog_sc + SCORE.Nontrans.Mom.NonCog_sc + SCORE.Trans.Dad.NonCog_sc + SCORE.Trans.Mom.NonCog_sc, data=datatrios))

# * 3.3 Trios analyses with CITO  ===========
summary(is.na(datatrios$CITO_sc)) # data for 1526

datatriosCITO <- datatrios[!is.na(datatrios$CITO_sc),]
datatriosCITO <- datatriosCITO[order(datatriosCITO$FamilyNumber),] 

# * * 3.3.1 Analyses CITO with PGS from both parents pulled together -------

#with gee #not used for the bootstrapped 
CITO_bothparents <- gee(CITO_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
                        sex + yob + sex*yob +
                        Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, id = FamilyNumber, data=datatriosCITO, corstr = "exchangeable")


summary(CITO_bothparents)
CITO_bothparents_coef <- as.data.frame(summary(CITO_bothparents)$coef)
head(CITO_bothparents_coef)
colnames(CITO_bothparents_coef) <- c("Estimate", "naive_SE", "naive_Z", "robust_SE", "robust_Z")
CITO_bothparents_coef$Pval <- 2*pnorm(-abs(CITO_bothparents_coef$robust_Z))
# Estimate    naive_SE    naive_Z   robust_SE   robust_Z         Pval
# (Intercept)              28.12470779 45.89728140  0.6127750 45.84512542  0.6134722 5.395642e-01
# SCORE.Nontrans.Cog_sc     0.01581185  0.02050492  0.7711248  0.02091467  0.7560171 4.496390e-01
# SCORE.Nontrans.NonCog_sc  0.02486621  0.02018873  1.2316879  0.02078278  1.1964818 2.315086e-01
# SCORE.Trans.Cog_sc        0.17341861  0.02096710  8.2709867  0.02109865  8.2194155 2.044973e-16
# SCORE.Trans.NonCog_sc     0.13380771  0.02021092  6.6205657  0.02061364  6.4912229 8.514240e-11

# with lme
CITO_bothparents_lme <- lme(CITO_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
                            sex + yob + sex*yob + 
                            Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FamilyNumber, method="ML", na.action=na.omit, data=datatriosCITO)

CITO_bothparents_coef <- summary(CITO_bothparents_lme)$tTable
CITO_bothparents_coef
# Value    Std.Error  DF      t-value      p-value
# (Intercept)               30.781855996  45.67037620 764  0.674000491 5.005149e-01
# SCORE.Nontrans.Cog_sc      0.015249917   0.02063496 743  0.739033069 4.601203e-01
# SCORE.Nontrans.NonCog_sc   0.023512527   0.02029745 743  1.158398112 2.470738e-01
# SCORE.Trans.Cog_sc         0.173913210   0.02108268 743  8.249104774 7.242193e-16
# SCORE.Trans.NonCog_sc      0.135007171   0.02031540 743  6.645557509 5.844153e-11

write.table(CITO_bothparents_coef, "Estimates_Trios_NTR_CITO_20200511.csv", quote=F )

# Extract estimates
total_NonCog <- CITO_bothparents_coef[5,1] # total is transmitted
total_Cog <- CITO_bothparents_coef[4,1]
indirect_NonCog <- CITO_bothparents_coef[3,1] #indirect is nontransmitted
indirect_Cog <- CITO_bothparents_coef[2,1]
direct_NonCog <- total_NonCog - indirect_NonCog 
direct_Cog <- total_Cog - indirect_Cog 
ratio_NonCog <- indirect_NonCog/direct_NonCog  
ratio_Cog <- indirect_Cog/direct_Cog  
ratio_tot_NonCog <- indirect_NonCog/total_NonCog 
ratio_tot_Cog <- indirect_Cog/total_Cog 

# * * * 3.3.1.1 Bootstrapping ------
nboot <- 10000

# Bootstrapping with gee is not really reliable, several iterations do not converge
# bootcoef<-function(data,index){
#   datx <- data[index,]
#   datx <- datx[order(datx$FamilyNumber),] 
#   mod <- gee(CITO_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
#               sex + yob + sex*yob +
#               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, id = FamilyNumber, data=datx, corstr = "exchangeable", silent=T)
#   mod$coefficients
# }

bootcoef<-function(data,index){
  datx<-data[index,]
  mod <- lme(CITO_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
               sex + yob + sex*yob + 
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FamilyNumber, method="ML", na.action=na.omit, data=datx, control=lmeControl(opt = "optim"))
  fixef(mod) #get fixed effects
}

boot.out<-boot(datatriosCITO, bootcoef, nboot, parallel = "multicore", ncpus=20) 
boot.out
saveRDS(boot.out, "bootstrapped_output_trios_NTR_CITO_20200511.Rda")

#plot to check bootstrapping
png("NTR.trios.CITO.bootstrap_lme.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out)
dev.off()

#Save t output of boot
bootoutput <- as.data.frame(boot.out$t)
colnames(bootoutput) <- rownames(as.data.frame(boot.out$t0))
head(bootoutput)
write.table(bootoutput, "Data_scores_trios_NTR_CITO_bootstrapped_20200511.csv", row.names=F, quote=F)

# get values out of boot.out for all estimates + create direct and ratio estimates
original <- as.data.frame(t(boot.out$t0)) # estimates of the original sample #best estimates of the effects

original$total_NonCog <- original$SCORE.Trans.NonCog_sc
original$total_Cog <- original$SCORE.Trans.Cog_sc
original$indirect_NonCog <- original$SCORE.Nontrans.NonCog_sc
original$indirect_Cog <- original$SCORE.Nontrans.Cog_sc
original$direct_NonCog <- original$total_NonCog - original$indirect_NonCog
original$direct_Cog <- original$total_Cog - original$indirect_Cog
original$ratio_NonCog <- original$indirect_NonCog / original$direct_NonCog
original$ratio_Cog <- original$indirect_Cog / original$direct_Cog
original$ratio_tot_NonCog <- original$indirect_NonCog / original$total_NonCog
original$ratio_tot_Cog <- original$indirect_Cog / original$total_Cog

bootoutput
bootoutput$total_NonCog <- bootoutput$SCORE.Trans.NonCog_sc
bootoutput$total_Cog <- bootoutput$SCORE.Trans.Cog_sc
bootoutput$indirect_NonCog <- bootoutput$SCORE.Nontrans.NonCog_sc
bootoutput$indirect_Cog <- bootoutput$SCORE.Nontrans.Cog_sc
bootoutput$direct_NonCog <- bootoutput$total_NonCog - bootoutput$indirect_NonCog
bootoutput$direct_Cog <- bootoutput$total_Cog - bootoutput$indirect_Cog
bootoutput$ratio_NonCog <- bootoutput$indirect_NonCog / bootoutput$direct_NonCog
bootoutput$ratio_Cog <- bootoutput$indirect_Cog / bootoutput$direct_Cog
bootoutput$ratio_tot_NonCog <- bootoutput$indirect_NonCog / bootoutput$total_NonCog
bootoutput$ratio_tot_Cog <- bootoutput$indirect_Cog / bootoutput$total_Cog

mean <- apply(bootoutput, 2, mean) # mean of the estimates of the bootstrap resamples
bias <- mean - original
se <- apply(bootoutput, 2, sd) #the standard deviation of the bootstrap estimates is the standard error of the sample estimates

error <- qnorm(0.975)*se
leftCI <- original - bias - error # normal Ci from boot.ci 
rightCI <- original - bias + error

statsoutput <- rbind(original, mean, bias, se, error, leftCI, rightCI)
statsoutput$Estimates <- c('original', 'mean', 'bias', 'se', 'error', 'leftCI', 'rightCI')
statsoutput
tot <- statsoutput[,c(ncol(statsoutput), 1:(ncol(statsoutput)-1))]
tot

write.table(tot, "summary_mean_CI_trios_NTR_CITO_20200511.csv", row.names=T, quote=F)


# * * 3.3.2 Analyses with PGS from parents separately  -------
CITO <- gee(cito_final ~ SCORE.Nontrans.Dad.Cog_sc + SCORE.Nontrans.Mom.Cog_sc + SCORE.Trans.Dad.Cog_sc + SCORE.Trans.Mom.Cog_sc + 
     SCORE.Nontrans.Dad.NonCog_sc + SCORE.Nontrans.Mom.NonCog_sc + SCORE.Trans.Dad.NonCog_sc + SCORE.Trans.Mom.NonCog_sc + 
     sex + yob + sex*yob + 
     Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, id = FamilyNumber, data=datatrios, corstr = "exchangeable")

summary(CITO)
CITO_coef <- as.data.frame(summary(CITO)$coef)
head(CITO_coef)
colnames(CITO_coef) <- c("Estimate", "naive_SE", "naive_Z", "robust_SE", "robust_Z")
CITO_coef$Pval <- 2*pnorm(-abs(CITO_coef$robust_Z))





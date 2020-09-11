
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
#setwd("C:/Users/Admin/Documents/GeneticNurtureNonCog/NTR")

# 1. Load Data ##################################

# Get phenotypic data from NTR
pheno <- fread('Phenotypic_data/NTR-DAC-1927_Demange_parental_env_EA_20200429.csv', # data was changed from spss to csv in spss
               header=T, colClasses=c("fisnumber"="character"))  # FISNumber is read as integer64 so fix to to be read as chr 
head(pheno)
str(pheno)
summary(pheno) #24984

names(pheno)[names(pheno) == 'fisnumber'] <- 'FISNumber'

# get variable and value labels
label <- read_spss('Phenotypic_data/NTR-DAC-1927_Demange_parental_env_EA_20200429.sav')
label.var <- get_label(label)
label.val <- get_labels(label,  values = "as.name")
label.val

#Load covariate file 
cov <- fread("Scores_NTR/MRG10_Pedigree_V2_NTR_only.sav/MRG10_Pedigree_V2_NTR_only.csv", 
             colClasses=c("ID"="character"))
head(cov)
str(cov)
cov <- as.data.frame(cov)



# * 1.1 Clean and adjust phenotypic data: EA + CITO ======================
# * * 1.1.1 Clean education data, recode to EduYears -----------------------------
# NTR contains two highest education variables: educat_a and educat_c. 
# Educat_a contains 4 categories: 1: primary school only, 2: lower vocational and lower
#    secondary school, 3: intermediate vocational and intermediate and higher 
#    secondary school, 4: higher vocational school and university  
# Educat_c  contains 7 categories: !: primary school only, 2: lower vocational,
#    3: lower secondary, 4: intermediate vocational, 5: intermadiate/higher secondary. 
#     6: higher vocational schooling, 7: university. 

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

# Recode educat_c
# recode educat_c in ISCED following Okbay et al 2016 Supplementary Table 1.3 
pheno$ISCED_c <- pheno$educat_c
pheno$ISCED_c[pheno$ISCED_c == 0] <- "ISCED0"
pheno$ISCED_c[pheno$ISCED_c == 1] <- "ISCED1"
pheno$ISCED_c[pheno$ISCED_c == 2] <- "ISCED2"
pheno$ISCED_c[pheno$ISCED_c == 3] <- "ISCED2"
pheno$ISCED_c[pheno$ISCED_c == 4] <- "ISCED3"
pheno$ISCED_c[pheno$ISCED_c == 5] <- "ISCED3"
pheno$ISCED_c[pheno$ISCED_c == 6] <- "ISCED5"
pheno$ISCED_c[pheno$ISCED_c == 7] <- "ISCED5"
pheno$ISCED_c[pheno$ISCED_c == 8] <- "ISCED6"
# recode in EduYears following Okbay et al 2016 Supplementary Table 1.6
pheno$Eduyears_c <-  pheno$ISCED_c
pheno$Eduyears_c[pheno$Eduyears_c == "ISCED0"] <- 1
pheno$Eduyears_c[pheno$Eduyears_c == "ISCED1"] <- 7
pheno$Eduyears_c[pheno$Eduyears_c == "ISCED2"] <- 10
pheno$Eduyears_c[pheno$Eduyears_c == "ISCED3"] <- 13
pheno$Eduyears_c[pheno$Eduyears_c == "ISCED5"] <- 19
pheno$Eduyears_c[pheno$Eduyears_c == "ISCED6"] <- 22
pheno$Eduyears_c <- as.numeric(pheno$Eduyears_c)
hist(pheno$Eduyears_c) # Only four categories: 7, 10, 13, and 19
#these 4 categories corresponds as well to the categories present in educat_a which has a higher sample size

# Recode Educat_a 
# recode educat_a in EduYears following above classification 
pheno$ISCED <- pheno$educat_a
pheno$ISCED[pheno$ISCED == 1] <- "ISCED1"
pheno$ISCED[pheno$ISCED == 2] <- "ISCED2"
pheno$ISCED[pheno$ISCED == 3] <- "ISCED3"
pheno$ISCED[pheno$ISCED == 4] <- "ISCED5"
pheno$Eduyears <- pheno$ISCED
pheno$Eduyears[pheno$Eduyears == "ISCED1"] <- 7
pheno$Eduyears[pheno$Eduyears == "ISCED2"] <- 10
pheno$Eduyears[pheno$Eduyears == "ISCED3"] <- 13
pheno$Eduyears[pheno$Eduyears == "ISCED5"] <- 19
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
onemember <- number[which(number$Freq == 1),] # identify family with only one member in this data #2003
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
# Get all multiple 
multiple <- datafamily[which(datafamily$Extension > 0 & datafamily$Extension <6),] #12002
label.val$twzyg
summary(multiple$twzyg) # 55 NAs
MZ <- multiple[which(multiple$twzyg == 1 | multiple$twzyg == 3),] #6654
DZ <- multiple[which(multiple$twzyg == 2 | multiple$twzyg > 3),] #5293
summary(DZ$multiple_type)

# Investigate NAs in zygosity
head(multiple[is.na(multiple$twzyg),])
summary(multiple[is.na(multiple$twzyg),]$famzyg) # all are MZ or DZ or NA
multiple[is.na(multiple$twzyg),][is.na(multiple[is.na(multiple$twzyg),]$famzyg)] 
# multiple with both zygosity and fam zygosity missing are in different family, except 2 triplets in fam 21783
# Get multiple with family zygosity being MZ (while zygosity is na)
MZ_multipleNa <- multiple[is.na(multiple$twzyg),][which(multiple[is.na(multiple$twzyg),]$famzyg == 1 |
                                                        multiple[is.na(multiple$twzyg),]$famzyg == 3),]
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

noncog <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_ALL_NONCOG_CEULDpred_inf.profile", 
                header=T, stringsAsFactors =F, colClasses=c("IID"="character"))
head(noncog)
cog <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_ALL_COG_CEULDpred_inf.profile", 
             header=T, stringsAsFactors =F, colClasses=c("IID"="character"))
head(cog)
head(allsib)
scores <- merge(noncog, cog, by= c("FID", "IID", "PHENO", "CNT", "CNT2"), 
                suffixes = c( ".NonCog", ".Cog"))
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
#write.table(datasib, "Data_siblings_NTR_20200531.csv", row.names=F, quote=F)
#datasib <- fread("Data_siblings_NTR_20200531.csv", colClasses=c("FISNumber"="character"))

# * 2.2 Sibling analysis with EA ==========================
# * * 2.2.1 Clean data  ------
# Keep only data with EA
datasibEA <- datasib[!is.na(datasib$Eduyears),] #3672
head(datasibEA)

number <- as.data.frame(table(datasibEA$FamilyNumber))
number
onemember <- number[which(number$Freq == 1),] #identify family with only one member in this data #695
list_onemember <- onemember$Var1 
datasibEA_sibonly<- datasibEA[!which(datasibEA$FamilyNumber %in% list_onemember), ] #2438
head(datasibEA_sibonly)

finalsib <- as.data.frame(datasibEA_sibonly)
# sample size = 3163

# * * 2.2.2 Sample descriptive -----------
fameasib <- unique(finalsib$FamilyNumber)
length(unique(finalsib$FamilyNumber)) #1309

summary(finalsib$sex)
nrow(finalsib[finalsib$sex==1,]) #1169 male
nrow(finalsib[finalsib$sex==2,]) #1994 female 
nrow(finalsib[finalsib$sex==1,]) / (nrow(finalsib[finalsib$sex==1,])+nrow(finalsib[finalsib$sex==2,]))
# 0.3695858

summary(finalsib$yob)
sd(finalsib$yob) # 12.99691
summary(finalsib$Eduyears)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.00   13.00   19.00   15.52   19.00   19.00 
sd(finalsib$Eduyears) #3.755471

summary(finalsib$SCORE.Cog)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -8.485e-07 -3.479e-07 -2.518e-07 -2.491e-07 -1.473e-07  3.244e-07 
summary(finalsib$SCORE.NonCog)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -2.917e-07 -1.372e-07 -9.382e-08 -9.349e-08 -5.049e-08  1.109e-07 
cor(finalsib$SCORE.NonCog,finalsib$SCORE.Cog) # -0.2999655

# Correlation between siblings

finalsib2 <- finalsib[order(finalsib$FamilyNumber),] 

groupfam<-group_by(finalsib,FamilyNumber) %>% summarize(m=cor())


# * * 2.2.3 Scale data  ------

finalsib[,c("EA_sc","scoreNonCog_sc", "scoreCog_sc")]<-apply(finalsib[,c("Eduyears",
                                                                         "SCORE.NonCog", 
                                                                         "SCORE.Cog")],
                                                             2,
                                                             scale)
hist(finalsib$EA_sc)

## Save data
#write.table(finalsib, "Data_siblings_NTR_EA_20200531.csv", row.names=F, quote=F)
#finalsib <- fread("Data_siblings_NTR_EA_20200531.csv", colClasses=c("FISNumber"="character"))
# head(finalsib)


# * * 2.2.4 ICC  -----
# ICC: functions written by Saskia Selzam, from Selzam et al. 2019
#calculate intraclass correlations
#i.e.  The ICC is the ratio of the between-family (i.e., random intercept) variance over the total variance 
#and is an estimate of how much of the total variation in the outcome is accounted for by family
ICCest <- function(model) {
  icc <- sqrt(diag(getVarCov(model)))^2 / (sqrt(diag(getVarCov(model)))^2 + model$sigma^2 )
  as.vector(icc)
}

# intercept model
m0 <- lme(EA_sc~1, 
          random=~1|FamilyNumber, 
          method="ML", 
          na.action=na.omit,
          data=finalsib)
ICCest(m0) #get ICC #0.4362072

m0 <- lme(SCORE.Cog~1, 
          random=~1|FamilyNumber, 
          method="ML", 
          na.action=na.omit,
          data=finalsib)
ICCest(m0) 


m0 <- lme(SCORE.NonCog~1, 
          random=~1|FamilyNumber, 
          method="ML", 
          na.action=na.omit,
          data=finalsib)
ICCest(m0) 



# * * 2.2.5 Simple linear model -----

global <- lm(EA_sc ~ scoreNonCog_sc + scoreCog_sc + sex + yob + sex*yob + 
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
             data=finalsib)
summary(global)

# * * 2.2.6 Create between-family and within-family estimates of PGS -----
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

cor.test(finalsib$GPS_W_NonCog, finalsib$GPS_B_NonCog) #-6.097644e-18 p=1
cor.test(finalsib$GPS_W_Cog, finalsib$GPS_B_Cog) #3.261838e-18  p=1


# * * 2.2.7 Run mixed model between-within regression -----
final <- lme(EA_sc~GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog + 
               sex + yob + sex*yob +
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
             random=~1|FamilyNumber, 
             method="ML", 
             na.action=na.omit,
             data=finalsib)
summary(final)


# Extract estimates 
names(summary(final))
summary(final)$tTable
# Value    Std.Error   DF    t-value      p-value
# (Intercept)  -14.238502310  8.259188949 1838 -1.7239589 8.488342e-02
# GPS_B_NonCog   0.259302293  0.023793454 1306 10.8980520 1.569952e-26
# GPS_B_Cog      0.234767971  0.023928715 1306  9.8111398 5.666971e-22
# GPS_W_NonCog   0.116504076  0.026178616 1838  4.4503528 9.087021e-06
# GPS_W_Cog      0.156129061  0.025660112 1838  6.0845044 1.418435e-09

#write.table(summary(final)$tTable, "Estimates_Siblings_NTR_EA_20200531.csv", quote=F )

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


# * * 2.2.8 Bootstrapping ------
nboot <- 10000
bootcoef<-function(data,index){
  datx<-data[index,]
  mod<-lme(EA_sc~GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog + 
               sex + yob + sex*yob +
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
           random=~1|FamilyNumber, 
           method="ML", 
           na.action=na.omit, 
           data=datx, 
           control=lmeControl(opt = "optim")) #need this optimizer to reach convergence
  fixef(mod) #get fixed effects
}

# Carry out bootstrap
boot.out<-boot(finalsib,bootcoef,nboot, parallel = "multicore", ncpus=20) 

#saveRDS(boot.out, "bootstrapped_output_siblings_NTR_EA_20200531.Rda")
#boot.out <- readRDS("bootstrapped_output_siblings_NTR_EA_20200531.Rda")

# Plot to check bootstrapping
png("NTR.sib.EA.bootstrap_20200531.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out)
dev.off()


# Save t output of boot
bootoutput <- as.data.frame(boot.out$t)
colnames(bootoutput) <- rownames(as.data.frame(boot.out$t0))
head(bootoutput)
#write.table(bootoutput, "Data_scores_siblings_NTR_EA_bootstrapped_20200531.csv", 
#              row.names=F, quote=F)

# Get values out of boot.out for all estimates + create indirect and ratio estimates
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

write.table(tot, "summary_mean_CI_siblings_NTR_EA_20200531.csv", row.names=T, quote=F)

# * * 2.2.9  Comparing estimates ------

diffcog <- original$direct_Cog - original$indirect_Cog 
diffnoncog <- original$direct_NonCog - original$indirect_NonCog
diffratio  <- original$ratio_tot_Cog - original$ratio_tot_NonCog

SD_sampling_diffcog <- sd(bootoutput$direct_Cog - bootoutput$indirect_Cog)
SD_sampling_diffnoncog <- sd(bootoutput$direct_NonCog - bootoutput$indirect_NonCog)
SD_sampling_diffratio <- sd(bootoutput$ratio_tot_Cog - bootoutput$ratio_tot_NonCog)

Z_diffcog <- diffcog/SD_sampling_diffcog
Z_diffnoncog <- diffnoncog/SD_sampling_diffnoncog
Z_diffratio <- diffratio/SD_sampling_diffratio

P_diffcog <- 2*pnorm(-abs(Z_diffcog))
P_diffnoncog <- 2*pnorm(-abs(Z_diffnoncog))
P_diffratio <- 2*pnorm(-abs(Z_diffratio))

compare <- cbind(Z_diffcog, P_diffcog, Z_diffnoncog, P_diffnoncog, Z_diffratio, P_diffratio)

#write.table(compare, "Ztests_sib_NTR_EA_20200531.csv", row.names=T, quote=F)


# * 2.3 Siblings analyses with CITO ==============================
# * * 2.3.1 Clean data and scale variables -----
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

# Sample descriptive 
famcitosib <- unique(finalsib$FamilyNumber)
length(unique(finalsib$FamilyNumber)) #757

summary(finalsib$sex)
nrow(finalsib[finalsib$sex==1,]) #734 male
nrow(finalsib[finalsib$sex==2,]) #897 female 
nrow(finalsib[finalsib$sex==1,]) / (nrow(finalsib[finalsib$sex==1,]) +nrow(finalsib[finalsib$sex==2,]))
# 0.4500307

summary(finalsib$yob)
sd(finalsib$yob) # 3.523733
summary(finalsib$cito_final)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 507.0   534.0   540.0   538.9   545.0   550.0 
sd(finalsib$cito_final) # 8.315999

summary(finalsib$SCORE.Cog)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -8.485e-07 -3.500e-07 -2.431e-07 -2.445e-07 -1.390e-07  2.482e-07 
summary(finalsib$SCORE.NonCog)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -2.856e-07 -1.358e-07 -9.383e-08 -9.493e-08 -5.272e-08  1.137e-07
cor(finalsib$SCORE.NonCog,finalsib$SCORE.Cog) # -0.2616856


# Scale variables
finalsib[,c("CITO_sc","scoreNonCog_sc", "scoreCog_sc")]<-apply(finalsib[,c("cito_final",
                                                                           "SCORE.NonCog", 
                                                                           "SCORE.Cog")],
                                                             2,
                                                             scale)
# check how data looks 
hist(finalsib$cito_final)
hist(finalsib$CITO_sc)
mean(finalsib$CITO_sc)
var(finalsib$CITO_sc)
hist(finalsib$SCORE.Cog)
hist(finalsib$scoreCog_sc)

## Save data
# write.table(finalsib, "Data_siblings_NTR_CITO_20200511.csv", row.names=F, quote=F)
#finalsib <- fread("Data_siblings_NTR_CITO_20200511.csv", colClasses=c("FISNumber"="character"))

# * * 2.3.2 ICC  -----
m0 <- lme(CITO_sc~1, 
          random=~1|FamilyNumber, 
          method="ML", 
          na.action=na.omit,
          data=finalsib)
ICCest(m0) #get ICC #0.4113066

m0 <- lme(SCORE.Cog~1, 
          random=~1|FamilyNumber, 
          method="ML", 
          na.action=na.omit,
          data=finalsib)
ICCest(m0) 


m0 <- lme(SCORE.NonCog~1, 
          random=~1|FamilyNumber, 
          method="ML", 
          na.action=na.omit,
          data=finalsib)
ICCest(m0) 

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

cor.test(finalsib$GPS_W_NonCog, finalsib$GPS_B_NonCog) #-2.231701e-18  p=1
cor.test(finalsib$GPS_W_Cog, finalsib$GPS_B_Cog) #8.14536e-18  p=1


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
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
           random=~1|FamilyNumber,
           method="ML", 
           na.action=na.omit, 
           data=datx, 
           control=lmeControl(opt = "optim"))
  fixef(mod) #get fixed effects
}

# Carry out bootstrap
boot.out<-boot(finalsib,bootcoef,nboot, parallel = "multicore", ncpus=20) 

#saveRDS(boot.out, "bootstrapped_output_siblings_NTR_CITO_20200511.Rda")
#boot.out <- readRDS("bootstrapped_output_siblings_NTR_CITO_20200511.Rda")

# Plot to check bootstrapping
png("NTR.sib.CITO.bootstrap_20200511.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out)
dev.off()

# Save t output of boot
bootoutput <- as.data.frame(boot.out$t)
colnames(bootoutput) <- rownames(as.data.frame(boot.out$t0))
head(bootoutput)
# write.table(bootoutput, "Data_scores_siblings_NTR_CITO_bootstrapped_20200511.csv", row.names=F, quote=F)

# Get values out of boot.out for all estimates + create indirect and ratio estimates
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

# write.table(tot, "summary_mean_CI_siblings_NTR_CITO_20200511.csv", row.names=T, quote=F)

# * * 2.3.6  Comparing estimates ------

diffcog <- original$direct_Cog - original$indirect_Cog 
diffnoncog <- original$direct_NonCog - original$indirect_NonCog
diffratio  <- original$ratio_tot_Cog - original$ratio_tot_NonCog

SD_sampling_diffcog <- sd(bootoutput$direct_Cog - bootoutput$indirect_Cog)
SD_sampling_diffnoncog <- sd(bootoutput$direct_NonCog - bootoutput$indirect_NonCog)
SD_sampling_diffratio <- sd(bootoutput$ratio_tot_Cog - bootoutput$ratio_tot_NonCog)

Z_diffcog <- diffcog/SD_sampling_diffcog
Z_diffnoncog <- diffnoncog/SD_sampling_diffnoncog
Z_diffratio <- diffratio/SD_sampling_diffratio

P_diffcog <- 2*pnorm(-abs(Z_diffcog))
P_diffnoncog <- 2*pnorm(-abs(Z_diffnoncog))
P_diffratio <- 2*pnorm(-abs(Z_diffratio))

compare <- cbind(Z_diffcog, P_diffcog, Z_diffnoncog, P_diffnoncog, Z_diffratio, P_diffratio)

#write.table(compare, "Ztests_sib_NTR_CITO_20200531.csv", row.names=T, quote=F)


# 3.  Trios analyses #####################################
# * 3.1 Clean data for trios analyses ===============================
# * * 3.1.1 Get transmitted and non-transmitted scores ----------
# only upload IID and score 
noncog_dad_trans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_DAD_Trans_NONCOG_CEULDpred_inf.profile", 
                          header=T, stringsAsFactors =F)[,c(2,6)]
noncog_dad_nontrans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_DAD_NonTrans_NONCOG_CEULDpred_inf.profile", 
                             header=T, stringsAsFactors =F)[,c(2,6)]
noncog_mom_trans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_MOM_Trans_NONCOG_CEULDpred_inf.profile", 
                          header=T, stringsAsFactors =F)[,c(2,6)]
noncog_mom_nontrans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_MOM_NonTrans_NONCOG_CEULDpred_inf.profile", 
                             header=T, stringsAsFactors =F)[,c(2,6)]

cog_dad_trans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_DAD_Trans_COG_CEULDpred_inf.profile", 
                       header=T, stringsAsFactors =F)[,c(2,6)]
cog_dad_nontrans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_DAD_NonTrans_COG_CEULDpred_inf.profile", 
                          header=T, stringsAsFactors =F)[,c(2,6)]
cog_mom_trans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_MOM_Trans_COG_CEULDpred_inf.profile", 
                       header=T, stringsAsFactors =F)[,c(2,6)]
cog_mom_nontrans <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_MOM_NonTrans_COG_CEULDpred_inf.profile", 
                          header=T, stringsAsFactors =F)[,c(2,6)]

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


noncog_dad <-  merge(noncog_dad_trans, noncog_dad_nontrans, 
                     by=  "IID", 
                     suffixes = c( ".Trans", ".Nontrans"))
noncog_mom <-  merge(noncog_mom_trans, noncog_mom_nontrans, 
                     by= "IID", 
                     suffixes = c( ".Trans", ".Nontrans"))
noncog <- merge(noncog_dad, noncog_mom, 
                by= "IID", 
                suffixes = c( ".Dad", ".Mom"))
head(noncog)


cog_dad <-  merge(cog_dad_trans, cog_dad_nontrans, 
                  by=  "IID", 
                  suffixes = c( ".Trans", ".Nontrans"))
cog_mom <-  merge(cog_mom_trans, cog_mom_nontrans, 
                  by= "IID", 
                  suffixes = c( ".Trans", ".Nontrans"))
cog <- merge(cog_dad, cog_mom, 
             by= "IID", 
             suffixes = c( ".Dad", ".Mom"))
head(cog)

scores<- merge(noncog, cog, 
               by="IID", 
               suffixes = c(".NonCog", ".Cog"))
str(scores)               

datatrios <- merge(pheno, scores, 
                   by.x = "FISNumber", by.y= "IID") #7448
head(datatrios)
datatrios <- as.data.frame(datatrios)

# add covariates
data <- merge(datatrios, cov, 
              by.x= "FISNumber", by.y= "ID")
head(data)
str(data)
datatrios <- data

# Reverse scores
f1 <- function(x) (x*-1)
datatrios[,c("SCORE.Nontrans.Dad.Cog", 
             "SCORE.Nontrans.Mom.Cog", 
             "SCORE.Trans.Dad.Cog", 
             "SCORE.Trans.Mom.Cog", 
             "SCORE.Nontrans.Dad.NonCog", 
             "SCORE.Nontrans.Mom.NonCog", 
             "SCORE.Trans.Dad.NonCog", 
             "SCORE.Trans.Mom.NonCog")]<-apply(datatrios[,c("SCORE.Nontrans.Dad.Cog", 
                                                            "SCORE.Nontrans.Mom.Cog", 
                                                            "SCORE.Trans.Dad.Cog", 
                                                            "SCORE.Trans.Mom.Cog", 
                                                            "SCORE.Nontrans.Dad.NonCog", 
                                                            "SCORE.Nontrans.Mom.NonCog", 
                                                            "SCORE.Trans.Dad.NonCog", 
                                                            "SCORE.Trans.Mom.NonCog")],
                                                  2,
                                                  f1)


# scale scores 
datatrios[,c("CITO_sc", 
             "EA_sc",
             "SCORE.Nontrans.Dad.Cog_sc", 
             "SCORE.Nontrans.Mom.Cog_sc", 
             "SCORE.Trans.Dad.Cog_sc", 
             "SCORE.Trans.Mom.Cog_sc", 
             "SCORE.Nontrans.Dad.NonCog_sc", 
             "SCORE.Nontrans.Mom.NonCog_sc", 
             "SCORE.Trans.Dad.NonCog_sc", 
             "SCORE.Trans.Mom.NonCog_sc")]<-apply(datatrios[,c("cito_final", 
                                                               "Eduyears",
                                                               "SCORE.Nontrans.Dad.Cog", 
                                                               "SCORE.Nontrans.Mom.Cog", 
                                                               "SCORE.Trans.Dad.Cog", 
                                                               "SCORE.Trans.Mom.Cog", 
                                                               "SCORE.Nontrans.Dad.NonCog", 
                                                               "SCORE.Nontrans.Mom.NonCog", 
                                                               "SCORE.Trans.Dad.NonCog", 
                                                               "SCORE.Trans.Mom.NonCog")],
                                                  2,
                                                  scale)

hist(datatrios$SCORE.Nontrans.Dad.Cog_sc)

# Get PGS tranmsitted and non-transmitted for both parents combined
datatrios$SCORE.Nontrans.Cog <- (datatrios$SCORE.Nontrans.Dad.Cog + datatrios$SCORE.Nontrans.Mom.Cog)/2
datatrios$SCORE.Trans.Cog <- (datatrios$SCORE.Trans.Dad.Cog + datatrios$SCORE.Trans.Mom.Cog)/2
datatrios$SCORE.Nontrans.NonCog <- (datatrios$SCORE.Nontrans.Dad.NonCog + datatrios$SCORE.Nontrans.Mom.NonCog)/2
datatrios$SCORE.Trans.NonCog <- (datatrios$SCORE.Trans.Dad.NonCog + datatrios$SCORE.Trans.Mom.NonCog)/2

datatrios$SCORE.Nontrans.Cog_sc <- scale(datatrios$SCORE.Nontrans.Cog)
datatrios$SCORE.Trans.Cog_sc <- scale(datatrios$SCORE.Trans.Cog)
datatrios$SCORE.Nontrans.NonCog_sc <-  scale(datatrios$SCORE.Nontrans.NonCog)
datatrios$SCORE.Trans.NonCog_sc <-  scale(datatrios$SCORE.Trans.NonCog)

#Investigate correlation of the transmitted and non-transmitted polygenic scores
cor.test(datatrios$SCORE.Nontrans.NonCog_sc, datatrios$SCORE.Trans.NonCog_sc) #0.02866575
cor.test(datatrios$SCORE.Nontrans.Cog_sc, datatrios$SCORE.Trans.Cog_sc)

cor.test(datatrios$SCORE.Nontrans.Dad.NonCog, datatrios$SCORE.Trans.Dad.NonCog)
cor.test(datatrios$SCORE.Nontrans.Dad.Cog, datatrios$SCORE.Trans.Dad.Cog)

cor.test(datatrios$SCORE.Nontrans.Mom.NonCog, datatrios$SCORE.Trans.Mom.NonCog)
cor.test(datatrios$SCORE.Nontrans.Mom.Cog, datatrios$SCORE.Trans.Mom.Cog)

# Order by FamilyNumber, necessary step for using GEE 
datatrios <- datatrios[order(datatrios$FamilyNumber),] 

# Save data trios
#write.table(datatrios, "Data_trios_NTR_20200531.csv", row.names=F, quote=F)
#datatrios <- fread("Data_trios_NTR_20200531.csv", header=T, colClasses=c("FISNumber"="character"))
#str(datatrios)

# * 3.2 Trios analyses with EA =========
summary(is.na(datatrios$EA_sc)) #data for 2535
hist(datatrios$EA_sc)

datatriosEA <- datatrios[!is.na(datatrios$EA_sc),]

# * * 3.2.1 Sample descriptive -------

summary(datatriosEA$sex) # there is one participant with -9 
datatriosEA <- datatriosEA[datatriosEA$sex>0, ]
#new sample size is 2534

fameatrio <- unique(datatriosEA$FamilyNumber)
length(unique(datatriosEA$FamilyNumber)) #1337

nrow(datatriosEA[datatriosEA$sex==1,]) #902 male
nrow(datatriosEA[datatriosEA$sex==2,]) #1632 female 
nrow(datatriosEA[datatriosEA$sex==1,])/ (nrow(datatriosEA[datatriosEA$sex==1,]) + nrow(datatriosEA[datatriosEA$sex==2,]))
#0.355959

summary(datatriosEA$yob)
sd(datatriosEA$yob) # 8.143513
summary(datatriosEA$Eduyears)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.00   13.00   19.00   16.16   19.00   19.00  
sd(datatriosEA$Eduyears) #3.425723

cor(datatriosEA$SCORE.Nontrans.NonCog, datatriosEA$SCORE.Nontrans.Cog) #-0.2821807
cor(datatriosEA$SCORE.Trans.NonCog, datatriosEA$SCORE.Trans.Cog)# -0.2749023

cor.test(datatriosEA$SCORE.Nontrans.NonCog, datatriosEA$SCORE.Trans.NonCog) # 0.02224433 p=0.263
cor.test(datatriosEA$SCORE.Nontrans.Cog, datatriosEA$SCORE.Trans.Cog)# 0.003770773  p=0.8495

# Get PGS for parents for correlation score between parents 

datatriosEA$SCORE.Dad.Cog <- (datatriosEA$SCORE.Nontrans.Dad.Cog 
                              + datatriosEA$SCORE.Trans.Dad.Cog)/2
datatriosEA$SCORE.Dad.NonCog <- (datatriosEA$SCORE.Nontrans.Dad.NonCog 
                                 + datatriosEA$SCORE.Trans.Dad.NonCog)/2
datatriosEA$SCORE.Mom.Cog <- (datatriosEA$SCORE.Nontrans.Mom.Cog 
                              + datatriosEA$SCORE.Trans.Mom.Cog)/2
datatriosEA$SCORE.Mom.NonCog <- (datatriosEA$SCORE.Nontrans.Mom.NonCog
                                 + datatriosEA$SCORE.Trans.Mom.NonCog)/2

cor.test(datatriosEA$SCORE.Dad.Cog, datatriosEA$SCORE.Mom.Cog) # 0.00358434 p = 0.8569
cor.test(datatriosEA$SCORE.Dad.NonCog, datatriosEA$SCORE.Mom.NonCog) # 0.03194523 p= 0.1079

# Phenotypic correlation of parental EA 
# data available for both parents is educat_c 
hist(datatriosEA$dad_educat_c)
summary(datatriosEA$dad_educat_c) #missing 238
summary(datatriosEA$mum_educat_c) #missing 273
cor.test(datatriosEA$dad_educat_c, datatriosEA$mum_educat_c) #cor=0.5964049 p< 2.2e-16

# Overlap with sibling subset
# be careful the finalsib data loaded is EA only and not CITO
overlap <- merge(finalsib, datatriosEA, by="FISNumber") #1374

# * * 3.2.2 Analyses EA with PGS from both parents pulled together -------

#With gee: might lead to issues when bootstrapping
EA_bothparents <- gee(EA_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + 
                        SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
                        sex + yob + sex*yob +
                        Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                      id = FamilyNumber, 
                      data=datatriosEA, 
                      corstr = "exchangeable")

summary(EA_bothparents)
EA_bothparents_coef <- as.data.frame(summary(EA_bothparents)$coef)
head(EA_bothparents_coef)
colnames(EA_bothparents_coef) <- c("Estimate", "naive_SE", "naive_Z", "robust_SE", "robust_Z")
EA_bothparents_coef$Pval <- 2*pnorm(-abs(EA_bothparents_coef$robust_Z))
# Estimate    naive_SE    naive_Z   robust_SE   robust_Z         Pval
# (Intercept)              -17.1803402 17.90987855 -0.9592662 19.27228798 -0.8914531 3.726862e-01
# SCORE.Nontrans.Cog_sc      0.1119145  0.02098053  5.3342052  0.02129393  5.2556976 1.474644e-07
# SCORE.Nontrans.NonCog_sc   0.1096250  0.02163807  5.0663019  0.02221917  4.9338028 8.064387e-07
# SCORE.Trans.Cog_sc         0.2178107  0.02121058 10.2689668  0.02106147 10.3416690 4.565168e-25
# SCORE.Trans.NonCog_sc      0.2080600  0.02129834  9.7688339  0.02137421  9.7341598 2.155918e-22
# sex                       -9.3119390  9.94507687 -0.9363366 10.74898606 -0.8663086 3.863210e-01


# With lme 
EA_bothparents_lme <- lme(EA_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + 
                            SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
                            sex + yob + sex*yob + 
                            Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                          random=~1|FamilyNumber, 
                          method="ML", 
                          na.action=na.omit, 
                          data=datatriosEA)


EA_bothparents_coef <- summary(EA_bothparents_lme)$tTable
# Value    Std.Error   DF     t-value      p-value
# (Intercept)              -16.737263109 1.789880e+01 1336 -0.93510538 3.499030e-01
# SCORE.Nontrans.Cog_sc      0.111772484 2.104074e-02 1179  5.31219328 1.293994e-07
# SCORE.Nontrans.NonCog_sc   0.109510120 2.169886e-02 1179  5.04681384 5.199286e-07
# SCORE.Trans.Cog_sc         0.218014744 2.127095e-02 1179 10.24941439 1.132108e-23
# SCORE.Trans.NonCog_sc      0.208268342 2.136085e-02 1179  9.75000307 1.182153e-21

#write.table(EA_bothparents_coef, "Estimates_Trios_NTR_EA_20200531.csv", quote=F )

total_NonCog <- EA_bothparents_coef[5,1] # total is transmitted
total_Cog <- EA_bothparents_coef[4,1]
indirect_NonCog <- EA_bothparents_coef[3,1] #indirect is nontransmitted
indirect_Cog <- EA_bothparents_coef[2,1]
direct_NonCog <- total_NonCog - indirect_NonCog 
direct_Cog <- total_Cog - indirect_Cog 
ratio_NonCog <- indirect_NonCog/direct_NonCog  
ratio_Cog <- indirect_Cog/direct_Cog  
ratio_tot_NonCog <- indirect_NonCog/total_NonCog 
ratio_tot_Cog <- indirect_Cog/total_Cog 


# * * * 3.2.2.1 Bootstrapping ------
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
  mod <- lme(EA_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + 
               SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
               sex + yob + sex*yob + 
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
             random=~1|FamilyNumber, 
             method="ML", 
             na.action=na.omit, 
             data=datx, 
             control=lmeControl(opt = "optim"))
  fixef(mod) #get fixed effects
}

boot.out<-boot(datatriosEA, bootcoef, nboot, parallel = "multicore", ncpus=20) 
boot.out

#saveRDS(boot.out, "bootstrapped_output_trios_NTR_EA_20200531.Rda")
#boot.out <- readRDS("bootstrapped_output_trios_NTR_EA_20200531.Rda")

# Plot to check bootstrapping
png("NTR.trios.EA.bootstrap_lme_20200531.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out)
dev.off()

# Save t output of boot
bootoutput <- as.data.frame(boot.out$t)
colnames(bootoutput) <- rownames(as.data.frame(boot.out$t0))
head(bootoutput)
#write.table(bootoutput, "Data_scores_trios_NTR_EA_bootstrapped_20200531.csv", row.names=F, quote=F)

# Get values out of boot.out for all estimates + create direct and ratio estimates
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

#write.table(tot, "summary_mean_CI_trios_NTR_EA_20200531.csv", row.names=T, quote=F)

# * * * 3.2.2.2  Comparing estimates ------

diffcog <- original$direct_Cog - original$indirect_Cog 
diffnoncog <- original$direct_NonCog - original$indirect_NonCog
diffratio  <- original$ratio_tot_Cog - original$ratio_tot_NonCog

SD_sampling_diffcog <- sd(bootoutput$direct_Cog - bootoutput$indirect_Cog)
SD_sampling_diffnoncog <- sd(bootoutput$direct_NonCog - bootoutput$indirect_NonCog)
SD_sampling_diffratio <- sd(bootoutput$ratio_tot_Cog - bootoutput$ratio_tot_NonCog)

Z_diffcog <- diffcog/SD_sampling_diffcog
Z_diffnoncog <- diffnoncog/SD_sampling_diffnoncog
Z_diffratio <- diffratio/SD_sampling_diffratio

P_diffcog <- 2*pnorm(-abs(Z_diffcog))
P_diffnoncog <- 2*pnorm(-abs(Z_diffnoncog))
P_diffratio <- 2*pnorm(-abs(Z_diffratio))

compare <- cbind(Z_diffcog, P_diffcog, Z_diffnoncog, P_diffnoncog, Z_diffratio, P_diffratio)

#write.table(compare, "Ztests_trio_NTR_EA_20200531.csv", row.names=T, quote=F)


# * * 3.2.3 Difference between parents -----

summary(lm(EA_sc ~ SCORE.Nontrans.Dad.Cog_sc + SCORE.Nontrans.Mom.Cog_sc + 
             SCORE.Trans.Dad.Cog_sc + SCORE.Trans.Mom.Cog_sc + 
             SCORE.Nontrans.Dad.NonCog_sc + SCORE.Nontrans.Mom.NonCog_sc + 
             SCORE.Trans.Dad.NonCog_sc + SCORE.Trans.Mom.NonCog_sc, 
           data=datatrios))

# * 3.3 Trios analyses with CITO  ===========
summary(is.na(datatrios$CITO_sc)) # data for 1526

datatriosCITO <- datatrios[!is.na(datatrios$CITO_sc),]
datatriosCITO <- datatriosCITO[order(datatriosCITO$FamilyNumber),] 

# * * 3.3.1 Sample descriptive --------------
summary(datatriosCITO$sex)

length(unique(datatriosCITO$FamilyNumber)) #765
nrow(datatriosCITO[datatriosCITO$sex==1,]) #674 male
nrow(datatriosCITO[datatriosCITO$sex==2,]) #852 female 
nrow(datatriosCITO[datatriosCITO$sex==1,])/ (nrow(datatriosCITO[datatriosCITO$sex==1,]) + nrow(datatriosCITO[datatriosCITO$sex==2,]))
#0.4416776 

summary(datatriosCITO$yob) 
sd(datatriosCITO$yob) # 3.61102
summary(datatriosCITO$cito_final)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 507.0   534.0   541.0   539.1   545.0   550.0 
sd(datatriosCITO$cito_final) #8.050127


cor(datatriosCITO$SCORE.Nontrans.NonCog, datatriosCITO$SCORE.Nontrans.Cog) #-0.2452256
cor(datatriosCITO$SCORE.Trans.NonCog, datatriosCITO$SCORE.Trans.Cog) # -0.26376

cor.test(datatriosCITO$SCORE.Nontrans.NonCog, datatriosCITO$SCORE.Trans.NonCog) # 0.07442533 p=0.003626
cor.test(datatriosCITO$SCORE.Nontrans.Cog, datatriosCITO$SCORE.Trans.Cog) # 0.05357893 p=0.03637


# Get PGS for parents for correlation score between parents 

datatriosCITO$SCORE.Dad.Cog <- (datatriosCITO$SCORE.Nontrans.Dad.Cog 
                              + datatriosCITO$SCORE.Trans.Dad.Cog)/2
datatriosCITO$SCORE.Dad.NonCog <- (datatriosCITO$SCORE.Nontrans.Dad.NonCog 
                                 + datatriosCITO$SCORE.Trans.Dad.NonCog)/2
datatriosCITO$SCORE.Mom.Cog <- (datatriosCITO$SCORE.Nontrans.Mom.Cog 
                              + datatriosCITO$SCORE.Trans.Mom.Cog)/2
datatriosCITO$SCORE.Mom.NonCog <- (datatriosCITO$SCORE.Nontrans.Mom.NonCog
                                 + datatriosCITO$SCORE.Trans.Mom.NonCog)/2

cor.test(datatriosCITO$SCORE.Dad.Cog, datatriosCITO$SCORE.Mom.Cog) # 0.01854552 p= 0.4691
cor.test(datatriosCITO$SCORE.Dad.NonCog, datatriosCITO$SCORE.Mom.NonCog) # 0.02997036 p= 0.242

# Phenotypic correlation of parental EA 
# data available for both parents is educat_c 
hist(datatriosCITO$dad_educat_c)
summary(datatriosCITO$dad_educat_c) #missing 115
summary(datatriosCITO$mum_educat_c) #missing 87
cor.test(datatriosCITO$dad_educat_c, datatriosCITO$mum_educat_c) #cor=0.448819  p< 2.2e-16

# Overlap with sibling subset
# be careful the finalsib data loaded is CITO only and not EA
overlap <- merge(finalsib, datatriosCITO, by="FISNumber") #823

# * * 3.3.2 Analyses CITO with PGS from both parents pulled together -------

#with gee #not used for the bootstrapped 
CITO_bothparents <- gee(CITO_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + 
                          SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
                          sex + yob + sex*yob +
                          Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                        id = FamilyNumber, 
                        data=datatriosCITO, 
                        corstr = "exchangeable")


summary(CITO_bothparents)
CITO_bothparents_coef <- as.data.frame(summary(CITO_bothparents)$coef)
head(CITO_bothparents_coef)
colnames(CITO_bothparents_coef) <- c("Estimate", "naive_SE", "naive_Z", "robust_SE", "robust_Z")
CITO_bothparents_coef$Pval <- 2*pnorm(-abs(CITO_bothparents_coef$robust_Z))

# Estimate    naive_SE    naive_Z   robust_SE   robust_Z         Pval
# (Intercept)              28.15338316 45.90337193  0.6133184 45.84667645  0.6140769 5.391645e-01
# SCORE.Nontrans.Cog_sc     0.02223799  0.02899614  0.7669292  0.02956270  0.7522313 4.519120e-01
# SCORE.Nontrans.NonCog_sc  0.03532489  0.02868183  1.2316122  0.02952472  1.1964511 2.315206e-01
# SCORE.Trans.Cog_sc        0.24648284  0.02988042  8.2489763  0.03004625  8.2034480 2.335888e-16
# SCORE.Trans.NonCog_sc     0.19129179  0.02891837  6.6148883  0.02948181  6.4884688 8.671311e-11
# sex                      -3.14529483 27.40242138 -0.1147816 28.68876163 -0.1096351 9.126988e-01

# with lme
CITO_bothparents_lme <- lme(CITO_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + 
                              SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
                              sex + yob + sex*yob + 
                              Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                            random=~1|FamilyNumber, 
                            method="ML", 
                            na.action=na.omit, 
                            data=datatriosCITO)

CITO_bothparents_coef <- summary(CITO_bothparents_lme)$tTable
CITO_bothparents_coef
# Value    Std.Error  DF     t-value      p-value
# (Intercept)               30.825940595  45.67605629 764  0.67488183 4.999550e-01
# SCORE.Nontrans.Cog_sc      0.021426282   0.02918063 743  0.73426386 4.630196e-01
# SCORE.Nontrans.NonCog_sc   0.033392411   0.02883703 743  1.15796993 2.472484e-01
# SCORE.Trans.Cog_sc         0.247163802   0.03004616 743  8.22613736 8.634499e-16
# SCORE.Trans.NonCog_sc      0.193004072   0.02906860 743  6.63960644 6.071255e-11

#write.table(CITO_bothparents_coef, "Estimates_Trios_NTR_CITO_20200531.csv", quote=F )

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

# * * * 3.3.2.1 Bootstrapping ------
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
  mod <- lme(CITO_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + 
               SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
               sex + yob + sex*yob + 
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
             random=~1|FamilyNumber, 
             method="ML", 
             na.action=na.omit, 
             data=datx, 
             control=lmeControl(opt = "optim"))
  fixef(mod) #get fixed effects
}

boot.out<-boot(datatriosCITO, bootcoef, nboot, parallel = "multicore", ncpus=20) 
boot.out
#saveRDS(boot.out, "bootstrapped_output_trios_NTR_CITO_20200531.Rda")
boot.out <- readRDS("bootstrapped_output_trios_NTR_CITO_20200531.Rda")

# Plot to check bootstrapping
png("NTR.trios.CITO.bootstrap_lme_2020531.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out)
dev.off()

# Save t output of boot
bootoutput <- as.data.frame(boot.out$t)
colnames(bootoutput) <- rownames(as.data.frame(boot.out$t0))
head(bootoutput)
#write.table(bootoutput, "Data_scores_trios_NTR_CITO_bootstrapped_20200531.csv", 
#            row.names=F, quote=F)

# Get values out of boot.out for all estimates + create direct and ratio estimates
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

#write.table(tot, "summary_mean_CI_trios_NTR_CITO_20200531.csv", row.names=T, quote=F)

# * * * 3.3.2.2  Comparing estimates ------

diffcog <- original$direct_Cog - original$indirect_Cog 
diffnoncog <- original$direct_NonCog - original$indirect_NonCog
diffratio  <- original$ratio_tot_Cog - original$ratio_tot_NonCog

SD_sampling_diffcog <- sd(bootoutput$direct_Cog - bootoutput$indirect_Cog)
SD_sampling_diffnoncog <- sd(bootoutput$direct_NonCog - bootoutput$indirect_NonCog)
SD_sampling_diffratio <- sd(bootoutput$ratio_tot_Cog - bootoutput$ratio_tot_NonCog)

Z_diffcog <- diffcog/SD_sampling_diffcog
Z_diffnoncog <- diffnoncog/SD_sampling_diffnoncog
Z_diffratio <- diffratio/SD_sampling_diffratio

P_diffcog <- 2*pnorm(-abs(Z_diffcog))
P_diffnoncog <- 2*pnorm(-abs(Z_diffnoncog))
P_diffratio <- 2*pnorm(-abs(Z_diffratio))


compare <- cbind(Z_diffcog, P_diffcog, Z_diffnoncog, P_diffnoncog, Z_diffratio, P_diffratio)
#write.table(compare, "Ztests_trio_NTR_CITO_20200531.csv", row.names=T, quote=F)


# * * 3.3.3 Analyses with PGS from parents separately  -------
CITO <- gee(cito_final ~ SCORE.Nontrans.Dad.Cog_sc + SCORE.Nontrans.Mom.Cog_sc + 
              SCORE.Trans.Dad.Cog_sc + SCORE.Trans.Mom.Cog_sc + 
              SCORE.Nontrans.Dad.NonCog_sc + SCORE.Nontrans.Mom.NonCog_sc + 
              SCORE.Trans.Dad.NonCog_sc + SCORE.Trans.Mom.NonCog_sc + 
              sex + yob + sex*yob + 
              Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
            id = FamilyNumber, 
            data=datatrios, 
            corstr = "exchangeable")

summary(CITO)
CITO_coef <- as.data.frame(summary(CITO)$coef)
head(CITO_coef)
colnames(CITO_coef) <- c("Estimate", "naive_SE", "naive_Z", "robust_SE", "robust_Z")
CITO_coef$Pval <- 2*pnorm(-abs(CITO_coef$robust_Z))

# * 3.4 Sibling effects using trios data  =========
#datatrios <- fread("Data_trios_NTR_20200531.csv", header=T, colClasses=c("FISNumber"="character"))
head(datatrios) #7448

# * * 3.4.1 Remove MZ: we can not use pairs that are only MZ because their PGS is identical ------------------------------------
# Get all multiple 
multiple <- datatrios[which(datatrios$Extension > 0 & datatrios$Extension <6),] #5758
label.val$twzyg
summary(multiple$twzyg) # 22 NAs
MZ <- multiple[which(multiple$twzyg == 1 | multiple$twzyg == 3),] #3125
DZ <- multiple[which(multiple$twzyg == 2 | multiple$twzyg > 3),] #2611
summary(DZ$multiple_type)

# to look if we want to include all sib not only dz 

# # Investigate NAs in zygosity
# # Get multiple with family zygosity being MZ (while zygosity is na)
# MZ_multipleNa <- multiple[is.na(multiple$twzyg),][which(multiple[is.na(multiple$twzyg),]$famzyg == 1 |
#                                                           multiple[is.na(multiple$twzyg),]$famzyg == 3),]
# head(MZ_multipleNa) # all of them are part of triplets # 8
# 
# MZ <- rbind(MZ, MZ_multipleNa) #6673
# 
# # Keep one MZ per pair
# 
# oneMZ <- MZ %>% group_by(FamilyNumber) %>% sample_n(1) # take at random one of MZ per family 
# oneMZ <- as.data.frame(oneMZ)
# a <- unique(oneMZ$FamilyNumber) #1600, same as oneMZ so only one per family 
# pickoneMZ <- oneMZ$FISNumber
# datawithrestMZ <- MZ[!which(MZ$FISNumber %in% pickoneMZ),] #1533
# datawithrestMZvec <- datawithrestMZ$FISNumber
# dataoneMZ <- datatrios[!which(datatrios$FISNumber %in% datawithrestMZvec),] # 5915
# 
# number <- as.data.frame(table(dataoneMZ$FamilyNumber))
# number
# onemember <- number[which(number$Freq == 1),] #identify family with only one member in this data #1375
# head(onemember)
# list_onemember <- onemember$Var1 
# dataoneMZ_sibonly<- dataoneMZ[!which(dataoneMZ$FamilyNumber %in% list_onemember), ] #4540
# head(dataoneMZ_sibonly)
# 
# datatrios_sib <- dataoneMZ_sibonly #4540
# 
# 
# # Keep only two sibling per family 
#   
# number <- as.data.frame(table(datatrios_DZ$FamilyNumber))
# 
# onemember <- number[which(number$Freq == 2),]
# as.data.frame(table(datatrios_sib$FamilyNumber))
# 
# 

# use only DZ to investigate sibling effect

datatrios_DZ <- DZ
number <- as.data.frame(table(datatrios_DZ$FamilyNumber))
twomember <- number[which(number$Freq == 2),]
list_twomember <- twomember$Var1
data2DZ <- datatrios_DZ[which(datatrios_DZ$FamilyNumber %in% list_twomember),] #2352

head(data2DZ)
summary(data2DZ$mult_ext)
First <- data2DZ[data2DZ$mult_ext == 1,]
Second <- data2DZ[data2DZ$mult_ext == 2,]

# merge with pgs of sb 
PGS_sib1 <- select(First, FamilyNumber,SCORE.Trans.Cog_sc, SCORE.Trans.NonCog_sc)
PGS_sib2 <- select(Second, FamilyNumber,SCORE.Trans.Cog_sc, SCORE.Trans.NonCog_sc)

datasib1 <- merge(First, PGS_sib2, by="FamilyNumber", suffixes= c(".sib1", ".sib2"))
datasib2 <- merge(Second, PGS_sib1, by="FamilyNumber", suffixes= c(".sib1", ".sib2"))

datatrios_sibef <- rbind(datasib1, datasib2) #2352



# * * 3.4.2 Analyses for CITO  ------------------------
datatrios_sibCITO <- datatrios_sibef[!is.na(datatrios_sibef$CITO_sc),] #N=675



CITO_sib <- lme(CITO_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + 
                              SCORE.Trans.Cog_sc.sib1 + SCORE.Trans.NonCog_sc.sib1 +
                              SCORE.Trans.Cog_sc.sib2 + SCORE.Trans.NonCog_sc.sib2 +
                              sex + yob + sex*yob + 
                              Platform + PC1 + PC2 + PC3 + 
                              PC4 + PC5 + PC6 + PC7 + PC8 +
                              PC9 + PC10, 
                            random=~1|FamilyNumber, 
                            method="ML", 
                            na.action=na.omit, 
                            data=datatrios_sibCITO)
summary(CITO_sib)


# Value Std.Error  DF   t-value p-value
# (Intercept)                  65.42895  66.46980 350  0.984341  0.3256
# SCORE.Nontrans.Cog_sc         0.05489   0.05995 304  0.915637  0.3606
# SCORE.Nontrans.NonCog_sc      0.11810   0.06136 304  1.924593  0.0552
# SCORE.Trans.Cog_sc.sib1       0.28979   0.04913 304  5.898455  0.0000
# SCORE.Trans.NonCog_sc.sib1    0.22796   0.04734 304  4.815766  0.0000
# SCORE.Trans.Cog_sc.sib2      -0.03062   0.05899 304 -0.518987  0.6041
# SCORE.Trans.NonCog_sc.sib2   -0.12360   0.05987 304 -2.064666  0.0398
# sex                         -55.27748  40.65144 304 -1.359791  0.1749
# yob                          -0.03598   0.03336 350 -1.078440  0.2816
# Platform                      0.04088   0.02270 304  1.800638  0.0728
# PC1                         -95.39372 215.21511 304 -0.443248  0.6579
# PC2                        -181.04460 119.22414 304 -1.518523  0.1299
# PC3                        -134.67633  72.72497 304 -1.851858  0.0650
# PC4                           8.19976  61.11062 304  0.134179  0.8933
# PC5                           0.15258  16.01305 304  0.009528  0.9924
# PC6                          -3.14360  30.15814 304 -0.104237  0.9170
# PC7                         -16.90077  29.78232 304 -0.567477  0.5708
# PC8                         -11.33701  25.53134 304 -0.444043  0.6573
# PC9                         -19.17351  19.03751 304 -1.007144  0.3147
# PC10                        -23.94537  14.93810 304 -1.602973  0.1100
# sex:yob                       0.02766   0.02041 304  1.355171  0.1764


# * * 3.4.2 Analyses for EA  ------------------------
datatrios_sibEA <- datatrios_sibef[!is.na(datatrios_sibef$EA_sc),] #N=788

EA_sib <- lme(EA_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + 
                  SCORE.Trans.Cog_sc.sib1 + SCORE.Trans.NonCog_sc.sib1 +
                  SCORE.Trans.Cog_sc.sib2 + SCORE.Trans.NonCog_sc.sib2 +
                  sex + yob + sex*yob + 
                  Platform + PC1 + PC2 + PC3 + 
                  PC4 + PC5 + PC6 + PC7 + PC8 +
                  PC9 + PC10, 
                random=~1|FamilyNumber, 
                method="ML", 
                na.action=na.omit, 
                data=datatrios_sibCITO)
summary(EA_sib)

# Value Std.Error  DF    t-value p-value
# (Intercept)                 -93.54892  349.5560 109 -0.2676221  0.7895
# SCORE.Nontrans.Cog_sc         0.09688    0.0973  38  0.9961848  0.3255
# SCORE.Nontrans.NonCog_sc     -0.11171    0.0933  38 -1.1972855  0.2386
# SCORE.Trans.Cog_sc.sib1       0.22343    0.1021  38  2.1893293  0.0348
# SCORE.Trans.NonCog_sc.sib1    0.25410    0.0876  38  2.9019036  0.0061
# SCORE.Trans.Cog_sc.sib2       0.00993    0.1136  38  0.0874778  0.9308
# SCORE.Trans.NonCog_sc.sib2   -0.10560    0.1025  38 -1.0298739  0.3096
# sex                         191.53139  209.0543  38  0.9161800  0.3654
# yob                           0.05069    0.1761 109  0.2879094  0.7740
# Platform                      0.02472    0.0516  38  0.4787616  0.6349
# PC1                         300.42995  362.1784  38  0.8295082  0.4120
# PC2                         290.61600  207.4233  38  1.4010770  0.1693
# PC3                        -122.73918  126.8437  38 -0.9676411  0.3393
# PC4                         -96.52209  123.1582  38 -0.7837246  0.4381
# PC5                          18.13065   25.7377  38  0.7044395  0.4855
# PC6                          50.23529   48.0467  38  1.0455505  0.3024
# PC7                         -59.20936   53.9678  38 -1.0971243  0.2795
# PC8                          77.94357   47.4549  38  1.6424757  0.1087
# PC9                          27.49711   36.4454  38  0.7544732  0.4552
# PC10                         21.05410   23.8657  38  0.8821894  0.3832
# sex:yob                      -0.09637    0.1051  38 -0.9167282  0.3651


# 4.  Siblings analyses #####################################


## Run '1. Load Data' again ##

head(pheno)


# * 4.1 Get data for MZ and DZ ===============================

# Add scores 
noncog <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_ALL_NONCOG_CEULDpred_inf.profile", 
                header=T, stringsAsFactors =F, colClasses=c("IID"="character"))

cog <- fread("Scores_NTR/ForPerlineNonAndCogInfMRG10Scores/MRG10_ALL_COG_CEULDpred_inf.profile", 
             header=T, stringsAsFactors =F, colClasses=c("IID"="character"))

scores <- merge(noncog, cog, by= c("FID", "IID", "PHENO", "CNT", "CNT2"), 
                suffixes = c( ".NonCog", ".Cog"))


datasib <- merge(pheno, scores, by.x = "FISNumber", by.y= "IID")
datasib <- merge(datasib, cov, by.x= "FISNumber", by.y= "ID")
datasib <- as.data.frame(datasib)

datasib$SCORE.Cog <- -1 * datasib$SCORE.Cog
datasib$SCORE.NonCog <- -1 * datasib$SCORE.NonCog

# Scale data
datasib[,c("EA_sc",
            "CITO_sc", 
            "scoreNonCog_sc", 
            "scoreCog_sc")]<-apply(datasib[,c("Eduyears", 
                                               "cito_final", 
                                               "SCORE.NonCog", 
                                               "SCORE.Cog")],
                                    2,
                                    scale)



# Get MZ only and Dz only 
multiple <- datasib[which(datasib$Extension > 0 & datasib$Extension <6),] #12701
MZ <- multiple[which(multiple$twzyg == 1 | multiple$twzyg == 3),] #6969
DZ <- multiple[which(multiple$twzyg == 2 | multiple$twzyg > 3),] #5657

# Get one MZ and DZ per family 

oneMZ <- MZ %>% group_by(FamilyNumber) %>% sample_n(1) # take at random one of MZ per family 
oneMZ <- as.data.frame(oneMZ)
head(oneMZ)
oneDZ <- DZ %>% group_by(FamilyNumber) %>% sample_n(1) # take at random one of DZ per family 
oneDZ <- as.data.frame(oneDZ)



# * 4.2 Analyses =====
# * * 4.2.1 EA -----------------
MZ_EA <- oneMZ[!is.na(oneMZ$EA_sc), ] #1600
DZ_EA <- oneDZ[!is.na(oneDZ$EA_sc), ] #1369


# * * 4.2.1.1 Descriptive -----------

hist(MZ_EA$Eduyears)
summary(MZ_EA$Eduyears)
summary(MZ_EA$dad_educat_c)
summary(MZ_EA$mum_educat_c)
hist(DZ_EA$Eduyears)
summary(DZ_EA$Eduyears)
summary(DZ_EA$dad_educat_c)
summary(DZ_EA$mum_educat_c)

t.test(MZ_EA$Eduyears, DZ_EA$Eduyears)
# Welch Two Sample t-test
# 
# data:  MZ_EA$Eduyears and DZ_EA$Eduyears
# t = 0.96679, df = 2913.5, p-value = 0.3337
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.1342649  0.3954446
# sample estimates:
#   mean of x mean of y 
# 15.79750  15.66691 


# * * 4.2.1.2 Effect of PGS ---------------

EA_MZlm <- lm(EA_sc ~ scoreCog_sc + scoreNonCog_sc + 
                sex + yob + sex*yob + 
                Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data=MZ_EA)
summary(EA_MZlm)


# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     16.900012  11.935505   1.416   0.1570    
# scoreCog_sc      0.195187   0.022696   8.600  < 2e-16 ***
#   scoreNonCog_sc   0.228498   0.023061   9.908  < 2e-16 ***
#   sex            -31.064434   6.665799  -4.660 3.42e-06 ***
#   yob             -0.006231   0.005967  -1.044   0.2965    
# Platform         0.020286   0.010878   1.865   0.0624 .  
# PC1            232.152563 106.424236   2.181   0.0293 *  
#   PC2             51.788403  52.025719   0.995   0.3197    
# PC3             10.070868  36.866514   0.273   0.7848    
# PC4             25.284642  32.528683   0.777   0.4371    
# PC5             11.792327   8.486487   1.390   0.1649    
# PC6            -19.678610  16.902976  -1.164   0.2445    
# PC7             25.455806  16.448196   1.548   0.1219    
# PC8            -13.492659  13.677079  -0.987   0.3240    
# PC9             23.413032  11.163926   2.097   0.0361 *  
#   PC10            -0.055084   8.122348  -0.007   0.9946    
# sex:yob          0.015670   0.003377   4.640 3.76e-06 ***
# 


EA_DZlm <- lm(EA_sc ~ scoreCog_sc + scoreNonCog_sc + 
                sex + yob + sex*yob + 
                Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
              data=DZ_EA)
summary(EA_DZlm)

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    -1.627e+01  1.351e+01  -1.204   0.2289    
# scoreCog_sc     1.920e-01  2.400e-02   7.999 2.68e-15 ***
#   scoreNonCog_sc  2.066e-01  2.461e-02   8.397  < 2e-16 ***
#   sex            -1.932e+01  7.866e+00  -2.456   0.0142 *  
#   yob             7.671e-03  6.749e-03   1.137   0.2559    
# Platform        8.783e-03  1.114e-02   0.789   0.4305    
# PC1            -1.309e+02  1.133e+02  -1.155   0.2482    
# PC2            -2.691e+01  5.044e+01  -0.534   0.5938    
# PC3             5.132e+01  3.782e+01   1.357   0.1750    
# PC4             1.208e+01  3.396e+01   0.356   0.7221    
# PC5            -1.650e+01  9.068e+00  -1.819   0.0691 .  
# PC6            -2.140e+00  1.764e+01  -0.121   0.9035    
# PC7            -1.609e+01  1.738e+01  -0.925   0.3550    
# PC8             1.289e+00  1.463e+01   0.088   0.9298    
# PC9             1.775e+01  1.124e+01   1.580   0.1144    
# PC10            1.341e+01  8.386e+00   1.599   0.1100    
# sex:yob         9.765e-03  3.986e-03   2.450   0.0144 *  


# * * 4.2.1.2 Bootstrap difference ---------------

nboot <- 10000
bootmzdzea<-function(data,index){
  datx<-data[index,]
  mod<-lm(EA_sc ~ scoreCog_sc + scoreNonCog_sc + 
            sex + yob + sex*yob + 
            Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=datx)
  return(mod$coefficients)
}


# carry out bootstrap
boot.out.mzea<-boot(MZ_EA,bootmzdzea,nboot, parallel = "multicore", ncpus=20)
boot.out.dzea<-boot(DZ_EA,bootmzdzea,nboot,parallel = "multicore", ncpus=20)
plot(boot.out.mzea)
plot(boot.out.dzea)

# Get values out of boot.out for all estimates + create difference
original.mzea <- as.data.frame(t(boot.out.mzea$t0)) # estimates of the original sample #best estimates of the effects
original.dzea <- as.data.frame(t(boot.out.dzea$t0))
bootoutput.mzea <- as.data.frame(boot.out.mzea$t)
bootoutput.dzea <- as.data.frame(boot.out.dzea$t)
colnames(bootoutput.mzea) <- rownames(as.data.frame(boot.out.mzea$t0))
colnames(bootoutput.dzea) <- rownames(as.data.frame(boot.out.dzea$t0))

diff_cog <- original.mzea$scoreCog_sc - original.dzea$scoreCog_sc
diff_noncog <- original.mzea$scoreNonCog_sc - original.dzea$scoreNonCog_sc
se_cog <- sd(bootoutput.mzea$scoreCog_sc - bootoutput.dzea$scoreCog_sc)
se_noncog <- sd(bootoutput.mzea$scoreNonCog_sc - bootoutput.dzea$scoreNonCog_sc)
z_cog <- diff_cog/se_cog
z_noncog <- diff_noncog/se_noncog
P_diffcog <- 2*pnorm(-abs(z_cog))
P_diffnoncog <- 2*pnorm(-abs(z_noncog))

compare_ea <- as.data.frame(cbind(diff_cog, se_cog, z_cog, P_diffcog, diff_noncog, se_noncog, z_noncog, P_diffnoncog))
compare_ea$pheno <- "EA"
compare_ea

# * * 4.2.2 CITO -----------------

MZ_CITO <- oneMZ[!is.na(oneMZ$CITO_sc), ] #818
DZ_CITO <- oneDZ[!is.na(oneDZ$CITO_sc), ] #865

# * * 4.2.2.1 Descriptive -----------

hist(MZ_CITO$cito_final)
summary(MZ_CITO$cito_final)
summary(MZ_CITO$dad_educat_c)
summary(MZ_CITO$mum_educat_c)
hist(DZ_CITO$cito_final)
summary(DZ_CITO$cito_final)
summary(DZ_CITO$dad_educat_c)
summary(DZ_CITO$mum_educat_c)
t.test(MZ_CITO$cito_final, DZ_CITO$cito_final)
# Welch Two Sample t-test
# 
# data:  MZ_CITO$cito_final and DZ_CITO$cito_final
# t = 0.74469, df = 1680.9, p-value = 0.4566
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.4963828  1.1040231
# sample estimates:
#   mean of x mean of y 
# 539.0171  538.7133 

# * * 4.2.2.2 Effect of PGS ---------------


CITO_MZlm <- lm(CITO_sc ~ scoreCog_sc + scoreNonCog_sc + 
                  sex + yob + sex*yob + 
                  Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                data=MZ_CITO)
summary(CITO_MZlm)

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     -70.55319   65.79847  -1.072   0.2839    
# scoreCog_sc       0.23299    0.03679   6.333 4.01e-10 ***
#   scoreNonCog_sc    0.19573    0.03523   5.555 3.78e-08 ***
#   sex              51.32309   38.79680   1.323   0.1863    
# yob               0.03622    0.03293   1.100   0.2717    
# Platform          0.03069    0.01617   1.897   0.0582 .  
# PC1            -140.44562  166.96190  -0.841   0.4005    
# PC2              79.87908   69.22207   1.154   0.2489    
# PC3             -21.97102   57.76944  -0.380   0.7038    
# PC4              19.20374   50.97728   0.377   0.7065    
# PC5               1.69520   13.81694   0.123   0.9024    
# PC6             -25.85026   24.72567  -1.045   0.2961    
# PC7             -31.36204   25.30551  -1.239   0.2156    
# PC8             -20.42639   21.37789  -0.955   0.3396    
# PC9              23.61597   16.25843   1.453   0.1467    
# PC10            -18.83579   12.13648  -1.552   0.1211    
# sex:yob          -0.02584    0.01948  -1.326   0.1851    


CITO_DZlm <- lm(CITO_sc ~ scoreCog_sc + scoreNonCog_sc + 
                  sex + yob + sex*yob + 
                  Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
                data=DZ_CITO)
summary(CITO_DZlm)

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    -131.53825   63.19349  -2.082  0.03769 *  
#   scoreCog_sc       0.32234    0.03549   9.083  < 2e-16 ***
#   scoreNonCog_sc    0.22503    0.03407   6.604 7.05e-11 ***
#   sex              74.07174   38.85744   1.906  0.05696 .  
# yob               0.06482    0.03171   2.044  0.04125 *  
#   Platform          0.04339    0.01618   2.681  0.00747 ** 
#   PC1            -117.27314  149.69781  -0.783  0.43361    
# PC2             -31.56328   73.72531  -0.428  0.66867    
# PC3             -66.09438   56.93804  -1.161  0.24604    
# PC4              -4.40009   48.22875  -0.091  0.92733    
# PC5             -10.97225   13.80965  -0.795  0.42711    
# PC6              25.77026   25.58175   1.007  0.31404    
# PC7             -19.33197   25.65405  -0.754  0.45132    
# PC8               1.31936   21.91119   0.060  0.95200    
# PC9               5.05150   16.37696   0.308  0.75781    
# PC10            -14.00959   12.34852  -1.135  0.25690    
# sex:yob          -0.03722    0.01951  -1.908  0.05675 .  

# * * 4.2.1.2 Bootstrap difference ---------------

nboot <- 10000
bootmzdzcito<-function(data,index){
  datx<-data[index,]
  mod<-lm(CITO_sc ~ scoreCog_sc + scoreNonCog_sc + 
            sex + yob + sex*yob + 
            Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=datx)
  return(mod$coefficients)
}


# carry out bootstrap
boot.out.mzcito<-boot(MZ_CITO,bootmzdzcito,nboot, parallel = "multicore", ncpus=20)
boot.out.dzcito<-boot(DZ_CITO,bootmzdzcito,nboot,parallel = "multicore", ncpus=20)
plot(boot.out.mzcito)
plot(boot.out.dzcito)

# Get values out of boot.out for all estimates + create difference
original.mzcito <- as.data.frame(t(boot.out.mzcito$t0)) # estimates of the original sample #best estimates of the effects
original.dzcito <- as.data.frame(t(boot.out.dzcito$t0))
bootoutput.mzcito <- as.data.frame(boot.out.mzcito$t)
bootoutput.dzcito <- as.data.frame(boot.out.dzcito$t)
colnames(bootoutput.mzcito) <- rownames(as.data.frame(boot.out.mzcito$t0))
colnames(bootoutput.dzcito) <- rownames(as.data.frame(boot.out.dzcito$t0))

diff_cog <- original.mzcito$scoreCog_sc - original.dzcito$scoreCog_sc
diff_noncog <- original.mzcito$scoreNonCog_sc - original.dzcito$scoreNonCog_sc
se_cog <- sd(bootoutput.mzcito$scoreCog_sc - bootoutput.dzcito$scoreCog_sc)
se_noncog <- sd(bootoutput.mzcito$scoreNonCog_sc - bootoutput.dzcito$scoreNonCog_sc)
z_cog <- diff_cog/se_cog
z_noncog <- diff_noncog/se_noncog
P_diffcog <- 2*pnorm(-abs(z_cog))
P_diffnoncog <- 2*pnorm(-abs(z_noncog))

compare_cito <- as.data.frame(cbind(diff_cog, se_cog, z_cog, P_diffcog, diff_noncog, se_noncog, z_noncog, P_diffnoncog))
compare_cito$pheno <- "CITO"
compare_cito
compare <- rbind(compare_ea, compare_cito)
#write.table(compare, "Ztests_MZDZ_sib_NTR_20200819.csv", row.names=T, quote=F)


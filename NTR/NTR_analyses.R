
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
#colnames(pheno)[colnames(pheno)=="ï..FISNumber"] <- "FISNumber" # all first columns are not weel named when exporting csv from spss in csv 

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
datasibEA <- datasib[!is.na(datasib$Eduyears),] #3144
head(datasibEA)

number <- as.data.frame(table(datasibEA$FamilyNumber))
number
onemember <- number[which(number$Freq == 1),] #identify family with only one member in this data #695
list_onemember <- onemember$Var1 
datasibEA_sibonly<- datasibEA[!which(datasibEA$FamilyNumber %in% list_onemember), ] #2449
head(datasibEA_sibonly)

finalsib <- as.data.frame(datasibEA_sibonly)
# sample size = 2449

# scale variables
finalsib[,c("EA_sc","scoreNonCog_sc", "scoreCog_sc")]<-apply(finalsib[,c("Eduyears","SCORE.NonCog", "SCORE.Cog")],
                                                             2,
                                                             scale)
hist(finalsib$EA_sc)

# save data
write.table(finalsib, "Data_siblings_NTR_EA__20200429.csv", row.names=F, quote=F)
finalsib <- fread("Data_siblings_NTR_EA__20200429.csv")
head(finalsib)

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
#                  Value Std.Error   DF   t-value p-value
# (Intercept)   -12.70919   9.12855 1392 -1.392247  0.1641
# GPS_B_NonCog    0.27790   0.02662 1038 10.440241  0.0000
# GPS_B_Cog       0.26365   0.02670 1038  9.873543  0.0000
# GPS_W_NonCog    0.11059   0.02872 1392  3.850907  0.0001
# GPS_W_Cog       0.14032   0.02814 1392  4.986136  0.0000
write.table(summary(final)$tTable, "Estimates_Siblings_NTR_EA_20200501.csv", quote=F )

direct_NonCog <- summary(final)$tTable[4,1]# direct effect is beta within 
direct_Cog <- summary(final)$tTable[5,1]
total_NonCog <- summary(final)$tTable[2,1] # total is beta between 
total_Cog <- summary(final)$tTable[3,1]
indirect_NonCog <- total_NonCog - direct_NonCog #0.1673022
indirect_Cog <- total_Cog - direct_Cog #0.1233267
ratio_NonCog <- indirect_NonCog/direct_NonCog  #1.512754
ratio_Cog <- indirect_Cog/direct_Cog  #0.8788987
ratio_tot_NonCog <- indirect_NonCog/total_NonCog #0.6020303
ratio_tot_Cog <- indirect_Cog/total_Cog #0.4677733

# * * 2.2.5 Bootstrapping ------
nboot <- 10000
bootcoef<-function(data,index){
  datx<-finalsib[index,]
  mod<-lme(EA_sc~GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog + 
               sex + yob + sex*yob +
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FamilyNumber, method="ML", na.action=na.omit, data=datx, control=lmeControl(opt = "optim"))
  fixef(mod) #get fixed effects
}

# carry out bootstrap
boot.out<-boot(finalsib,bootcoef,nboot, parallel = "multicore", ncpus=20) 

#plot to check bootstrapping
png("NTR.sib.EA.bootstrap.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out)
dev.off()

# get CI 
boot.ci(boot.out, type = c("norm", "basic"))

#Intervals : 
#Level      Normal              Basic         
#95%   (-33.09,   7.16 )   (-33.17,   7.33 )  

bootoutput <- as.data.frame(boot.out$t)
colnames(bootoutput) <- rownames(as.data.frame(boot.out$t0))
write.table(bootoutput, "Data_scores_siblings_NTR_EA_bootstrapped_2020429.csv", row.names=F, quote=F)

bmain <- bootoutput[,2:5] #save variables of interest
head(bmain)


# Create indirect and ratio variables
bmain$direct_NonCog <- bmain$GPS_W_NonCog
bmain$direct_Cog <- bmain$GPS_W_Cog
bmain$total_NonCog <- bmain$GPS_B_NonCog
bmain$total_Cog <- bmain$GPS_B_Cog
bmain$indirect_NonCog <- bmain$total_NonCog - bmain$direct_NonCog
bmain$indirect_Cog <- bmain$total_Cog - bmain$direct_Cog
bmain$ratio_NonCog <- bmain$indirect_NonCog / bmain$direct_NonCog
bmain$ratio_Cog <- bmain$indirect_Cog / bmain$direct_Cog
bmain$ratio_tot_NonCog <- bmain$indirect_NonCog / bmain$total_NonCog
bmain$ratio_tot_Cog <- bmain$indirect_Cog / bmain$total_Cog
bmain <- bmain[,5:14]
write.table(bmain, "Data_scores_siblings_NTR_EA_bootstrapped_ratio_2020429.csv", row.names=F, quote=F)

#Perform  t-tests
t.test(bmain$indirect_NonCog, bmain$direct_NonCog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_NonCog and bmain$direct_NonCog
# t = 127.92, df = 19592, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  0.06048688 0.06236936
# sample estimates:
# mean of x mean of y 
# 0.1696953 0.1082671 

t.test(bmain$indirect_Cog, bmain$direct_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_Cog and bmain$direct_Cog
# t = -23.984, df = 19669, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.01221518 -0.01036947
# sample estimates:
# mean of x mean of y 
# 0.1252292 0.1365215 

t.test(bmain$indirect_NonCog, bmain$indirect_Cog)
# 	Welch Two Sample t-test
# 
# data:  bmain$indirect_NonCog and bmain$indirect_Cog
# t = 87.704, df = 19984, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  0.04347229 0.04545982
# sample estimates:
# mean of x mean of y 
# 0.1696953 0.1252292 

t.test(bmain$ratio_NonCog, bmain$ratio_Cog)
# 	Welch Two Sample t-test
# 
# data:  bmain$ratio_NonCog and bmain$ratio_Cog
# t = 34.219, df = 11381, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  0.7836599 0.8788965
# sample estimates:
# mean of x mean of y 
#  1.869414  1.038136 

t.test(bmain$ratio_tot_NonCog, bmain$ratio_tot_Cog)
# 	Welch Two Sample t-test
# 
# data:  bmain$ratio_tot_NonCog and bmain$ratio_tot_Cog
# t = 77.776, df = 19931, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  0.1292206 0.1359021
# sample estimates:
# mean of x mean of y 
# 0.6087715 0.4762102 

# Get CI of all 
meanall <- apply(bmain, 2, mean)
sdall <- apply(bmain, 2, sd)
n <- nboot
error <- qnorm(0.975)*sdall/sqrt(n)
leftCI <- meanall-error
rightCI <- meanall+error

tot <- rbind(meanall, sdall, error, leftCI, rightCI)
#         direct_NonCog  direct_Cog total_NonCog    total_Cog indirect_NonCog indirect_Cog ratio_NonCog  ratio_Cog
# meanall  0.1082671479 0.136521536 0.2779624138 0.2617507461    0.1696952659 0.1252292101   1.86941398 1.03813577
# sdall    0.0314173948 0.031066743 0.0177187112 0.0175320353    0.0363165649 0.0353783351   2.34912323 0.61893200
# error    0.0006157696 0.000608897 0.0003472804 0.0003436216    0.0007117916 0.0006934026   0.04604197 0.01213084
# leftCI   0.1076513783 0.135912639 0.2776151334 0.2614071246    0.1689834743 0.1245358074   1.82337201 1.02600493
# rightCI  0.1088829175 0.137130433 0.2783096941 0.2620943677    0.1704070575 0.1259226127   1.91545595 1.05026662
#         ratio_tot_NonCog ratio_tot_Cog
# meanall      0.608771515   0.476210164
# sdall        0.116972301   0.123963847
# error        0.002292615   0.002429647
# leftCI       0.606478900   0.473780517
# rightCI      0.611064130   0.478639811

write.table(tot, "summary_mean_CI_siblings_NTR_EA_2020429.csv", row.names=T, quote=F)

# * 2.3 Siblings analyses with CITO ==============================
#* * 2.3.1 Clean data and scale variables -----
datasibcito <- datasib[!is.na(datasib$cito_final),] #1851

# Process again to exclude individual without sib 
number <- as.data.frame(table(datasibcito$FamilyNumber))
onemember <- number[which(number$Freq == 1),] #identify family with only one member in this data #222
head(onemember)
list <- onemember$Var1 
datasibcito <- datasibcito[!which(datasibcito$FamilyNumber %in% list), ] #1629
head(datasibcito)

finalsib <- as.data.frame(datasibcito)

summary(is.na(finalsib$cito_final))  #1629: final sample size


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
write.table(finalsib, "Data_siblings_NTR_CITO_20200430.csv", row.names=F, quote=F)
finalsib <- fread("Data_siblings_NTR_CITO_20200430.csv")

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
#                   Value Std.Error  DF   t-value p-value
# (Intercept)    16.33028  42.99540 857  0.379815  0.7042
# GPS_B_NonCog    0.17335   0.03410 753  5.083048  0.0000
# GPS_B_Cog       0.31725   0.03385 753  9.371423  0.0000
# GPS_W_NonCog    0.17727   0.03794 857  4.672696  0.0000
# GPS_W_Cog       0.23336   0.03877 857  6.019014  0.0000
write.table(summary(final)$tTable, "Estimates_Siblings_NTR_CITO_20200501.csv", quote=F )


direct_NonCog <- summary(final)$tTable[4,1]# direct effect is beta within 
direct_Cog <- summary(final)$tTable[5,1]
total_NonCog <- summary(final)$tTable[2,1] # total is beta between 
total_Cog <- summary(final)$tTable[3,1]
indirect_NonCog <- total_NonCog - direct_NonCog #- 0.003914381
indirect_Cog <- total_Cog - direct_Cog #0.08388637
ratio_NonCog <- indirect_NonCog/direct_NonCog  #-0.0220817
ratio_Cog <- indirect_Cog/direct_Cog  #0.359468
ratio_tot_NonCog <- indirect_NonCog/total_NonCog #-0.02258031
ratio_tot_Cog <- indirect_Cog/total_Cog #0.2644181

# * * 2.3.5 Bootstrapping ------------
nboot <- 10000
bootcoef<-function(data,index){
  datx<-finalsib[index,]
  mod<-lme(CITO_sc~GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog + 
               sex + yob + sex*yob +
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FamilyNumber, method="ML", na.action=na.omit, data=datx, control=lmeControl(opt = "optim"))
  fixef(mod) #get fixed effects
}

# carry out bootstrap
boot.out<-boot(finalsib,bootcoef,nboot, parallel = "multicore", ncpus=20) 


#plot to check bootstrapping
png("NTR.sib.CITO.bootstrap.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out)
dev.off()

# get CI 
boot.ci(boot.out, type = c("norm", "basic"))
# Intervals : 
# Level      Normal              Basic         
# 95%   (-84.83,  94.75 )   (-83.70,  96.74 )  
# Calculations and Intervals on Original Scale

bootoutput <- as.data.frame(boot.out$t)
colnames(bootoutput) <- rownames(as.data.frame(boot.out$t0))
write.table(bootoutput, "Data_scores_siblings_NTR_CITO_bootstrapped_20200430.csv", row.names=F, quote=F)

bmain <- bootoutput[,2:5] #save variables of interest
head(bmain)


# Create indirect and ratio variables
bmain$direct_NonCog <- bmain$GPS_W_NonCog
bmain$direct_Cog <- bmain$GPS_W_Cog
bmain$total_NonCog <- bmain$GPS_B_NonCog
bmain$total_Cog <- bmain$GPS_B_Cog
bmain$indirect_NonCog <- bmain$total_NonCog - bmain$direct_NonCog
bmain$indirect_Cog <- bmain$total_Cog - bmain$direct_Cog
bmain$ratio_NonCog <- bmain$indirect_NonCog / bmain$direct_NonCog
bmain$ratio_Cog <- bmain$indirect_Cog / bmain$direct_Cog
bmain$ratio_tot_NonCog <- bmain$indirect_NonCog / bmain$total_NonCog
bmain$ratio_tot_Cog <- bmain$indirect_Cog / bmain$total_Cog
bmain <- bmain[,5:14]
write.table(bmain, "Data_scores_siblings_NTR_CITO_bootstrapped_ratio_20200430.csv", row.names=F, quote=F)

#Perform  t-tests
t.test(bmain$indirect_NonCog, bmain$direct_NonCog)
# 	Welch Two Sample t-test
# 
# data:  bmain$indirect_NonCog and bmain$direct_NonCog
# t = -293.35, df = 19739, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.1883413 -0.1858411
# sample estimates:
#    mean of x    mean of y 
# -0.005805308  0.181285872 

t.test(bmain$indirect_Cog, bmain$direct_Cog)
# data:  bmain$indirect_Cog and bmain$direct_Cog
# t = -229.13, df = 19722, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.1491094 -0.1465799
# sample estimates:
#  mean of x  mean of y 
# 0.08715159 0.23499623 

t.test(bmain$indirect_NonCog, bmain$indirect_Cog)
# 	Welch Two Sample t-test
# 
# data:  bmain$indirect_NonCog and bmain$indirect_Cog
# t = -137.14, df = 19994, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.09428551 -0.09162830
# sample estimates:
#    mean of x    mean of y 
# -0.005805308  0.087151593 

t.test(bmain$ratio_NonCog, bmain$ratio_Cog)
# 	Welch Two Sample t-test
# 
# data:  bmain$ratio_NonCog and bmain$ratio_Cog
# t = -85.946, df = 19865, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.3984510 -0.3806821
# sample estimates:
#  mean of x  mean of y 
# 0.03212295 0.42168947 

t.test(bmain$ratio_tot_NonCog, bmain$ratio_tot_Cog)
# 	Welch Two Sample t-test
# 
# data:  bmain$ratio_tot_NonCog and bmain$ratio_tot_Cog
# t = -99.888, df = 14811, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.3233772 -0.3109301
# sample estimates:
#   mean of x   mean of y 
# -0.05012655  0.26702711 


# Get CI of all 
meanall <- apply(bmain, 2, mean)
sdall <- apply(bmain, 2, sd)
n <- nboot
error <- qnorm(0.975)*sdall/sqrt(n)
leftCI <- meanall-error
rightCI <- meanall+error

tot <- rbind(meanall, sdall, error, leftCI, rightCI)
#         direct_NonCog   direct_Cog total_NonCog    total_Cog indirect_NonCog indirect_Cog ratio_NonCog   ratio_Cog
# meanall  0.1812858722 0.2349962332 0.1754805641 0.3221478257    -0.005805308 0.0871515925   0.03212295 0.421689472
# sdall    0.0424346118 0.0428436392 0.0223781997 0.0220291433     0.047610316 0.0482473182   0.33334696 0.307135712
# error    0.0008317031 0.0008397199 0.0004386047 0.0004317633     0.000933145 0.0009456301   0.00653348 0.006019749
# leftCI   0.1804541691 0.2341565133 0.1750419594 0.3217160624    -0.006738453 0.0862059625   0.02558947 0.415669722
# rightCI  0.1821175753 0.2358359531 0.1759191687 0.3225795890    -0.004872163 0.0880972226   0.03865643 0.427709221
#         ratio_tot_NonCog ratio_tot_Cog
# meanall     -0.050126550   0.267027111
# sdall        0.283259592   0.143444329
# error        0.005551786   0.002811457
# leftCI      -0.055678336   0.264215654
# rightCI     -0.044574764   0.269838568

write.table(tot, "summary_mean_CI_siblings_NTR_CITO_20200430.csv", row.names=T, quote=F)

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
#write.table(datatrios, "Data_trios_NTR_20200430.csv", row.names=F, quote=F)
#datatrios <- fread("Data_trios_NTR_20200430.csv", header=T, colClasses=c("FISNumber"="character"))
#str(datatrios)

# * 3.2 Trios analyses with EA =========
summary(is.na(datatrios$EA_sc)) #data for 2207
hist(datatrios$EA_sc)

datatriosEA <- datatrios[!is.na(datatrios$EA_sc),]

# * * 3.2.1 Analyses with PGS from both parents pulled together -------

# With gee: might lead to issues when bootstrapping 
# EA_bothparents <- gee(EA_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
#               sex + yob + sex*yob + 
#               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, id = FamilyNumber, data=datatrios, corstr = "exchangeable")
# 
# summary(EA_bothparents)
# EA_bothparents_coef <- as.data.frame(summary(EA_bothparents)$coef)
# head(EA_bothparents_coef)
# colnames(EA_bothparents_coef) <- c("Estimate", "naive_SE", "naive_Z", "robust_SE", "robust_Z")
# EA_bothparents_coef$Pval <- 2*pnorm(-abs(EA_bothparents_coef$robust_Z))
# #                              Estimate    naive_SE    naive_Z   robust_SE   robust_Z         Pval
# # (Intercept)              -12.21069663 18.57714156 -0.6572969 20.44943695 -0.5971165 5.504296e-01
# # SCORE.Nontrans.Cog_sc      0.08353486  0.01570860  5.3177779  0.01579574  5.2884430 1.233620e-07
# # SCORE.Nontrans.NonCog_sc   0.08391856  0.01616183  5.1923932  0.01672194  5.0184689 5.208491e-07
# # SCORE.Trans.Cog_sc         0.15818493  0.01581062 10.0049768  0.01552264 10.1905973 2.184180e-24
# # SCORE.Trans.NonCog_sc      0.14257480  0.01571810  9.0707429  0.01606072  8.8772364 6.854205e-19
# 
# 
# # Extract estimates 
# total_NonCog <- EA_bothparents_coef[5,1] # total is transmitted
# total_Cog <- EA_bothparents_coef[4,1]
# indirect_NonCog <- EA_bothparents_coef[3,1] #indirect is nontransmitted
# indirect_Cog <- EA_bothparents_coef[2,1]
# direct_NonCog <- total_NonCog - indirect_NonCog #0.05865624
# direct_Cog <- total_Cog - indirect_Cog #0.07465008
# ratio_NonCog <- indirect_NonCog/direct_NonCog  #1.430684
# ratio_Cog <- indirect_Cog/direct_Cog  #1.119019
# ratio_tot_NonCog <- indirect_NonCog/total_NonCog #0.5885932
# ratio_tot_Cog <- indirect_Cog/total_Cog #0.5280835


# With lme 
EA_bothparents_lme <- lme(EA_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
                            sex + yob + sex*yob + 
                            Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FamilyNumber, method="ML", na.action=na.omit, data=datatriosEA)

# Extract estimates 

EA_bothparents_coef <- summary(EA_bothparents_lme)$tTable
# Value    Std.Error   DF    t-value      p-value
# (Intercept)              -11.772613691 1.854943e+01 1253 -0.6346617 5.257649e-01
# SCORE.Nontrans.Cog_sc      0.083633521 1.576822e-02  935  5.3039285 1.415182e-07
# SCORE.Nontrans.NonCog_sc   0.083983778 1.622382e-02  935  5.1765731 2.766792e-07
# SCORE.Trans.Cog_sc         0.158383632 1.587046e-02  935  9.9797773 2.341669e-22
# SCORE.Trans.NonCog_sc      0.142455775 1.577912e-02  935  9.0281204 9.762234e-19
# sex                      -11.402430101 1.032674e+01  935 -1.1041653 2.698055e-01
# yob                        0.006218973 9.323855e-03  935  0.6669959 5.049393e-01
# Platform                  -0.015858136 1.028448e-02  935 -1.5419481 1.234246e-01
# PC1                      -26.955376465 1.102013e+02  935 -0.2446013 8.068189e-01
# PC2                       14.859082155 5.606952e+01  935  0.2650118 7.910588e-01
# PC3                       -9.697984835 3.626914e+01  935 -0.2673894 7.892283e-01
# PC4                       28.527360194 3.172429e+01  935  0.8992278 3.687629e-01
# PC5                        3.579629973 8.501116e+00  935  0.4210777 6.737952e-01
# PC6                       -2.853055113 1.684066e+01  935 -0.1694147 8.655071e-01
# PC7                        8.900582981 1.556831e+01  935  0.5717115 5.676549e-01
# PC8                       22.418600236 1.332097e+01  935  1.6829551 9.271762e-02
# PC9                       -1.347515560 1.070564e+01  935 -0.1258697 8.998622e-01
# PC10                      -5.354564703 7.865695e+00  935 -0.6807491 4.961988e-01
# sex:yob                    0.005743111 5.222603e-03  935  1.0996644 2.717614e-01

write.table(EA_bothparents_coef, "Estimates_Trios_NTR_EA_20200501.csv", quote=F )

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
# bootcoef<-function(data,index){
#   datx<-datatrios[index,]
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

boot.out<-boot(datatriosEA,bootcoef,nboot, parallel = "multicore", ncpus=20) 
head(boot.out)


#plot to check bootstrapping
png("NTR.trios.EA.bootstrap_lme.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out)
dev.off()

# get CI 
boot.ci(boot.out, type = c("norm", "basic"))
# Intervals with gee: 
# Level      Normal              Basic         
# 95%   (-43.23,  33.70 )   (-43.72,  33.49 )  

# Intervals : 
#   Level      Normal              Basic         
# 95%   (-61.86,  23.97 )   (-61.50,  23.54 )  

bootoutput <- as.data.frame(boot.out$t)
colnames(bootoutput) <- rownames(as.data.frame(boot.out$t0))
#write.table(bootoutput, "Data_scores_trios_NTR_EA_bootstrapped_20200430.csv", row.names=F, quote=F)
write.table(bootoutput, "Data_scores_trios_NTR_EA_bootstrapped_20200501.csv", row.names=F, quote=F)


bmain <- bootoutput[,2:5] #save variables of interest
head(bmain)

# Create direct and ratio variables
bmain$total_NonCog <- bmain$SCORE.Trans.NonCog_sc
bmain$total_Cog <- bmain$SCORE.Trans.Cog_sc
bmain$indirect_NonCog <- bmain$SCORE.Nontrans.NonCog_sc
bmain$indirect_Cog <- bmain$SCORE.Nontrans.Cog_sc
bmain$direct_NonCog <- bmain$total_NonCog - bmain$indirect_NonCog
bmain$direct_Cog <- bmain$total_Cog - bmain$indirect_Cog
bmain$ratio_NonCog <- bmain$indirect_NonCog / bmain$direct_NonCog
bmain$ratio_Cog <- bmain$indirect_Cog / bmain$direct_Cog
bmain$ratio_tot_NonCog <- bmain$indirect_NonCog / bmain$total_NonCog
bmain$ratio_tot_Cog <- bmain$indirect_Cog / bmain$total_Cog
bmain <- bmain[,5:14]

#write.table(bmain, "Data_scores_trios_NTR_EA_bootstrapped_ratio_20200430.csv", row.names=F, quote=F)
write.table(bmain, "Data_scores_trios_NTR_EA_bootstrapped_ratio_20200501.csv", row.names=F, quote=F)

# Get CI of all 
meanall <- apply(bmain, 2, mean)
sdall <- apply(bmain, 2, sd)
n <- nboot
error <- qnorm(0.975)*sdall/sqrt(n)
leftCI <- meanall-error
rightCI <- meanall+error

tot <- rbind(meanall, sdall, error, leftCI, rightCI)
# total_NonCog   total_Cog indirect_NonCog indirect_Cog direct_NonCog   direct_Cog ratio_NonCog ratio_Cog ratio_tot_NonCog ratio_tot_Cog
# meanall 0.1427658660 0.158462921    0.0848940274 0.0845116680   0.057871839 0.0739512529     2.326123 1.4675715      0.606917266   0.541910532
# sdall   0.0152555048 0.015046040    0.0154951208 0.0145911894   0.026029816 0.0248108570    56.223616 3.5681933      0.153429948   0.125143465
# error   0.0002990024 0.000294897    0.0003036988 0.0002859821   0.000510175 0.0004862839     1.101963 0.0699353      0.003007172   0.002452767
# leftCI  0.1424668637 0.158168024    0.0845903286 0.0842256860   0.057361664 0.0734649691     1.224160 1.3976362      0.603910094   0.539457765
# rightCI 0.1430648684 0.158757818    0.0851977262 0.0847976501   0.058382014 0.0744375368     3.428085 1.5375068      0.609924438   0.544363299

#write.table(tot, "summary_mean_CI_trios_NTR_EA_20200430.csv", row.names=T, quote=F)
write.table(tot, "summary_mean_CI_trios_NTR_EA_20200501.csv", row.names=T, quote=F)


# * * 3.2.2 Difference between parents -----

summary(lm(EA_sc ~ SCORE.Nontrans.Dad.Cog_sc + SCORE.Nontrans.Mom.Cog_sc + SCORE.Trans.Dad.Cog_sc + SCORE.Trans.Mom.Cog_sc + 
             SCORE.Nontrans.Dad.NonCog_sc + SCORE.Nontrans.Mom.NonCog_sc + SCORE.Trans.Dad.NonCog_sc + SCORE.Trans.Mom.NonCog_sc, data=datatrios))

# * 3.3 Trios analyses with CITO  ===========
summary(is.na(datatrios$CITO_sc)) # data for 1526

datatriosCITO <- datatrios[!is.na(datatrios$CITO_sc),]

# * * 3.3.1 Analyses with PGS from both parents pulled together -------
#with gee
# CITO_bothparents <- gee(CITO_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
#                         sex + yob + sex*yob + 
#                         Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, id = FamilyNumber, data=datatriosCITO, corstr = "exchangeable")
# 
# 
# summary(CITO_bothparents)
# CITO_bothparents_coef <- as.data.frame(summary(CITO_bothparents)$coef)
# head(CITO_bothparents_coef)
# colnames(CITO_bothparents_coef) <- c("Estimate", "naive_SE", "naive_Z", "robust_SE", "robust_Z")
# CITO_bothparents_coef$Pval <- 2*pnorm(-abs(CITO_bothparents_coef$robust_Z))
# 
# #                             Estimate    naive_SE    naive_Z   robust_SE   robust_Z         Pval
# # (Intercept)              28.12470779 45.89728140  0.6127750 45.84512542  0.6134722 5.395642e-01
# # SCORE.Nontrans.Cog_sc     0.01581185  0.02050492  0.7711248  0.02091467  0.7560171 4.496390e-01
# # SCORE.Nontrans.NonCog_sc  0.02486621  0.02018873  1.2316879  0.02078278  1.1964818 2.315086e-01
# # SCORE.Trans.Cog_sc        0.17341861  0.02096710  8.2709867  0.02109865  8.2194155 2.044973e-16
# # SCORE.Trans.NonCog_sc     0.13380771  0.02021092  6.6205657  0.02061364  6.4912229 8.514240e-1
# 
# # Extract estimates 
# total_NonCog <- CITO_bothparents_coef[5,1] # total is transmitted
# total_Cog <- CITO_bothparents_coef[4,1]
# indirect_NonCog <- CITO_bothparents_coef[3,1] #indirect is nontransmitted
# indirect_Cog <- CITO_bothparents_coef[2,1]
# direct_NonCog <- total_NonCog - indirect_NonCog #0.1089415
# direct_Cog <- total_Cog - indirect_Cog #0.1576068
# ratio_NonCog <- indirect_NonCog/direct_NonCog  #0.2282529
# ratio_Cog <- indirect_Cog/direct_Cog  #0.1003247
# ratio_tot_NonCog <- indirect_NonCog/total_NonCog #0.1858354
# ratio_tot_Cog <- indirect_Cog/total_Cog #0.0911773

# wiht lme
CITO_bothparents_lme <- lme(CITO_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
                            sex + yob + sex*yob + 
                            Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FamilyNumber, method="ML", na.action=na.omit, data=datatrios)

CITO_bothparents_coef <- summary(CITO_bothparents_lme)$tTable
# Value    Std.Error  DF     t-value      p-value
# (Intercept)               30.781869851  45.67037662 764  0.67400079 5.005147e-01
# SCORE.Nontrans.Cog_sc      0.015249914   0.02063496 743  0.73903288 4.601204e-01
# SCORE.Nontrans.NonCog_sc   0.023512520   0.02029745 743  1.15839770 2.470740e-01
# SCORE.Trans.Cog_sc         0.173913213   0.02108268 743  8.24910436 7.242216e-16
# SCORE.Trans.NonCog_sc      0.135007177   0.02031540 743  6.64555739 5.844158e-11
# sex                       -3.992045812  27.18243599 743 -0.14686122 8.832814e-01
# yob                       -0.016807625   0.02291623 743 -0.73343770 4.635229e-01
# Platform                   0.036139464   0.01395023 743  2.59059945 9.768755e-03
# PC1                       -1.403532215 143.54479456 743 -0.00977766 9.922013e-01
# PC2                      -84.914487985  78.41541985 743 -1.08287998 2.792130e-01
# PC3                      -97.696503748  47.71903599 743 -2.04732769 4.097711e-02
# PC4                        4.494302461  40.03968844 743  0.11224619 9.106585e-01
# PC5                        9.043699488  11.16136715 743  0.81026808 4.180456e-01
# PC6                       -4.281199074  20.09808240 743 -0.21301530 8.313735e-01
# PC7                       -6.687229311  20.68168983 743 -0.32334057 7.465283e-01
# PC8                        1.509112893  16.81104173 743  0.08976915 9.284949e-01
# PC9                       11.074650686  12.90338483 743  0.85827485 3.910175e-01
# PC10                     -19.811338149   9.80530222 743 -2.02047196 4.369276e-02
# sex:yob                    0.001933724   0.01364669 743  0.14169916 8.873561e-01

write.table(CITO_bothparents_coef, "Estimates_Trios_NTR_CITO_20200501.csv", quote=F )

# Extract estimates
total_NonCog <- CITO_bothparents_coef[5,1] # total is transmitted
total_Cog <- CITO_bothparents_coef[4,1]
indirect_NonCog <- CITO_bothparents_coef[3,1] #indirect is nontransmitted
indirect_Cog <- CITO_bothparents_coef[2,1]
direct_NonCog <- total_NonCog - indirect_NonCog #0.1114947
direct_Cog <- total_Cog - indirect_Cog #0.1586633
ratio_NonCog <- indirect_NonCog/direct_NonCog  #0.2108847
ratio_Cog <- indirect_Cog/direct_Cog  #0.09611494
ratio_tot_NonCog <- indirect_NonCog/total_NonCog #0.1741576
ratio_tot_Cog <- indirect_Cog/total_Cog #0.08768692

# * * * 3.3.1.1 Bootstrapping ------
nboot <- 10000
# bootcoef<-function(data,index){
#   datx<-datatrios[index,]
#   mod<-gee(CITO_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
#               sex + yob + sex*yob + 
#               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, id = FamilyNumber, data=datx, corstr = "exchangeable")  
#   mod$coefficients
# }

bootcoef<-function(data,index){
  datx<-data[index,]
  mod <- lme(CITO_sc ~ SCORE.Nontrans.Cog_sc + SCORE.Nontrans.NonCog_sc + SCORE.Trans.Cog_sc + SCORE.Trans.NonCog_sc +
               sex + yob + sex*yob + 
               Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FamilyNumber, method="ML", na.action=na.omit, data=datx, control=lmeControl(opt = "optim"))
  fixef(mod) #get fixed effects
}

boot.out<-boot(datatriosCITO,bootcoef,nboot, parallel = "multicore", ncpus=20) 
head(boot.out)


#plot to check bootstrapping
png("NTR.trios.CITO.bootstrap_lme.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out)
dev.off()

# get CI 
boot.ci(boot.out, type = c("norm", "basic"))
# Intervals : 
#   Level      Normal              Basic         
# 95%   (-98.78, 113.32 )   (-98.07, 112.90 )  
# Calculations and Intervals on Original Scale

bootoutput <- as.data.frame(boot.out$t)
colnames(bootoutput) <- rownames(as.data.frame(boot.out$t0))
head(bootoutput)
write.table(bootoutput, "Data_scores_trios_NTR_CITO_bootstrapped_20200501.csv", row.names=F, quote=F)

bmain <- bootoutput[,2:5] #save variables of interest
head(bmain)

# Create direct and ratio variables
bmain$total_NonCog <- bmain$SCORE.Trans.NonCog_sc
bmain$total_Cog <- bmain$SCORE.Trans.Cog_sc
bmain$indirect_NonCog <- bmain$SCORE.Nontrans.NonCog_sc
bmain$indirect_Cog <- bmain$SCORE.Nontrans.Cog_sc
bmain$direct_NonCog <- bmain$total_NonCog - bmain$indirect_NonCog
bmain$direct_Cog <- bmain$total_Cog - bmain$indirect_Cog
bmain$ratio_NonCog <- bmain$indirect_NonCog / bmain$direct_NonCog
bmain$ratio_Cog <- bmain$indirect_Cog / bmain$direct_Cog
bmain$ratio_tot_NonCog <- bmain$indirect_NonCog / bmain$total_NonCog
bmain$ratio_tot_Cog <- bmain$indirect_Cog / bmain$total_Cog
bmain <- bmain[,5:14]

write.table(bmain, "Data_scores_trios_NTR_CITO_bootstrapped_ratio_20200501.csv", row.names=F, quote=F)
#bmain <- fread("Data_scores_trios_NTR_CITO_bootstrapped_ratio_20200501.csv")

# Get CI of all 
meanall <- apply(bmain, 2, mean)
sdall <- apply(bmain, 2, sd)
n <- nboot
error <- qnorm(0.975)*sdall/sqrt(n)
leftCI <- meanall-error
rightCI <- meanall+error

tot <- rbind(meanall, sdall, error, leftCI, rightCI)
# total_NonCog    total_Cog indirect_NonCog indirect_Cog direct_NonCog   direct_Cog ratio_NonCog   ratio_Cog ratio_tot_NonCog ratio_tot_Cog
# meanall 0.1429005593 0.1766437892    0.0166190641 0.0149600877  0.1262814952 0.1616837015  0.185701120 0.121179949      0.127761312   0.092056308
# sdall   0.0190448502 0.0201216026    0.0188627591 0.0187305407  0.0328073088 0.0335937158  0.306075655 0.162733235      0.146985320   0.114428568
# error   0.0003732722 0.0003943762    0.0003697033 0.0003671119  0.0006430114 0.0006584247  0.005998973 0.003189513      0.002880859   0.002242759
# leftCI  0.1425272871 0.1762494130    0.0162493608 0.0145929758  0.1256384837 0.1610252768  0.179702147 0.117990436      0.124880452   0.089813549
# rightCI 0.1432738315 0.1770381654    0.0169887674 0.0153271995  0.1269245066 0.1623421262  0.191700093 0.124369462      0.130642171   0.094299067

write.table(tot, "summary_mean_CI_trios_NTR_CITO_20200501.csv", row.names=T, quote=F)


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





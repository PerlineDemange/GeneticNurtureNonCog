# Rosa Cheesman & Perline Demange 
# PRS analysis of adoptees and non-adopted in UKBiobank
# 31-03-2020

###########################################
# Make file with FID, polygenic scores and covariates 
###########################################
module load pre2019
module load R/3.4.3

R
library(data.table)
library(psych)

## data of adoptees 
######################

adop<- fread("../../EA/adoptees.csv", data.table=F, header=T) 
adop <- adop[,c(1,3)]
colnames(adop)<-c("ID1","EA")

polyNC <- fread('/home/pdemange/UKB/PGS/NonCog/scores/NONCOG_LDpred-inf_scores.profile', data.table=F, header=T) 
polyNC <- polyNC[,c(2,6)]
colnames(polyNC) <- c('ID1', 'scoreNonCog')
finaladop1 <- merge(adop, polyNC, by='ID1') 
polyC <- fread('/home/pdemange/UKB/PGS/Cog/scores/COG_LDpred-inf_scores.profile', data.table=F, header=T)
polyC <- polyC[,c(2,6)]
colnames(polyC) <- c('ID1', 'scoreCog')
finaladop <- merge(finaladop1, polyC, by='ID1') 

# score were based on wrong effect alelles so sign need to be change
finaladop$scoreNonCogrev <- -1*finaladop$scoreNonCog
finaladop$scoreCogrev <- -1*finaladop$scoreCog

#scale variables 
finaladop[,c("EA_sc","scoreNonCog_sc", "scoreCog_sc")]<-apply(finaladop[,c("EA","scoreNonCogrev", "scoreCogrev")],
                                                             2,
                                                             scale)
# add covariates
sex <- read.table("/project/ukbaumc/UKBGWAS/phenotypes/sex.array.cov") 
colnames(sex) <- c("ID1", "IID", "sex", "array")

age <- read.table("/project/ukbaumc/UKBGWAS/phenotypes/age.25PCs.qcov") 
colnames(age) <- c("ID1", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13",
                   "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25", "Age")

finaladop <- merge(finaladop, sex, by='ID1')
finaladop <- merge(finaladop, age, by=c("ID1", "IID")) #6407

## data of non-adopted controls
##################################

control <- fread("../../EA/nonadopted_6.5k.csv", data.table=F, header=T) 
control <- control[,c(1,3)]
colnames(control)<-c("ID1","EA")

finalcontrol <- merge(control, polyNC, by='ID1') 
finalcontrol <- merge(finalcontrol, polyC, by='ID1') 

finalcontrol$scoreNonCogrev <- -1*finalcontrol$scoreNonCog
finalcontrol$scoreCogrev <- -1*finalcontrol$scoreCog

finalcontrol[,c("EA_sc","scoreNonCog_sc", "scoreCog_sc")]<-apply(finalcontrol[,c("EA","scoreNonCogrev", "scoreCogrev")],
                                                              2,
                                                              scale)

finalcontrol <- merge(finalcontrol, sex, by='ID1')
finalcontrol <- merge(finalcontrol, age, by=c("ID1", "IID")) #6500


####################
## Run PGS analyses
####################
library(boot)

# create function to be able to bootstrap
# DO WE USE IT? NEED TO CHANGE?
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
} 

## For adoptees
###############

# Run simple lm 
simple <- lm(EA_sc~scoreNonCog_sc + scoreCog_sc + sex + array + Age + sex*Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=finaladop)
summary(simple)

# bootstrapping with 1000 replications 
results <- boot(data=finaladop, statistic=rsq, 
                R=1000, formula=EA_sc~scoreNonCog_sc + scoreCog_sc + sex + array + Age + sex*Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)

# view results
results 

## For Non-adoptees 
###################

# run simple model 
simplecontrol <- lm(EA_sc~scoreNonCog_sc + scoreCog_sc + sex + array + Age + sex*Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=finalcontrol)
summary(simplecontrol)

# bootstrapping with 1000 replications 
results <- boot(data=finalcontrol, statistic=rsq, 
                R=1000, formula=Edu~PRS_std)
results 


## Comparison adoptees - non-adopted
#####################################
names(summary(simple))

direct_NonCog <- summary(simple)$coefficients[2,1] #beta from adopted 
direct_Cog <- summary(simple)$coefficients[3,1]
total_NonCog <- summary(simplecontrol)$coefficients[2,1] #beta from non-adopted 
total_Cog <- summary(simplecontrol)$coefficients[3,1] 
indirect_NonCog <- total_NonCog - direct_NonCog #0.02133141
indirect_Cog <- total_Cog - direct_Cog #0.07886815
ratio_NonCog <- indirect_NonCog / direct_NonCog #0.1137116
ratio_Cog <- indirect_Cog / direct_Cog #0.4330012

# Save data
#####################
write.table(finaladop, file="Data_scores_adop_UKB_20200310.csv", row.names=F, quote=F) 
write.table(finalcontrol, file="Data_scores_nonadop_UKB_20200310.csv", row.names=F, quote=F) 

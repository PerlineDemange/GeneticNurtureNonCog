# Rosa Cheesman & Perline Demange 
# PRS analysis of adoptees and non-adopted in UKBiobank
# 31-03-2020

###########################################
# Make file with FID, polygenic scores and covariates 
###########################################
module load 2019
module load R/3.5.1-foss-2019b

R
library(data.table)
library(psych)
set.seed(42)

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


## Sample descriptive
##################################

summary(finaladop$sex)
nrow(finaladop[finaladop$sex==1,]) #3057 male
nrow(finalsib[finalsib$sex==0,]) #3350 female 
nrow(finaladop[finaladop$sex==1,]) / (nrow(finaladop[finaladop$sex==1,]) + nrow(finaladop[finaladop$sex==0,])) #0.47713

summary(finaladop$Age)
sd(finaladop$Age) #8.512444

summary(finaladop$EA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 7.0    10.0    13.0    13.5    20.0    20.0

sd(finaladop$EA) #5.0796

summary(finaladop$scoreNonCogrev)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -3.450e-07 -1.420e-07 -1.001e-07 -1.003e-07 -5.839e-08  1.374e-07

summary(finaladop$scoreCogrev)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -8.309e-07 -3.894e-07 -2.845e-07 -2.840e-07 -1.772e-07  2.223e-07

cor(finaladop$scoreNonCogrev,finaladop$scoreCogrev) # -0.2858899


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

## Sample descriptive
##################################

summary(finalcontrol$sex)
nrow(finalcontrol[finalcontrol$sex==1,]) #2961 male
nrow(finalcontrol[finalcontrol$sex==0,]) #3539 female 
nrow(finalcontrol[finalcontrol$sex==1,]) / (nrow(finalcontrol[finalcontrol$sex==1,]) + nrow(finalcontrol[finalcontrol$sex==0,])) #0.455538

summary(finalcontrol$Age)
sd(finalcontrol$Age) #7.851899

summary(finalcontrol$EA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 7.00   10.00   13.00   13.93   20.00   20.00

sd(finalcontrol$EA) #5.126729

summary(finalcontrol$scoreNonCogrev)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -3.507e-07 -1.370e-07 -9.544e-08 -9.419e-08 -5.166e-08  1.611e-07

summary(finalcontrol$scoreCogrev)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -8.228e-07 -3.656e-07 -2.626e-07 -2.617e-07 -1.585e-07  2.849e-07

cor(finalcontrol$scoreNonCogrev,finalcontrol$scoreCogrev) #-0.2625942

# Save data
#####################
write.table(finaladop, file="Data_scores_adop_UKB_20200529.csv", row.names=F, quote=F) 
write.table(finalcontrol, file="Data_scores_nonadop_UKB_20200529.csv", row.names=F, quote=F) 

####################
## Run PGS analyses
####################

## For adoptees
###############

# Run simple lm 
simple <- lm(EA_sc~scoreNonCog_sc + scoreCog_sc + sex + array + Age + sex*Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=finaladop)
summary(simple)


## For Non-adoptees 
###################

# run simple model 
simplecontrol <- lm(EA_sc~scoreNonCog_sc + scoreCog_sc + sex + array + Age + sex*Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=finalcontrol)
summary(simplecontrol)


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

# save results 
adopcoef <- summary(simple)$coefficients
nonadopcoef <- summary(simplecontrol)$coefficients
colnames(adopcoef) <- paste(colnames(adopcoef), "adop", sep = "_")
colnames(nonadopcoef) <- paste(colnames(nonadopcoef), "nonadop", sep = "_")
resultsadoption <- cbind(adopcoef, nonadopcoef)
#write.table(resultsadoption, file="Results_lm_adoption_UKB_20200401.csv", row.names=T, quote=F) #results are the same with data creted on the 0330 and0529

#################
# Bootstrap 
#################

library(boot)
nboot <- 10000
bootadop<-function(data,index){
  datx<-data[index,]
  mod<-lm(EA_sc~scoreNonCog_sc + scoreCog_sc + sex + array + Age + sex*Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=datx)
  return(mod$coefficients)
}
bootcontrol<-function(data,index){
  datx<-data[index,]
  mod<-lm(EA_sc~scoreNonCog_sc + scoreCog_sc + sex + array + Age + sex*Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=datx)
  return(mod$coefficients) #get fixed effects
}

# carry out bootstrap
boot.out.adop<-boot(finaladop,bootadop,nboot, parallel = "multicore", ncpus=20)
boot.out.control<-boot(finalcontrol,bootcontrol,nboot,parallel = "multicore", ncpus=20)

#saveRDS(boot.out.adop, "bootstrapped_output_adop_UKB_EA_20200529.Rda")
#saveRDS(boot.out.control, "bootstrapped_output_nonadop_UKB_EA_20200529.Rda")

#plot to check bootstrapping
options(bitmapType='cairo')
png("UKB.adop.bootstrap_20200529.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out.adop)
dev.off()

png("UKB.nonadop.bootstrap_)20200529.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out.control)
dev.off()


# save results for each bootstrap in data frame
bootoutput.adop <- as.data.frame(boot.out.adop$t)
bootoutput.nonadop <- as.data.frame(boot.out.control$t)
colnames(bootoutput.adop) <- rownames(as.data.frame(boot.out.adop$t0))
colnames(bootoutput.nonadop) <- rownames(as.data.frame(boot.out.control$t0))
#write.table(bootoutput.adop, "Data_scores_adop_UKB_bootstrapped_20200529.csv", row.names=F, quote=F)
#write.table(bootoutput.nonadop, "Data_scores_nonadop_UKB_bootstrapped_20200529.csv", row.names=F, quote=F)

bmain.adop <- bootoutput.adop[,2:3] #Save variables of interest
bmain.nonadop <- bootoutput.nonadop[,2:3] #Save variables of interest

# Create indirect and ratio variables
bmain <- bmain.adop

bmain$direct_NonCog <- bmain.adop$scoreNonCog_sc 
bmain$direct_Cog <- bmain.adop$scoreCog_sc 
bmain$total_NonCog <- bmain.nonadop$scoreNonCog_sc 
bmain$total_Cog <- bmain.nonadop$scoreCog_sc 
bmain$indirect_NonCog <- bmain$total_NonCog - bmain$direct_NonCog
bmain$indirect_Cog <- bmain$total_Cog - bmain$direct_Cog
bmain$ratio_NonCog <- bmain$indirect_NonCog / bmain$direct_NonCog
bmain$ratio_Cog <- bmain$indirect_Cog / bmain$direct_Cog
bmain$ratio_tot_NonCog <- bmain$indirect_NonCog / bmain$total_NonCog
bmain$ratio_tot_Cog <- bmain$indirect_Cog / bmain$total_Cog

bmain <- bmain[, 3:12]


# Get values out of boot.out for all estimates + create indirect and ratio estimates
original.adop <- as.data.frame(t(boot.out.adop$t0)) # estimates of the original sample #best estimates of the effects
original.nonadop <- as.data.frame(t(boot.out.control$t0))

direct_NonCog <- original.adop$scoreNonCog_sc 
direct_Cog <- original.adop$scoreCog_sc 
total_NonCog <- original.nonadop$scoreNonCog_sc 
total_Cog <- original.nonadop$scoreCog_sc
original <- as.data.frame(cbind(direct_NonCog, direct_Cog, total_NonCog, total_Cog))
original$indirect_NonCog <- original$total_NonCog - original$direct_NonCog
original$indirect_Cog <- original$total_Cog - original$direct_Cog
original$ratio_NonCog <- original$indirect_NonCog / original$direct_NonCog
original$ratio_Cog <- original$indirect_Cog / original$direct_Cog
original$ratio_tot_NonCog <- original$indirect_NonCog / original$total_NonCog
original$ratio_tot_Cog <- original$indirect_Cog / original$total_Cog

direct_NonCog <- bootoutput.adop$scoreNonCog_sc 
direct_Cog <- bootoutput.adop$scoreCog_sc 
total_NonCog <- bootoutput.nonadop$scoreNonCog_sc 
total_Cog <- bootoutput.nonadop$scoreCog_sc 
bootoutput <- as.data.frame(cbind(direct_NonCog, direct_Cog, total_NonCog, total_Cog))
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

#write.table(tot, "summary_mean_CI_adoption_UKB_20200529.csv", row.names=T, quote=F)


### Compare estimates 
######################


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

#write.table(compare, "Ztests_adoption_UKB_20200529.csv", row.names=T, quote=F)

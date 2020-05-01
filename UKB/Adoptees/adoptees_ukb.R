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

# Save data
#####################
write.table(finaladop, file="Data_scores_adop_UKB_20200310.csv", row.names=F, quote=F) 
write.table(finalcontrol, file="Data_scores_nonadop_UKB_20200310.csv", row.names=F, quote=F) 

# save results 
adopcoef <- summary(simple)$coefficients
nonadopcoef <- summary(simplecontrol)$coefficients
colnames(adopcoef) <- paste(colnames(adopcoef), "adop", sep = "_")
colnames(nonadopcoef) <- paste(colnames(nonadopcoef), "nonadop", sep = "_")
resultsadoption <- cbind(adopcoef, nonadopcoef)
write.table(resultsadoption, file="Results_lm_adoption_UKB_20200401.csv", row.names=T, quote=F) 

#################
# Bootstrap 
#################

library(boot)
nboot <- 10000
bootadop<-function(data,index){
  datx<-finaladop[index,]
  mod<-lm(EA_sc~scoreNonCog_sc + scoreCog_sc + sex + array + Age + sex*Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=datx)
  return(mod$coefficients)
}
bootcontrol<-function(data,index){
  datx<-finalcontrol[index,]
  mod<-lm(EA_sc~scoreNonCog_sc + scoreCog_sc + sex + array + Age + sex*Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=datx)
  return(mod$coefficients) #get fixed effects
}

# carry out bootstrap
boot.out.adop<-boot(finaladop,bootadop,nboot, parallel = "multicore", ncpus=20)
boot.out.control<-boot(finalcontrol,bootcontrol,nboot,parallel = "multicore", ncpus=20)

#plot to check bootstrapping
png("UKB.adop.bootstrap.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out.adop)
dev.off()

png("UKB.nonadop.bootstrap.png",
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
write.table(bootoutput.adop, "Data_scores_adop_UKB_bootstrapped_20200401.csv", row.names=F, quote=F)
write.table(bootoutput.nonadop, "Data_scores_nonadop_UKB_bootstrapped_20200401.csv", row.names=F, quote=F)
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


# perform t.tests

t.test(bmain$indirect_NonCog, bmain$direct_NonCog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_NonCog and bmain$direct_NonCog
# t = -778.98, df = 18162, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.1666354 -0.1657989
# sample estimates:
#   mean of x  mean of y
# 0.02142506 0.18764224


t.test(bmain$indirect_Cog, bmain$direct_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_Cog and bmain$direct_Cog
# t = -497.17, df = 18263, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.1039558 -0.1031393
# sample estimates:
#   mean of x  mean of y
# 0.07874064 0.18228820


t.test(bmain$indirect_NonCog, bmain$indirect_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_NonCog and bmain$indirect_Cog
# t = -237.22, df = 19982, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.05778915 -0.05684200
# sample estimates:
#   mean of x  mean of y
# 0.02142506 0.07874064


t.test(bmain$ratio_NonCog, bmain$ratio_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$ratio_NonCog and bmain$ratio_Cog
# t = -207.89, df = 19507, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.3222513 -0.3162313
# sample estimates:
#   mean of x mean of y
# 0.1191830 0.4384243

t.test(bmain$ratio_tot_NonCog, bmain$ratio_tot_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$ratio_tot_NonCog and bmain$ratio_tot_Cog
# t = -206.17, df = 18000, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.2026915 -0.1988738
# sample estimates:
#   mean of x  mean of y
# 0.09947362 0.30025630


# Get CI of all 
meanall <- apply(bmain, 2, mean)
sdall <- apply(bmain, 2, sd)
n <- nboot
error <- qnorm(0.975)*sdall/sqrt(n)
leftCI <- meanall-error
rightCI <- meanall+error

tot <- rbind(meanall, sdall, error, leftCI, rightCI)
# direct_NonCog  direct_Cog total_NonCog    total_Cog indirect_NonCog
# meanall  0.1876422423 0.182288201 0.2090673041 0.2610288382    0.0214250618
# sdall    0.0124611240 0.012248745 0.0120412232 0.0118393837    0.0173210753
# error    0.0002442335 0.000240071 0.0002360036 0.0002320477    0.0003394868
# leftCI   0.1873980088 0.182048130 0.2088313005 0.2607967905    0.0210855749
# rightCI  0.1878864759 0.182528272 0.2093033077 0.2612608858    0.0217645486
# indirect_Cog ratio_NonCog   ratio_Cog ratio_tot_NonCog ratio_tot_Cog
# meanall 0.0787406367  0.119183009 0.438424293      0.099473624   0.300256298
# sdall   0.0168446013  0.099604507 0.116882129      0.079509067   0.056233639
# error   0.0003301481  0.001952212 0.002290848      0.001558349   0.001102159
# leftCI  0.0784104886  0.117230797 0.436133445      0.097915275   0.299154139
# rightCI 0.0790707848  0.121135222 0.440715141      0.101031973   0.301358457

write.table(tot, "summary_mean_CI_adoption_UKB_20200430.csv", row.names=T, quote=F)


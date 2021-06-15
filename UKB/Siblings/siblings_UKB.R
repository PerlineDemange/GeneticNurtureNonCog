# Rosa Cheesman & Perline Demange 
# PRS analysis of sib comparison in UKBiobank
# 24.02.20 - 30-03-2020 - 19-05-2021

# Get interactive node on LISA
#ssh pdemange@login4.lisa.surfsara.nl
#srun -N 1 -t 24:00:00 --pty bash

#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////
####################
# Make sibling id
####################
# take dataset with just sibs
# this greps loads of different files each containing occurrences of each ID from Sib_id_uniques:
awk '{print $1,$2}' Sib > Sib_ids

awk '{print $1}' Sib_ids > Sib1

awk '!seen[$1]++' Sib1 > Sib_id_uniques #only for id 1!!


mkdir rels
#for every unique id1, print out the rows this id appears in
# id1fam4
# 3 4
# 4 5
# 6 4
while read list; do grep $list Sib_ids > rels/id1fam${list}; done < Sib_id_uniques

#edit the contents of each family file resulting from above, so that ids are sorted, non duplicated, all in one column
#split up two cols from above
while read list; do awk '{print $1}' rels/id1fam${list} > rels/id1fam${list}_1; done < Sib_id_uniques
while read list; do awk '{print $2}' rels/id1fam${list} > rels/id1fam${list}_2; done < Sib_id_uniques

#append 1 and 2 and remove duplicates --this makes sure none of the pairs are duplicated--then sort
while read list; do cat rels/id1fam${list}_1 rels/id1fam${list}_2 > rels/id1fam${list}_3; done < Sib_id_uniques
while read list; do awk '!seen[$1]++' rels/id1fam${list}_3 > rels/id1fam${list}_4; done < Sib_id_uniques
while read list; do sort -k 1 rels/id1fam${list}_4 > rels/id1fam${list}_5; done < Sib_id_uniques

#now create string variable containing ids of everyone in family
#this doesnt work, it just prints the ID of the first person in that family
while read list; do awk '{ for (i=1;i<=NF;i++ ) printf $i "_" }' rels/id1fam${list}_5 > rels/famid${list}; done < Sib_id_uniques
while read list; do id=`cat rels/famid${list} | awk '{print $1}'`
awk '{print $1,'$id'}' rels/id1fam${list}_5 > rels/id1fam${list}_6; done < Sib_id_uniques

cd rels
cat id1fam*_6 > sibs_famid #should be id fid
sort -k 2 sibs_famid > sibs_famid2
awk '!seen[$1]++' sibs_famid2 > sibs_famid3

awk '{print $2}' sibs_famid3 > fams #this will be the no. of families
awk -F '\t' '{print $1}' fams | sort | uniq -c | sort -nr > no_of_fams_and_fids #do this to know what max.no of sibs is


###########################################
# Make file with FID, polygenic scores and covariates 
###########################################
module load 2019
module load R/3.5.1-foss-2019b

R
library(data.table)
library(psych)
library(nlme)
library(boot)
library(dplyr)

set.seed(42)

b <- fread("../rels/sibs_famid3",h=F, data.table=F, verbose=T) #41498
colnames(b)[1:2]<-c("ID1","FID")
##siblings
pheno <- fread("../../EA/siblings.csv", data.table=F, verbose=T, header=T) #39500
pheno <- pheno[,c(1,3)]
colnames(pheno)<-c("ID1","EA")
dat2<-merge(pheno, b ,by = "ID1") #39500
polyNC <- fread('/home/pdemange/UKB/PGS/NonCog/scores/NONCOG_LDpred-inf_scores.profile', data.table=F, header=T) 
polyNC <- polyNC[,c(2,6)]
colnames(polyNC) <- c('ID1', 'scoreNonCog')
finalsib1 <- merge(dat2, polyNC, by='ID1') #39500
polyC <- fread('/home/pdemange/UKB/PGS/Cog/scores/COG_LDpred-inf_scores.profile', data.table=F, header=T)
polyC <- polyC[,c(2,6)]
colnames(polyC) <- c('ID1', 'scoreCog')
finalsib <- merge(finalsib1, polyC, by='ID1') #39500

# add covariates

sex <- read.table("/project/ukbaumc/UKBGWAS/phenotypes/sex.array.cov") 
colnames(sex) <- c("ID1", "IID", "sex", "array")

age <- read.table("/project/ukbaumc/UKBGWAS/phenotypes/age.25PCs.qcov") 
colnames(age) <- c("ID1", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13",
                   "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25", "Age")

finalsib <- merge(finalsib, sex, by='ID1')
finalsib <- merge(finalsib, age, by=c("ID1", "IID"))


# Reverse scores because wrong effect allele was taken in LDPred
finalsib$scoreNoNCogrev <- -1 * finalsib$scoreNonCog
finalsib$scoreCogrev <- -1 * finalsib$scoreCog


###########################################
# Scale variables and average per family 
###########################################

# scale variables
finalsib[,c("EA_sc","scoreNonCog_sc", "scoreCog_sc")]<-apply(finalsib[,c("EA","scoreNoNCogrev", "scoreCogrev")],
                          2,
                          scale)


nrow(finalsib) # 39500

# Save data
#####################
#write.table(finalsib, file="Data_scores_siblings_UKB_20200529.csv", row.names=F, quote=F) 
finalsib <- read.table("Data_scores_siblings_UKB_20200529.csv", header=T, stringsAsFactors = F)
## Sample descriptive
##################################
colnames(finalsib)

length(unique(finalsib$FID)) #19389

summary(finalsib$sex)
nrow(finalsib[finalsib$sex==1,]) #16627 male
nrow(finalsib[finalsib$sex==0,]) #22873 female 
nrow(finalsib[finalsib$sex==1,]) / (nrow(finalsib[finalsib$sex==1,]) + nrow(finalsib[finalsib$sex==0,])) #0.4209367

summary(finalsib$Age)
sd(finalsib$Age) #7.338735

summary(finalsib$EA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 7.00   10.00   13.00   13.62   20.00   20.00

sd(finalsib$EA) #5.1000257

summary(finalsib$scoreNoNCogrev)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -3.656e-07 -1.368e-07 -9.468e-08 -9.421e-08 -5.193e-08  1.740e-07

summary(finalsib$scoreCogrev)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -9.590e-07 -3.631e-07 -2.582e-07 -2.570e-07 -1.513e-07  4.107e-07

cor(finalsib$scoreNoNCogrev,finalsib$scoreCogrev) # -0.2613102

############################################
# Mixed model between-within family
############################################

# linear model
global <- lm(EA_sc ~ scoreNonCog_sc + scoreCog_sc + sex + array + Age + sex*Age + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=finalsib)
summary(global)


# ICC: functions written by Saskia Selzam, from Selzam et al. 2019
#################

#calculate intraclass correlations
#i.e.  The ICC is the ratio of the between-family (i.e., random intercept) variance over the total variance 
#and is an estimate of how much of the total variation in the outcome is accounted for by family
ICCest <- function(model) {
  icc <- sqrt(diag(getVarCov(model)))^2 / (sqrt(diag(getVarCov(model)))^2 + model$sigma^2 )
  as.vector(icc)
}

# intercept model
m0 <- lme(EA_sc~1, random=~1|FID, method="ML", na.action=na.omit,data=finalsib)
ICCest(m0) #get ICC #0.3258229

m0 <- lme(scoreNonCog_sc~1, 
          random=~1|FID, 
          method="ML", 
          na.action=na.omit,
          data=finalsib, control=lmeControl(opt = "optim"))
ICCest(m0) # 0.5255507  #when using standardized measure, found 0.5255509


m0 <- lme(scoreCog_sc~1, 
          random=~1|FID, 
          method="ML", 
          na.action=na.omit,
          data=finalsib) #convergence issue when not standardized, use standardized 
ICCest(m0) #0.5292896

# Carry out bootstrap
nboot <- 1000

bootcoef_EA<-function(data,index){
  datx<-data[index,]
  mod <- lme(EA_sc~1, 
             random=~1|FID, 
             method="ML", 
             na.action=na.omit,
             data=datx, 
             control=lmeControl(opt = "optim")) #problem with convergence otherwise
  ICCest(mod) #get icc 
}

bootcoef_scoreNonCog<-function(data,index){
  datx<-data[index,]
  mod <- lme(scoreNonCog_sc~1, 
             random=~1|FID, 
             method="ML", 
             na.action=na.omit,
             data=datx, 
             control=lmeControl(opt = "optim")) #problem with convergence otherwise
  ICCest(mod) #get icc 
}

bootcoef_scoreCog<-function(data,index){
  datx<-data[index,]
  mod <- lme(scoreCog_sc~1, 
             random=~1|FID, 
             method="ML", 
             na.action=na.omit,
             data=datx, 
             control=lmeControl(opt = "optim")) #problem with convergence otherwise
  ICCest(mod) #get icc 
}

## try other optimizers

boot.ea <- boot(finalsib,bootcoef_EA,nboot, parallel = "multicore", ncpus=14) 
boot.ea$t0 # 0.3258231  estimate is different than previous from 5th decimal because of the control optim
sd(boot.ea$t) #0.004590459 the standard deviation of the bootstrap estimates is the standard error of the sample estimates
boot.noncog <- boot(finalsib,bootcoef_scoreNonCog,nboot, parallel = "multicore", ncpus=14) 
boot.noncog$t0 #0.5256068
sd(boot.noncog$t) # 0.003528444 
boot.cog <- boot(finalsib,bootcoef_scoreCog,nboot, parallel = "multicore", ncpus=14) 
boot.cog$t0 #0.5293449
sd(boot.cog$t) #0.003449969




# Create between family estimates of PGS: average across family member
######################################################################


meanNC<-group_by(finalsib,FID) %>% summarize(m=mean(scoreNonCog_sc))
colnames(meanNC) <- c("FID", "GPS_B_NonCog")
meanC<-group_by(finalsib,FID) %>% summarize(m=mean(scoreCog_sc)) #19389 rows
colnames(meanC) <- c("FID", "GPS_B_Cog")
finalsib<-merge(finalsib,meanNC,by="FID")
finalsib<-merge(finalsib,meanC,by="FID")

# Create within-family variable
################################

finalsib$GPS_W_NonCog <- finalsib$scoreNonCog_sc  - finalsib$GPS_B_NonCog  
finalsib$GPS_W_Cog <- finalsib$scoreCog_sc  - finalsib$GPS_B_Cog  

cor.test(finalsib$GPS_W_NonCog, finalsib$GPS_B_NonCog) #-7.066594e-19 p=1
cor.test(finalsib$GPS_W_Cog, finalsib$GPS_B_Cog) #7.562188e-19  p=1


# 1. Analyses with  indirect = bteween - within ##################
######################################################################

# Mixed effects model
######################

final <- lme(EA_sc~GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog +
               sex + array + Age + sex*Age + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
             random=~1|FID, method="ML", na.action=na.omit,data=finalsib)
summary(final)
# Value Std.Error    DF   t-value p-value
# (Intercept)   1.0804198 0.0520973 20095  20.73848  0.0000
# GPS_B_NonCog  0.2674314 0.0064073 19386  41.73850  0.0000
# GPS_B_Cog     0.2813860 0.0063180 19386  44.53717  0.0000
# GPS_W_NonCog  0.1101560 0.0088315 20095  12.47315  0.0000
# GPS_W_Cog     0.1234973 0.0088651 20095  13.93070  0.0000

final <- lmer(EA_sc~GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog +
               sex + array + Age + sex*Age + 
               PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                (1|FID), na.action=na.omit,data=finalsib)
# Estimate Std. Error t value
# (Intercept)   1.0804035  0.0520984  20.738
# GPS_B_NonCog  0.2674312  0.0064079  41.735
# GPS_B_Cog     0.2813859  0.0063186  44.533
# GPS_W_NonCog  0.1101561  0.0088307  12.474
# GPS_W_Cog     0.1234973  0.0088644  13.932
# 

final <- lm(EA_sc~GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog +
                sex + array + Age + sex*Age + 
                PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, na.action=na.omit,data=finalsib)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   1.0983530  0.0507654  21.636  < 2e-16 ***
#   GPS_B_NonCog  0.2678471  0.0056923  47.054  < 2e-16 ***
#   GPS_B_Cog     0.2816461  0.0056110  50.195  < 2e-16 ***
#   GPS_W_NonCog  0.1099937  0.0102424  10.739  < 2e-16 ***
#   GPS_W_Cog     0.1234489  0.0102824  12.006  < 2e-16 ***
#   sex          -0.4082619  0.0749141  -5.450 5.07e-08 ***
  


# 2. Analyses with  indirect = pop - within ##################
######################################################################

nboot <- 10000

both_lm<-function(data,index){
  datx<-data[index,]
  modpop <- lm(EA_sc ~ scoreNonCog_sc + scoreCog_sc +
              sex + array + Age + sex*Age + 
              PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, na.action=na.omit, data=datx)
  modwithin <-lm(EA_sc ~ GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog +
              sex + array + Age + sex*Age + 
              PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, na.action=na.omit, data=datx)
  results <- c(coef(modpop), coef(modwithin))
  results
}

both_lm_boot<-boot(finalsib, both_lm, nboot, parallel = "multicore", ncpus=4) 
saveRDS(both_lm_boot, "bootstrapped_output_UKB_siblings_both_lm_20210519.Rda")

options(bitmapType='cairo')
png("UKB.sib.bootstrap_both_lm_202105019.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(both_lm_boot)
dev.off()


bootoutput <- as.data.frame(both_lm_boot$t)
colnames(bootoutput) <- colnames(as.data.frame(t(both_lm_boot$t0)))

# Get values out of boot.out for all estimates + create indirect and ratio estimates
original <- as.data.frame(t(both_lm_boot$t0))
original$direct_NonCog <- original$GPS_W_NonCog
original$direct_Cog <- original$GPS_W_Cog
original$pop_NonCog <- original$scoreNonCog_sc
original$pop_Cog <- original$scoreCog_sc
original$indirect_NonCog <- original$pop_NonCog - original$direct_NonCog
original$indirect_Cog <- original$pop_Cog - original$direct_Cog
original$ratio_NonCog <- original$indirect_NonCog / original$direct_NonCog
original$ratio_Cog <- original$indirect_Cog / original$direct_Cog
original$ratio_pop_NonCog <- original$indirect_NonCog / original$pop_NonCog
original$ratio_pop_Cog <- original$indirect_Cog / original$pop_Cog

bootoutput$direct_NonCog <- bootoutput$GPS_W_NonCog
bootoutput$direct_Cog <- bootoutput$GPS_W_Cog
bootoutput$pop_NonCog <- bootoutput$scoreNonCog_sc
bootoutput$pop_Cog <- bootoutput$scoreCog_sc
bootoutput$indirect_NonCog <- bootoutput$pop_NonCog - bootoutput$direct_NonCog
bootoutput$indirect_Cog <- bootoutput$pop_Cog - bootoutput$direct_Cog
bootoutput$ratio_NonCog <- bootoutput$indirect_NonCog / bootoutput$direct_NonCog
bootoutput$ratio_Cog <- bootoutput$indirect_Cog / bootoutput$direct_Cog
bootoutput$ratio_pop_NonCog <- bootoutput$indirect_NonCog / bootoutput$pop_NonCog
bootoutput$ratio_pop_Cog <- bootoutput$indirect_Cog / bootoutput$pop_Cog


mean <- apply(bootoutput, 2, mean) # mean of the estimates of the bootstrap resamples
bias <- mean - original
se <- apply(bootoutput, 2, sd) #the standard deviation of the bootstrap estimates is the standard error of the sample estimates

error <- qnorm(0.975)*se
leftCI <- original - bias - error # normal Ci from boot.ci 
rightCI <- original - bias + error

statsoutput <- rbind(original, mean, bias, se, error, leftCI, rightCI)
statsoutput$Estimates <- c('original', 'mean', 'bias', 'se', 'error', 'leftCI', 'rightCI')
tot <- statsoutput[,(ncol(statsoutput)-10):ncol(statsoutput)] # get only summary statistics 
tot <- tot[,c(ncol(tot), 1:(ncol(tot)-1))]
tot

write.table(tot, "summary_mean_CI_siblings_UKB_pop_lm_20210519_BOTH.csv", row.names=T, quote=F)
write.table(statsoutput, "full_summary_mean_CI_siblings_UKB_pop_lm_20210519_BOTH.csv", row.names=T, quote=F)

### Compare estimates 
######################

diffcog <- original$direct_Cog - original$indirect_Cog 
diffnoncog <- original$direct_NonCog - original$indirect_NonCog
diffratio  <- original$ratio_pop_Cog - original$ratio_pop_NonCog

SD_sampling_diffcog <- sd(bootoutput$direct_Cog - bootoutput$indirect_Cog)
SD_sampling_diffnoncog <- sd(bootoutput$direct_NonCog - bootoutput$indirect_NonCog)
SD_sampling_diffratio <- sd(bootoutput$ratio_pop_Cog - bootoutput$ratio_pop_NonCog)

Z_diffcog <- diffcog/SD_sampling_diffcog
Z_diffnoncog <- diffnoncog/SD_sampling_diffnoncog
Z_diffratio <- diffratio/SD_sampling_diffratio

P_diffcog <- 2*pnorm(-abs(Z_diffcog))
P_diffnoncog <- 2*pnorm(-abs(Z_diffnoncog))
P_diffratio <- 2*pnorm(-abs(Z_diffratio))

compare <- cbind(Z_diffcog, P_diffcog, Z_diffnoncog, P_diffnoncog, Z_diffratio, P_diffratio)

write.table(compare, "Ztests_sib_UKB_pop_lm_20210519.csv", row.names=T, quote=F)

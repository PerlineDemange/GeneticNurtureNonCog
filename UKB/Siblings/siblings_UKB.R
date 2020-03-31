# Rosa Cheesman & Perline Demange 
# PRS analysis of sib comparison in UKBiobank
# 24.02.20 - 30-03-2020

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
module load pre2019
module load R/3.4.3

R
library(data.table)
library(psych)

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




############################################
# Mixed model between-within family
############################################

# linear model
global <- lm(EA_sc ~ scoreNonCog_sc + scoreCog_sc + sex + array + Age + sex*Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=finalsib)
summary(global)


# libraries
library(nlme)
library(boot)


# ICC and total effect: functions written by Saskia Selzam, from Selzam et al. 2019
#################

#calculate intraclass correlations
#i.e.  The ICC is the ratio of the between-family (i.e., random intercept) variance over the total variance 
#and is an estimate of how much of the total variation in the outcome is accounted for by family
ICCest <- function(model) {
  icc <- sqrt(diag(getVarCov(model)))^2 / (sqrt(diag(getVarCov(model)))^2 + model$sigma^2 )
  as.vector(icc)
}

#calculate total effect based on between- & within-family estimate and ICC
totaleffect <-  function(fmodel,imodel){
  coef(summary(fmodel))[2,1] * ICCest(imodel) + coef(summary(fmodel))[3,1]* (1- ICCest(imodel))
}


# Create between family estimates of PGS: average across family member
######################################################################

library(dplyr)
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

# Mixed effects model
######################

# intercept model
m0 <- lme(EA_sc~1, random=~1|FID, method="ML", na.action=na.omit,data=finalsib)
ICCest(m0) #get ICC #0.3258229

# include within and between family effect
m1 <- lme(EA_sc~GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog, random=~1|FID, method="ML", na.action=na.omit,data=finalsib)
totaleffect(m1, m0) #get total effect #0.2745393 
# check if this makes sense

summary(m1)

# Covariates in regression

final <- lme(EA_sc~GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog + sex + array + Age + sex*Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FID, method="ML", na.action=na.omit,data=finalsib)
summary(final)


# Extract estimates 
######################
names(summary(final))
summary(final)$tTable

direct_NonCog <- summary(final)$tTable[4,1]# direct effect is beta within 
direct_Cog <- summary(final)$tTable[5,1]
total_NonCog <- summary(final)$tTable[2,1] # total is beta between 
total_Cog <- summary(final)$tTable[3,1]
indirect_NonCog <- total_NonCog - direct_NonCog #0.1572753
indirect_Cog <- total_Cog - direct_Cog #0.1578887
ratio_NonCog <- indirect_NonCog/direct_NonCog #1.42775
ratio_Cog <- indirect_Cog/direct_Cog #1.278479

  
# Save data
#####################
write.table(finalsib, file="Data_scores_siblings_UKB_20200310.csv", row.names=F, quote=F) 


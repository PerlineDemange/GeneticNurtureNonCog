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

# save results
resultsib <- summary(final)$tTable
write.table(resultsib, file="Results_lme_siblings_UKB_20200401.csv", row.names=T, quote=F) 

#################
# Bootstrap 
#################

library(boot)
nboot <- 10000
bootcoef<-function(data,index){
  datx<-finalsib[index,]
  mod<-lme(EA_sc~GPS_B_NonCog + GPS_B_Cog + GPS_W_NonCog + GPS_W_Cog + sex + array + Age + sex*Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=~1|FID, method="ML", na.action=na.omit, data=datx,control=lmeControl(opt = "optim"))
  fixef(mod) #get fixed effects
}

# carry out bootstrap
boot.out<-boot(finalsib,bootcoef,nboot, parallel = "multicore", ncpus=20) 

#plot to check bootstrapping
png("UKB.sib.bootstrap.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out)
dev.off()

# get CI 
#boot.ci(boot.out)
#Error in bca.ci(boot.out, conf, index[1L], L = L, t = t.o, t0 = t0.o,  :
#                  estimated adjustment 'a' is NA
                
boot.ci(boot.out, type = c("norm", "basic"))

#Intervals :
#  Level      Normal              Basic
#95%   ( 1.025,  1.231 )   ( 1.014,  1.228 )


# save results for each bootstrap in data frame
bootoutput <- as.data.frame(boot.out$t)
colnames(bootoutput) <- rownames(as.data.frame(boot.out$t0))
write.table(bootoutput, "Data_scores_siblings_UKB_bootstrapped_20203031.csv", row.names=F, quote=F)
bmain <- bootoutput[,2:5] #Same variables of interest

# Create indirect and ratio variables

bmain$direct_NonCog <- bmain$GPS_W_NonCog
bmain$direct_Cog <- bmain$GPS_W_Cog
bmain$total_NonCog <- bmain$GPS_B_NonCog
bmain$total_Cog <- bmain$GPS_B_Cog
bmain$indirect_NonCog <- bmain$total_NonCog - bmain$direct_NonCog
bmain$indirect_Cog <- bmain$total_Cog - bmain$direct_Cog
bmain$ratio_NonCog <- bmain$indirect_NonCog / bmain$direct_NonCog
bmain$ratio_Cog <- bmain$indirect_Cog / bmain$direct_Cog
bmain <- bmain[,5:13]
write.table(bmain, "Data_scores_siblings_UKB_bootstrapped_ratio_20203031.csv", row.names=F, quote=F)

# perform t.tests

t.test(bmain$indirect_NonCog, bmain$direct_NonCog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_NonCog and bmain$direct_NonCog
# t = 332.04, df = 19776, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.04705787 0.04761675
# sample estimates:
#   mean of x mean of y
# 0.1573154 0.1099781


t.test(bmain$indirect_Cog, bmain$direct_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_Cog and bmain$direct_Cog
# t = 243.34, df = 19802, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.03494714 0.03551472
# sample estimates:
#   mean of x mean of y
# 0.1583185 0.1230876

t.test(bmain$indirect_NonCog, bmain$indirect_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$indirect_NonCog and bmain$indirect_Cog
# t = -6.6481, df = 19995, p-value = 3.045e-11
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.0012988103 -0.0007073346
# sample estimates:
#   mean of x mean of y
# 0.1573154 0.1583185


t.test(bmain$ratio_NonCog, bmain$ratio_Cog)
# Welch Two Sample t-test
# 
# data:  bmain$ratio_NonCog and bmain$ratio_Cog
# t = 50.843, df = 19493, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.1427075 0.1541520
# sample estimates:
#   mean of x mean of y
# 1.449193  1.300764


# Get CI of all 
meanall <- apply(bmain, 2, mean)
sdall <- apply(bmain, 2, sd)
n <- nboot
error <- qnorm(0.975)*sdall/sqrt(n)
leftCI <- meanall-error
rightCI <- meanall+error

tot <- rbind(meanall, sdall, error, leftCI, rightCI)
# direct_NonCog   direct_Cog total_NonCog    total_Cog indirect_NonCog
# meanall  0.1099780988 0.1230875527 2.672935e-01 2.814060e-01    0.1573154089
# sdall    0.0095313302 0.0097147855 4.564024e-03 4.524683e-03    0.0106021220
# error    0.0001868106 0.0001904063 8.945323e-05 8.868215e-05    0.0002077978
# leftCI   0.1097912881 0.1228971464 2.672041e-01 2.813174e-01    0.1571076111
# rightCI  0.1101649094 0.1232779590 2.673830e-01 2.814947e-01    0.1575232066
# indirect_Cog ratio_NonCog   ratio_Cog
# meanall 0.1583184813  1.449193289 1.300763521
# sdall   0.0107351602  0.222424405 0.189092945
# error   0.0002104053  0.004359438 0.003706154
# leftCI  0.1581080760  1.444833850 1.297057368
# rightCI 0.1585288866  1.453552727 1.304469675


write.table(tot, "summary_mean_CI_siblings_UKB_20200331.csv", row.names=T, quote=F)


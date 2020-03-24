# Rosa Cheesman
# PRs analysis of sib comparison in UKBiobank
# 24.02.20

#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////
####################
# make sibling id
####################
#take dataset with just sibs
#this greps loads of different files each containing occurrences of each ID from Sib_id_uniques:
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


#########################
# make file with FID
########################

R
library(data.table)
library(psych)

b <- fread("rels/sibs_famid3",h=F, data.table=F, verbose=T) #41493
colnames(b)[1:2]<-c("ID1","FID")
##siblings
c <- fread("../EA/siblings.csv", data.table=F, verbose=T, header=T) #39500
pheno <- c
pheno <- pheno[,c(1,3)]
colnames(pheno)<-c("ID1","EA")
dat2<-merge(pheno, b ,by = "ID1") #39500
poly <- fread('/project/ukbaumc/UKBGWAS/polygenic_scores/UKB.AMC.EA3_excl_23andMe_UK.HM3.EUR.SBLUP.10k.csv', data.table=F, header=T) 
colnames(poly) <- c('ID1', 'scoreEA')
finalsib <- merge(dat2, poly, by='ID1') #39500



#scale variables
finalsib[,c("EA_sc","scoreEA_sc")]<-apply(finalsib[,c("EA","scoreEA")],
                          2,
                          scale)



# create meanEA per fam 
library(dplyr)
mean<-group_by(finalsib,FID) %>% summarize(m=mean(scoreEA_sc))
#19389 rows
finalsibmeanprs<-merge(finalsib,mean,by="FID")




#now startsib difference model
library(nlme)
library(boot)

#################
#functions written 
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

#################
#Mixed model--- polygenic p factor predicting life experiences age 21 composite
#################

#create between-family variable
#DZ_scaled$GPS_B <- (DZ_scaled$PC_gps + DZ_scaled$PC_gps2)/2
#when there is more than one sibling in the family, the avg family GPS should reflect this..
#see qbove

finalsibmeanprs$GPS_B  <- finalsibmeanprs$m  
#create within-family variable
finalsibmeanprs$GPS_W <- finalsibmeanprs$scoreEA_sc  - finalsibmeanprs$GPS_B  

#Mixed effects model
#intercept model
#GPS not incl. yet...
m0 <- lme(EA_sc~1, random=~1|FID, method="ML", na.action=na.omit,data=finalsibmeanprs)
ICCest(m0) #get ICC

#include within and between family effect
m1 <- lme(EA_sc~GPS_B+GPS_W, random=~1|FID, method="ML", na.action=na.omit,data=finalsibmeanprs)
totaleffect(m1) #get total effect

summary(m1)



# add covariates and include in regression 












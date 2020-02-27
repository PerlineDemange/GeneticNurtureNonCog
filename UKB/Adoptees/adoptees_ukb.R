
#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////
#get phenotype for adoptees, nonadoptees, siblings
module add bioinformatics/R/3.3.3
R
library(data.table)
library(psych)

a <- fread("Allvars.txt",h=T, data.table=F)
#describe(a)

a[a==-818] <- NA
a[a==-121] <- NA
# a[a==-1] <- NA #dnno
# a[a==-3] <- NA #pref no ans
# a[a==-7] <- NA #none of above

x<-a$Qualifications.0.0
a$Edu<-ifelse(is.na(x), NA,
              ifelse(!is.na(x) & x == -7, 7,
                     ifelse(!is.na(x) & x == 4 | x == 3, 10,
                            ifelse(!is.na(x) & x == 2, 13,
                                   ifelse(!is.na(x) & x == 6, 15,
                                          ifelse(!is.na(x) & x ==5, 19,
                                                 ifelse(!is.na(x) & x ==1, 20, NA)))))))
#      7     10     13     15     19     20 
# 85291 132110  55331  25810  32734 161198 

#scale the phnenotype
dat3$PRS_std=scale(dat3$PRS)

#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////
#merge phenotype, PRS, and covariates, and family role for adoptees, nonadoptees, siblings


#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////
#run PRS analyses:
#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////
#subset to adoptees
adopteds<-subset(dat3,Adopted==1)
#scale the polygenic score within the subsample
dat3$PRS_std=scale(dat3$PRS)
library(boot)
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
} 
# bootstrapping with 1000 replications 
results <- boot(data=adopteds, statistic=rsq, 
                R=1000, formula=Edu~PRS_std)
# view results
results 
#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////
#subset to nonadoptees
nonadopteds<-subset(dat3,Adopted==0)
#scale the polygenic score within the subsample
dat3$PRS_std=scale(dat3$PRS)
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
} 
# bootstrapping with 1000 replications 
results <- boot(data=adopteds, statistic=rsq, 
                R=1000, formula=Edu~PRS_std)
results 

#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////
#subset to siblings
#create wide dataset grouping siblings by family






#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////
##saskia sib comparison script:
###nb this doesn't include covariates, they should be regressed out first...

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

#see max no. of siblings in a single family e.g. 5
# then sum over 5 columns, ignoring NAs.



#create within-family variable
DZ_scaled$GPS_W <- DZ_scaled$PC_gps  - DZ_scaled$GPS_B  

#Mixed effects model
#intercept model
#GPS not incl. yet...
m0 <- lme(exp1~1, random=~1|id_fam, method="ML", na.action=na.omit,data=DZ_scaled)
ICCest(m0) #get ICC

#include within and between family effect
m1 <- lme(exp1~GPS_B+GPS_W, random=~1|id_fam, method="ML", na.action=na.omit,data=DZ_scaled)
totaleffect(m1) #get total effect

summary(m1)
#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////#/////////////////////////////////////////////////////////////////////////




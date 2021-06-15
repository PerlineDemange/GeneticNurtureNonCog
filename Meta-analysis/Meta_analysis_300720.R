# july 2020
# meta-analysis of indirect and direct effects of cog/noncog pgs on educational outcomes
# Rosa Cheesman & Perline Demange 
# based on script from Margherita Malanchini

rm(list = ls())
require(metafor)
require(ggplot2)
require(ggsci)

# read in data (took summarised results supplementary file and split for cog non cog)
#need to redo if results are not up to date.
library(readxl)
NONCOG_PGS <- read_excel("NONCOG_summarised_results_20210519.xlsx")
COG_PGS <- read_excel("COG_summarised_results_20210519.xlsx")

##############################################################################################

# 1-- GET META ANALYTIC RESULTS OVERALL, DIRECT VS INDIRECT
### NONCOG PGS:
### direct
noncog_all_direct <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "direct"),]
meta_reg_noncog_all_direct <- rma.mv(yi = original,
                                     V = diag(se^2),
                                     random= ~ 1| Method,
                                     data=noncog_all_direct,
                                     slab = paste(Method, Cohort, Phenotype, sep=", "))
summary(meta_reg_noncog_all_direct)
forest(meta_reg_noncog_all_direct)


#Compared with:
# meta_reg_noncog_all_direct <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Cohort,data=noncog_all_direct,slab = paste(Method, Cohort, Phenotype))
# meta_reg_noncog_all_direct <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Cohort/Method,data=noncog_all_direct,slab = paste(Method, Cohort, Phenotype))
# meta_reg_noncog_all_direct <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Cohort/Method/Phenotype,data=noncog_all_direct,slab = paste(Method, Cohort, Phenotype))
# meta_reg_noncog_all_direct <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Phenotype,data=noncog_all_direct,slab = paste(Method, Cohort, Phenotype))
# #without any random effect
# meta_reg_noncog_all_direct <- rma.mv(yi = original,V = diag(se^2),data=noncog_all_direct,slab = paste(Method, Cohort, Phenotype))

#model with random intercept for methods has lowest AICc

### indirect
noncog_all_indirect <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "indirect"),]
meta_reg_noncog_all_indirect <- rma.mv(yi = original,V = diag(se^2),random= ~ 1|Method,data=noncog_all_indirect,slab = paste(Method, Cohort, Phenotype))
summary(meta_reg_noncog_all_indirect)
forest(meta_reg_noncog_all_indirect)

#### COG PGS:
### direct
cog_all_direct <- COG_PGS[c(COG_PGS$Type ==  "direct"),]
meta_reg_cog_all_direct <- rma.mv(yi = original, V = diag(se^2),random= ~ 1|Method,data=cog_all_direct,slab = paste(Method, Cohort, Phenotype, sep=", "))
summary(meta_reg_cog_all_direct)
forest(meta_reg_cog_all_direct)

### indirect
cog_all_indirect <- COG_PGS[c(COG_PGS$Type ==  "indirect"),]
meta_reg_cog_all_indirect <- rma.mv(yi = original, V = diag(se^2),random= ~ 1|Method,data=cog_all_indirect,slab = paste(Method, Cohort, Phenotype, sep=", "))
summary(meta_reg_cog_all_indirect)
forest(meta_reg_cog_all_indirect)

##############################################################################################

# 2-- GET META ANALYTIC RESULTS DIRECT VS INDIRECT per cohort
### NONCOG PGS:

###UKB
### direct
noncog_direct_UKB <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "direct"& NONCOG_PGS$Cohort ==  "UKB"),]
#doesnt make sense to have cohort random effect when that does not vary...
#so i put method as random effect, as that is the thing that always varies within cohort
meta_reg_noncog_direct_UKB <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=noncog_direct_UKB,slab = Method)
summary(meta_reg_noncog_direct_UKB)
forest(meta_reg_noncog_direct_UKB)
### indirect
noncog_indirect_UKB <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "indirect"& NONCOG_PGS$Cohort ==  "UKB"),]
meta_reg_noncog_indirect_UKB <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=noncog_indirect_UKB,slab = Method)
summary(meta_reg_noncog_indirect_UKB)
forest(meta_reg_noncog_indirect_UKB)

###TEDS
## no random effect as method does not vary
### direct
noncog_direct_TEDS <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "direct"& NONCOG_PGS$Cohort ==  "TEDS"),]
meta_reg_noncog_direct_TEDS <- rma.mv(yi = original,V = diag(se^2),data=noncog_direct_TEDS,slab = Method)
summary(meta_reg_noncog_direct_TEDS)
forest(meta_reg_noncog_direct_TEDS)
### indirect
noncog_indirect_TEDS <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "indirect"& NONCOG_PGS$Cohort ==  "TEDS"),]
meta_reg_noncog_indirect_TEDS <- rma.mv(yi = original,V = diag(se^2),data=noncog_indirect_TEDS,slab = Method)
summary(meta_reg_noncog_indirect_TEDS)
forest(meta_reg_noncog_indirect_TEDS)

###NTR
### direct
noncog_direct_NTR <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "direct"& NONCOG_PGS$Cohort ==  "NTR"),]
meta_reg_noncog_direct_NTR <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=noncog_direct_NTR,slab = Method)
summary(meta_reg_noncog_direct_NTR)
forest(meta_reg_noncog_direct_NTR)
### indirect
noncog_indirect_NTR <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "indirect"& NONCOG_PGS$Cohort ==  "NTR"),]
meta_reg_noncog_indirect_NTR <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=noncog_indirect_NTR,slab = Method)
summary(meta_reg_noncog_indirect_NTR)
forest(meta_reg_noncog_indirect_NTR)

#### COG PGS:

###UKB
### direct
cog_direct_UKB <- COG_PGS[c(COG_PGS$Type ==  "direct"& COG_PGS$Cohort ==  "UKB"),]
meta_reg_cog_direct_UKB <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=cog_direct_UKB,slab = Method)
summary(meta_reg_cog_direct_UKB)
### indirect
cog_indirect_UKB <- COG_PGS[c(COG_PGS$Type ==  "indirect"& COG_PGS$Cohort ==  "UKB"),]
meta_reg_cog_indirect_UKB <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=cog_indirect_UKB,slab = Method)
summary(meta_reg_cog_indirect_UKB)

###TEDS
### direct
cog_direct_TEDS <- COG_PGS[c(COG_PGS$Type ==  "direct"& COG_PGS$Cohort ==  "TEDS"),]
meta_reg_cog_direct_TEDS <- rma.mv(yi = original,V = diag(se^2),data=cog_direct_TEDS,slab = Method)
summary(meta_reg_cog_direct_TEDS)
### indirect
cog_indirect_TEDS <- COG_PGS[c(COG_PGS$Type ==  "indirect"& COG_PGS$Cohort ==  "TEDS"),]
meta_reg_cog_indirect_TEDS <- rma.mv(yi = original,V = diag(se^2),data=cog_indirect_TEDS,slab = Method)
summary(meta_reg_cog_indirect_TEDS)

###NTR
### direct
cog_direct_NTR <- COG_PGS[c(COG_PGS$Type ==  "direct"& COG_PGS$Cohort ==  "NTR"),]
meta_reg_cog_direct_NTR <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=cog_direct_NTR,slab = Method)
summary(meta_reg_cog_direct_NTR)
### indirect
cog_indirect_NTR <- COG_PGS[c(COG_PGS$Type ==  "indirect"& COG_PGS$Cohort ==  "NTR"),]
meta_reg_cog_indirect_NTR <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=cog_indirect_NTR,slab = Method)
summary(meta_reg_cog_indirect_NTR)

##############################################################################################

# 3-- GET META ANALYTIC RESULTS for DIRECT VS INDIRECT PER METHOD
### NONCOG PGS
### SIBLING DESIGN
### direct
#wihtout random for cohort lower aicc
noncog_direct_sib <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "direct"& NONCOG_PGS$Method ==  "Siblings"),]
meta_reg_noncog_direct_sib <- rma.mv(yi = original,V = diag(se^2),data=noncog_direct_sib,slab = Cohort)
summary(meta_reg_noncog_direct_sib)
forest(meta_reg_noncog_direct_sib)

### indirect
noncog_indirect_sib <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "indirect"& NONCOG_PGS$Method ==  "Siblings"),]
meta_reg_noncog_indirect_sib <- rma.mv(yi = original,V = diag(se^2),data=noncog_indirect_sib,slab = Cohort)
summary(meta_reg_noncog_indirect_sib)


### ADOPTION DESIGN --# only one cohort/outcome so no need

###TRIO DESIGN
### direct
noncog_direct_trio <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "direct"& NONCOG_PGS$Method ==  "Trios"),]
meta_reg_noncog_direct_trio <- rma.mv(yi = original,V = diag(se^2),data=noncog_direct_trio,slab = Phenotype)
summary(meta_reg_noncog_direct_trio)
### indirect
noncog_indirect_trio <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "indirect"& NONCOG_PGS$Method ==  "Trios"),]
meta_reg_noncog_indirect_trio <- rma.mv(yi = original,V = diag(se^2),data=noncog_indirect_trio,slab = Phenotype)
summary(meta_reg_noncog_indirect_trio)

#### COG PGS:
### SIBLING DESIGN
### direct
cog_direct_sib <- COG_PGS[c(COG_PGS$Type ==  "direct"& COG_PGS$Method ==  "Siblings"),]
meta_reg_cog_direct_sib <- rma.mv(yi = original,V = diag(se^2),data=cog_direct_sib,slab = Cohort)
summary(meta_reg_cog_direct_sib)
### indirect
cog_indirect_sib <- COG_PGS[c(COG_PGS$Type ==  "indirect"& COG_PGS$Method ==  "Siblings"),]
meta_reg_cog_indirect_sib <- rma.mv(yi = original,V = diag(se^2),data=cog_indirect_sib,slab = Cohort)
summary(meta_reg_cog_indirect_sib)

### ADOPTION DESIGN --# only one cohort/outcome so no need

###TRIO DESIGN
### direct
cog_direct_trio <- COG_PGS[c(COG_PGS$Type ==  "direct"& COG_PGS$Method ==  "Trios"),]
meta_reg_cog_direct_trio <- rma.mv(yi = original,V = diag(se^2),data=cog_direct_trio,slab = Phenotype)
summary(meta_reg_cog_direct_trio)
### indirect
cog_indirect_trio <- COG_PGS[c(COG_PGS$Type ==  "indirect"& COG_PGS$Method ==  "Trios"),]
meta_reg_cog_indirect_trio <- rma.mv(yi = original,V = diag(se^2),data=cog_indirect_trio,slab = Phenotype)
summary(meta_reg_cog_indirect_trio)

##############################################################################################

# 4-- GET META ANALYTIC RESULTS for DIRECT VS INDIRECT PER PHENOTYPIC OUTCOME
### NONCOG PGS
### EDUCATIONAL ATTAINMENT
### direct
noncog_direct_EA <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "direct"& NONCOG_PGS$Phenotype ==  "EA"),]
meta_reg_noncog_direct_EA <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method, data=noncog_direct_EA,slab = Cohort)
summary(meta_reg_noncog_direct_EA)
### indirect
noncog_indirect_EA <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "indirect"& NONCOG_PGS$Phenotype ==  "EA"),]
meta_reg_noncog_indirect_EA <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=noncog_indirect_EA,slab = Cohort)
summary(meta_reg_noncog_indirect_EA)

##CHILDHOOD ACHIEVEMENT
### direct
noncog_direct_child <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "direct"& NONCOG_PGS$Phenotype ==  "CITO"| NONCOG_PGS$Type ==  "direct"&NONCOG_PGS$Phenotype =="12yo"|NONCOG_PGS$Type ==  "direct"&NONCOG_PGS$Phenotype =="GCSE"),]
meta_reg_noncog_direct_child <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=noncog_direct_child,slab = Cohort)
summary(meta_reg_noncog_direct_child)
### indirect
noncog_indirect_child <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "indirect"& NONCOG_PGS$Phenotype ==  "CITO"| NONCOG_PGS$Type ==  "indirect"&NONCOG_PGS$Phenotype =="12yo"|NONCOG_PGS$Type ==  "indirect"&NONCOG_PGS$Phenotype =="GCSE"),]
meta_reg_noncog_indirect_child <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=noncog_indirect_child,slab = Cohort)
summary(meta_reg_noncog_indirect_child)

### COG PGS
### EDUCATIONAL ATTAINMENT
### direct
cog_direct_EA <- COG_PGS[c(COG_PGS$Type ==  "direct"& COG_PGS$Phenotype ==  "EA"),]
meta_reg_cog_direct_EA <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=cog_direct_EA,slab = Cohort)
summary(meta_reg_cog_direct_EA)
### indirect
cog_indirect_EA <- COG_PGS[c(COG_PGS$Type ==  "indirect"& COG_PGS$Phenotype ==  "EA"),]
meta_reg_cog_indirect_EA <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=cog_indirect_EA,slab = Cohort)
summary(meta_reg_cog_indirect_EA)

##CHILDHOOD ACHIEVEMENT
### direct
cog_direct_child <- COG_PGS[c(COG_PGS$Type ==  "direct"& COG_PGS$Phenotype ==  "CITO"| COG_PGS$Type ==  "direct"&COG_PGS$Phenotype =="12yo"|COG_PGS$Type ==  "direct"&COG_PGS$Phenotype =="GCSE"),]
meta_reg_cog_direct_child <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=cog_direct_child,slab = Cohort)
summary(meta_reg_cog_direct_child)
### indirect
cog_indirect_child <- COG_PGS[c(COG_PGS$Type ==  "indirect"& COG_PGS$Phenotype ==  "CITO"| COG_PGS$Type ==  "indirect"&COG_PGS$Phenotype =="12yo"|COG_PGS$Type ==  "indirect"&COG_PGS$Phenotype =="GCSE"),]
meta_reg_cog_indirect_child <- rma.mv(yi = original,V = diag(se^2),random= ~ 1| Method,data=cog_indirect_child,slab = Cohort)
summary(meta_reg_cog_indirect_child)


##############################################################################################

# 6-- GET META ANALYTIC RESULTS for DIRECT VS INDIRECT PER METHOD ONLY EA 

# Get only EA 
head(NONCOG_PGS)
NONCOG_PGS_EA <- NONCOG_PGS[NONCOG_PGS$Phenotype == "EA", ]
### NONCOG PGS
### SIBLING DESIGN
### direct
#wihtout random for cohort lower aicc
noncog_direct_sib_EA <- NONCOG_PGS_EA[c(NONCOG_PGS_EA$Type ==  "direct" & NONCOG_PGS_EA$Method ==  "Siblings"),]
meta_reg_noncog_direct_sib_EA <- rma.mv(yi = original,V = diag(se^2),data=noncog_direct_sib_EA,slab = Cohort)
summary(meta_reg_noncog_direct_sib_EA)
forest(meta_reg_noncog_direct_sib_EA)

### indirect
noncog_indirect_sib_EA <- NONCOG_PGS_EA[c(NONCOG_PGS_EA$Type ==  "indirect"& NONCOG_PGS_EA$Method ==  "Siblings"),]
meta_reg_noncog_indirect_sib_EA <- rma.mv(yi = original,V = diag(se^2),data=noncog_indirect_sib_EA,slab = Cohort)
summary(meta_reg_noncog_indirect_sib_EA)


### ADOPTION DESIGN --# only one cohort/outcome so no need

###TRIO DESIGN --# only one cohort/outcome

#### COG PGS:
COG_PGS_EA <- COG_PGS[COG_PGS$Phenotype == "EA", ]
### SIBLING DESIGN
### direct
cog_direct_sib_EA <- COG_PGS_EA[c(COG_PGS_EA$Type ==  "direct"& COG_PGS_EA$Method ==  "Siblings"),]
meta_reg_cog_direct_sib_EA <- rma.mv(yi = original,V = diag(se^2),data=cog_direct_sib_EA,slab = Cohort)
summary(meta_reg_cog_direct_sib_EA)
### indirect
cog_indirect_sib_EA <- COG_PGS_EA[c(COG_PGS_EA$Type ==  "indirect"& COG_PGS_EA$Method ==  "Siblings"),]
meta_reg_cog_indirect_sib_EA <- rma.mv(yi = original,V = diag(se^2),data=cog_indirect_sib_EA,slab = Cohort)
summary(meta_reg_cog_indirect_sib_EA)

### ADOPTION DESIGN --# only one cohort/outcome so no need

###TRIO DESIGN -- #only one cohort/outocme

##############################################################################################

# 7 -- GET META ANALYTIC RESULTS OVERALL, total effect
### NONCOG PGS:
noncog_all_tot <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "total"),]
meta_reg_noncog_all_tot <- rma.mv(yi = original,
                                     V = diag(se^2),
                                     random= ~ 1| Method,
                                     data=noncog_all_tot,
                                     slab = paste(Method, Cohort, Phenotype, sep=", "))
summary(meta_reg_noncog_all_tot)
forest(meta_reg_noncog_all_tot)

cog_all_tot <- COG_PGS[c(COG_PGS$Type ==  "total"),]
meta_reg_cog_all_tot <- rma.mv(yi = original,
                                  V = diag(se^2),
                                  random= ~ 1| Method,
                                  data=cog_all_tot,
                                  slab = paste(Method, Cohort, Phenotype, sep=", "))
summary(meta_reg_cog_all_tot)
forest(meta_reg_cog_all_tot)


###############################################################
# 8 -- GET META ANALYTIC RESULTS OVERALL, ratio

noncog_all_ratio <- NONCOG_PGS[c(NONCOG_PGS$Type ==  "ratio_tot"),]
meta_reg_noncog_all_ratio <- rma.mv(yi = original,
                                  V = diag(se^2),
                                  random= ~ 1| Method,
                                  data=noncog_all_ratio,
                                  slab = paste(Method, Cohort, Phenotype, sep=", "))
summary(meta_reg_noncog_all_ratio)
forest(meta_reg_noncog_all_ratio)

cog_all_ratio <- COG_PGS[c(COG_PGS$Type ==  "ratio_tot"),]
meta_reg_cog_all_ratio <- rma.mv(yi = original,
                               V = diag(se^2),
                               random= ~ 1| Method,
                               data=cog_all_ratio,
                               slab = paste(Method, Cohort, Phenotype, sep=", "))
summary(meta_reg_cog_all_ratio)
forest(meta_reg_cog_all_ratio)


##############################################################################################
# COLLATE RESULTS INTO A FILE FOR PLOTTING

# the following code is targeted towards stack plot.

#results all
meta_all<- as.data.frame(cbind(c(replicate(2, "NonCog"),replicate(2, "Cog")),
                                  rep(c('Direct', 'Indirect'), 2)))

colnames(meta_all)[1:2] <- c("PGS","Effect")

meta_all2 <- rbind(coef(summary(meta_reg_noncog_all_direct)), 
      coef(summary(meta_reg_noncog_all_indirect)),
      coef(summary(meta_reg_cog_all_direct)),
      coef(summary(meta_reg_cog_all_indirect)))

meta_all <- cbind(meta_all, meta_all2)
write.csv(meta_all, "meta_all_20210519.csv", row.names=F)

#2. results per sample
meta_sample<- as.data.frame(cbind(c(replicate(6, "NonCog"),replicate(6, "Cog")),
                    rep(c('Direct', 'Indirect'), 6),
                    rep(c('UKB', 'TEDS','NTR'), 2, each=2)))

colnames(meta_sample)[1:3] <- c("PGS","Effect","Sample")

#OR get everything with coef(summary(model))
meta <- rbind(coef(summary(meta_reg_noncog_direct_UKB)),
              coef(summary(meta_reg_noncog_indirect_UKB)), 
              coef(summary(meta_reg_noncog_direct_TEDS)), 
              coef(summary(meta_reg_noncog_indirect_TEDS)), 
              coef(summary(meta_reg_noncog_direct_NTR)), 
              coef(summary(meta_reg_noncog_indirect_NTR)),
              coef(summary(meta_reg_cog_direct_UKB)),
              coef(summary(meta_reg_cog_indirect_UKB)), 
              coef(summary(meta_reg_cog_direct_TEDS)), 
              coef(summary(meta_reg_cog_indirect_TEDS)), 
              coef(summary(meta_reg_cog_direct_NTR)), 
              coef(summary(meta_reg_cog_indirect_NTR))) 

meta_sample <- cbind(meta_sample, meta)
meta_sample
write.csv(meta_sample, "meta_sample_20210519.csv", row.names=F)


#then do the same stratifying by method:
meta_method<- as.data.frame(cbind(c(replicate(4, "NonCog"),replicate(4, "Cog")),
                                  rep(c('Direct', 'Indirect'), 4),
                                  rep(c('sib', 'trio'), 2, each=2)))

colnames(meta_method)[1:3] <- c("PGS","Effect","Sample")
meta_method2 <- rbind(coef(summary(meta_reg_noncog_direct_sib)),
                     coef(summary(meta_reg_noncog_indirect_sib)), 
                     coef(summary(meta_reg_noncog_direct_trio)), 
                     coef(summary(meta_reg_noncog_indirect_trio)), 
                     coef(summary(meta_reg_cog_direct_sib)),
                     coef(summary(meta_reg_cog_indirect_sib)), 
                     coef(summary(meta_reg_cog_direct_trio)), 
                     coef(summary(meta_reg_cog_indirect_trio))) 


meta_method <- cbind(meta_method, meta_method2)
meta_method
write.csv(meta_method, "meta_method_20210519.csv", row.names=F)


#4 per outcome (attainment adult vs child achievement)
meta_outcome<- as.data.frame(cbind(c(replicate(4, "NonCog"),replicate(4, "Cog")),
                                  rep(c('Direct', 'Indirect'), 4),
                                  rep(c('EA', 'child'), 2, each=2)))

colnames(meta_outcome)[1:3] <- c("PGS","Effect","Sample")
meta_outcome2 <- rbind(coef(summary(meta_reg_noncog_direct_EA)),
                      coef(summary(meta_reg_noncog_indirect_EA)), 
                      coef(summary(meta_reg_noncog_direct_child)), 
                      coef(summary(meta_reg_noncog_indirect_child)), 
                      coef(summary(meta_reg_cog_direct_EA)),
                      coef(summary(meta_reg_cog_indirect_EA)), 
                      coef(summary(meta_reg_cog_direct_child)), 
                      coef(summary(meta_reg_cog_indirect_child))) 

meta_outcome <- cbind(meta_outcome, meta_outcome2)
meta_outcome
write.csv(meta_outcome, "meta_outcome_20210519.csv", row.names=F)


# 5 per method only EA 
meta_outcome<- as.data.frame(cbind(c(replicate(2, "NonCog"),replicate(2, "Cog")),
                                   rep(c('Direct', 'Indirect'), 2)))
colnames(meta_outcome)[1:2] <- c("PGS","Effect")
meta_methodEA2 <- rbind(coef(summary(meta_reg_noncog_direct_sib_EA)),
coef(summary(meta_reg_noncog_indirect_sib_EA)),
coef(summary(meta_reg_cog_direct_sib_EA)),
coef(summary(meta_reg_cog_indirect_sib_EA))) 
meta_methodEA <- cbind(meta_outcome, meta_methodEA2)
meta_methodEA
write.csv(meta_methodEA, "meta_methodEA_20210519.csv", row.names=F)

# 6 total effect
meta_total_effect <- rbind(coef(summary(meta_reg_cog_all_tot )), coef(summary(meta_reg_noncog_all_tot)))
meta_total_effect$PGS <- rbind("Cog", "NonCog")
meta_total_effect$Effect <- "Total"
meta_total_effect
write.csv(meta_total_effect, "meta_total_effect_20210519.csv", row.names=F)

# 7 ratio
meta_ratio <- rbind(coef(summary(meta_reg_cog_all_ratio )), coef(summary(meta_reg_noncog_all_ratio)))
meta_ratio$PGS <- rbind("Cog", "NonCog")
meta_ratio$Effect <- "Ratio"
meta_ratio
write.csv(meta_ratio, "meta_ratio_20210519.csv", row.names=F)

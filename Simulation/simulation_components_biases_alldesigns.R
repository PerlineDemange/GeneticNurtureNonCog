# Code for the simulation 
# Please find the function to run the simulation from line 5
# We set parameters and run the function 100 times from line 1121 
# We plot results from line 1158


# Function #################

run_sim <- function(sample_size, number_snp, var_g, var_mother, var_father, var_sibling, var_prenatal, r.assortment, env_conf){
  # Function to run create simulated data with different components of indirect genetic effects and biases
  # and estimate the parental indirect genetic effect using 4 different designs.
  # Additionally we compare 3 different versions of the sibling design, which shows us that using the population effect is the appropriate version
  # Arguments of the function are: 
  # sample_size: sample size of target population 
  # number_snp: number of snps per genotype used in the simulation
  # var_g: variance explained by the individual phenotype
  # var_mother: variance explained by the mother genotype (after birth)
  # var_father: variance explained by the father genotype (after birth)
  # var_sibling: variance explained by the sibling genotype 
  # var_prenatal: variance explained by the mother genotype (before birth)
  # r.assortment: correlation between parents' phenotypes 
  # env_conf: environmental confounding, mean of the second population used for simulating population stratification (mean of the first population is 0)
  
  # 1. WITHOUT ASSORTMENT AND POPSTRAT ------
  ## 1.1 CREATE true and GWAS PGS --------
  ### 1.1.1 Simulate true and GWAS SNP effects #### 
  
  number_snp <- number_snp
  
  # Make 'true' effect sizes for the effect of the SNPs on the trait
  true_eff <- rnorm(number_snp)
  
  # Make GWAS effect sizes for the effect of the SNPs on the trait
  shrink_eff <- rnorm(number_snp)
  GWAS_eff <- sqrt(.2)*true_eff + sqrt(.8)*shrink_eff #GWAS effects are based on the true SNP effects with some error (shrinkage) 
  m.error <- cor(GWAS_eff,true_eff)
  
  # make true MAF's for number_snp SNPs
  maf <- runif(number_snp,.1,.5)
  
  ### 1.1.2 Create genotypes ###############
  # make bi-allelic SNP calls from number_snp SNPs in sample_size mothers, sample_size fathers
  # and their kids, a non-transmitted PRS and an adopted kid
  n <- sample_size
  
  # make the first "SNP" for n people:
  mothers <-  rbinom(n,size = 2,prob = maf[1])
  fathers <-  rbinom(n,size = 2,prob = maf[1])
  # make child & non-transmitted allele:
  ft <- rbinom(n,size = 1,prob = fathers/2)
  mt <- rbinom(n,size = 1,prob = mothers/2)
  children <- ft + mt
  ntc <- (fathers - ft) + (mothers - mt) 
  
  # make some sibs (for sib design)
  ft <- rbinom(n,size = 1,prob = fathers/2)
  mt <- rbinom(n,size = 1,prob = mothers/2)
  sibs <- ft + mt
  
  # Make some adoptees (i.e drawn fully independently from the rest):
  biological_mothers <- rbinom(n,size = 2,prob = maf[1])
  biological_fathers <- rbinom(n,size = 2,prob = maf[1])
  # make child & non-transmitted allele:
  bft <- rbinom(n,size = 1,prob = biological_fathers/2)
  bmt <- rbinom(n,size = 1,prob = biological_mothers/2)
  adoptees <- ft + mt
  

  # make SNP 2 to number_snp:
  for(i in 2:number_snp){ 
    
    mother <- rbinom(n,size = 2,prob = maf[i])
    father <- rbinom(n,size = 2,prob = maf[i])
    ft <- rbinom(n,size = 1,prob = father/2) # draw 1 allele from the father
    mt <- rbinom(n,size = 1,prob = mother/2) # draw 1 allele from the mother
    child  <- ft + mt # make the child SNP
    nt <- (father - ft) + (mother- mt) # make the non-transmitted genotype 
    
    mothers <- cbind(mothers,mother) # add SNP to file
    fathers <- cbind(fathers,father)
    children <- cbind(children,child)
    ntc <- cbind(ntc,nt)
    
    # For the sib design, make a SNP for a second child:
    ft <- rbinom(n,size = 1,prob = father/2)
    mt <- rbinom(n,size = 1,prob = mother/2)
    sib  <- ft + mt
    
    sibs <- cbind(sibs,sib) # add the data
    
    # For the adoptee design make a SNP for the adoptee:
    biological_mother <- rbinom(n,size = 2,prob = maf[i])
    biological_father <- rbinom(n,size = 2,prob = maf[i])
    
    bmt <- rbinom(n,size = 1,prob = biological_mother/2)
    bft <- rbinom(n,size = 1,prob = biological_father/2)
    
    adoptee <- bft + bmt
    adoptees <- cbind(adoptees,adoptee)
    biological_mothers <- cbind(biological_mothers,biological_mother) # add SNP to file
    biological_fathers <- cbind(biological_fathers,biological_father)
  }
  
  ### 1.1.3 Create 'true' genetic scores ##### 
  # multiply the beta and the SNPs and sum to a perfect PRS:
  mother_g <- true_eff %*% t(mothers)
  father_g<- true_eff %*% t(fathers)
  biological_mother_g <- true_eff %*% t(biological_mothers)
  biological_father_g<- true_eff %*% t(biological_fathers)
  child_g<- true_eff %*% t(children)
  ntc_g <-  true_eff %*% t(ntc)
  sib_g <- true_eff %*% t(sibs)
  adoptee_g <-  true_eff %*% t(adoptees)


  # scale genetic scores:
  mother_g <- scale(t(mother_g))
  father_g <- scale(t(father_g))
  biological_mother_g <- scale(t(biological_mother_g))
  biological_father_g <- scale(t(biological_father_g))
  child_g <- scale(t(child_g))
  ntc_g <- scale(t(ntc_g))
  sib_g <- scale(t(sib_g))
  adoptee_g <- scale(t(adoptee_g))

  ### 1.1.4 Create GWAS polygenetic scores - as if measured ##### 
  # multiply the beta and the SNPs and sum to a PRS:
  mother_prs_noam <- GWAS_eff %*% t(mothers)
  father_prs_noam <- GWAS_eff %*% t(fathers)
  child_prs_noam <- GWAS_eff %*% t(children)
  ntc_prs_noam <- GWAS_eff %*% t(ntc)
  sib_prs_noam <- GWAS_eff %*% t(sibs)
  adoptee_prs_noam <- GWAS_eff %*% t(adoptees)
  
  # scale genetic scores:
  mother_prs_noam <- scale(t(mother_prs_noam))
  father_prs_noam <- scale(t(father_prs_noam))
  child_prs_noam <- scale(t(child_prs_noam))
  ntc_prs_noam <- scale(t(ntc_prs_noam))
  sib_prs_noam <- scale(t(sib_prs_noam))
  adoptee_prs_noam <- scale(t(adoptee_prs_noam))
  
  ## 1.2. SIMULATE TRAITS -------
  # Traits are simulated using the true genetic scores 
  
  ### 1.2.1 No indirect effects: ----------------------------------------------------------------------------------------
  
  var_e <- 1 - var_g
  
  mother_trait_noind <- sqrt(var_g)*mother_g + sqrt(var_e)*rnorm(n)
  father_trait_noind <- sqrt(var_g)*father_g + sqrt(var_e)*rnorm(n)
  child_trait_noind <- sqrt(var_g)*child_g + sqrt(var_e)*rnorm(n) 
  sib_trait_noind <- sqrt(var_g)*sib_g + sqrt(var_e)*rnorm(n)
  adoptee_trait_noind <- sqrt(var_g)*adoptee_g + sqrt(var_e)*rnorm(n) 
  
  
  ### 1.2.2 With indirect effects: ------------------------------------------------------------------------------------------------------------------
  
  var_e <- 1 - var_g
  var_e_parental <- 1 - (var_father + var_mother)
  
  mother_trait_ind <- sqrt(var_g)*mother_g  + sqrt(var_e)*rnorm(n)
  father_trait_ind  <- sqrt(var_g)*father_g  + sqrt(var_e)*rnorm(n)
  child_trait_ind  <- sqrt(var_g)*child_g   + sqrt(var_mother)*mother_trait_ind + sqrt(var_father)*father_trait_ind + sqrt(var_e_parental)*rnorm(n)  
  sib_trait_ind     <- sqrt(var_g)*sib_g     + sqrt(var_mother)*mother_trait_ind + sqrt(var_father)*father_trait_ind + sqrt(var_e_parental)*rnorm(n)
  adoptee_trait_ind <- sqrt(var_g)*adoptee_g + sqrt(var_mother)*mother_trait_ind + sqrt(var_father)*father_trait_ind + sqrt(var_e_parental)*rnorm(n)
  

  
  ### 1.2.3. With sibling indirect effects: -------------------------------------------------------------------------------------------------------
  
  var_e <- 1 - var_g
  
  mother_trait_sib  <- sqrt(var_g)*mother_g  + sqrt(var_e)*rnorm(n)
  father_trait_sib  <- sqrt(var_g)*father_g  + sqrt(var_e)*rnorm(n)  
  sib_trait_sib     <- sqrt(var_g)*sib_g     + sqrt(var_e)*rnorm(n)
  child_trait_sib   <- sqrt(var_g)*child_g + sqrt(var_e)*rnorm(n)  
  adoptee_trait_sib <- sqrt(var_g)*adoptee_g + sqrt(var_e)*rnorm(n)
  
  
  # let's siblings influence each other:
  temp_scores <- cbind(child_trait_sib,sib_trait_sib,adoptee_trait_sib)
  
  # siblign indirect parameter
  s_e <- var_sibling
  sib_eff <- matrix(c(0,s_e,s_e,
                      s_e,0,s_e,
                      s_e,s_e,0),3,3,byrow=TRUE)
  identity <- matrix(c(1,0,0,
                       0,1,0,
                       0,0,1),3,3,byrow=TRUE)
  sib_eff <- solve(identity - sib_eff)

  post_sib_effect <- t(sib_eff %*% t(temp_scores))
  child_trait_sib <- post_sib_effect[,1]
  sib_trait_sib <- post_sib_effect[,2]
  adoptee_trait_sib <- post_sib_effect[,3]
  

  
  ### 1.2.4. With parental and sibling indirect effects:-----------------------------------------------------------------------------------------
  var_e <- 1 - var_g
  var_e_parental <- 1 - var_father - var_mother 
  
  mother_trait_ind_sib  <- sqrt(var_g)*mother_g  + sqrt(var_e)*rnorm(n)
  father_trait_ind_sib  <- sqrt(var_g)*father_g  + sqrt(var_e)*rnorm(n)  
  sib_trait_ind_sib    <- sqrt(var_g)*sib_g     + sqrt(var_mother)*mother_trait_ind_sib + sqrt(var_father)*father_trait_ind_sib + sqrt(var_e_parental)*rnorm(n)
  child_trait_ind_sib  <- sqrt(var_g)*child_g   + sqrt(var_mother)*mother_trait_ind_sib + sqrt(var_father)*father_trait_ind_sib + sqrt(var_e_parental)*rnorm(n)  
  adoptee_trait_ind_sib <- sqrt(var_g)*adoptee_g + sqrt(var_mother)*mother_trait_ind_sib + sqrt(var_father)*father_trait_ind_sib + sqrt(var_e_parental)*rnorm(n)
  
  # let's siblings influence each other:
  temp_scores <- cbind(child_trait_ind_sib,sib_trait_ind_sib,adoptee_trait_ind_sib)
  
  # siblign indirect parameter
  s_e <- var_sibling
  sib_eff <- matrix(c(0,s_e,s_e,
                      s_e,0,s_e,
                      s_e,s_e,0),3,3,byrow=TRUE)
  identity <- matrix(c(1,0,0,
                       0,1,0,
                       0,0,1),3,3,byrow=TRUE)
  sib_eff <- solve(identity - sib_eff)
  post_sib_effect <- t(sib_eff %*% t(temp_scores))
  child_trait_ind_sib <- post_sib_effect[,1]
  sib_trait_ind_sib <- post_sib_effect[,2]
  adoptee_trait_ind_sib <- post_sib_effect[,3]
  

  ### 1.2.5. With prenatal and postnatal indirect effects: ------------------------------------------------------------------------------------------------------------------
  
  var_e <- 1 - var_g
  var_e_parental_prenatal <- 1 - (var_father + var_mother + var_prenatal)
  
  mother_trait_ind_prenatal <- sqrt(var_g)*mother_g  + sqrt(var_e)*rnorm(n)
  father_trait_ind_prenatal  <- sqrt(var_g)*father_g  + sqrt(var_e)*rnorm(n)
  biological_mother_trait_ind_prenatal <- sqrt(var_g)*biological_mother_g  + sqrt(var_e)*rnorm(n)
  
  child_trait_ind_prenatal  <- sqrt(var_g)*child_g  + sqrt(var_prenatal)*mother_trait_ind_prenatal + sqrt(var_mother)*mother_trait_ind_prenatal + sqrt(var_father)*father_trait_ind_prenatal + sqrt(var_e_parental)*rnorm(n)  
  sib_trait_ind_prenatal     <- sqrt(var_g)*sib_g    + sqrt(var_prenatal)*mother_trait_ind_prenatal + sqrt(var_mother)*mother_trait_ind_prenatal + sqrt(var_father)*father_trait_ind_prenatal + sqrt(var_e_parental)*rnorm(n)
  adoptee_trait_ind_prenatal <- sqrt(var_g)*adoptee_g + sqrt(var_prenatal)*biological_mother_trait_ind_prenatal + sqrt(var_mother)*mother_trait_ind_prenatal + sqrt(var_father)*father_trait_ind_prenatal + sqrt(var_e_parental)*rnorm(n)
  
  
  
  # 2. With assortative mating: ------------------------------------------------------------------------------------------------------------------
  ## 2.1 CREATE true and GWAS PGS ------
  
  ### 2.1.2 Create genotype ######
  # make bi-allelic SNP calls from 100 SNPs in mothers, fathers
  # and their kids, a non-transmitted PRS and an adopted kids (which have their own biological parents simulated # make the first "SNP" for n people:
  mothers <- rbinom(n,size = 2,prob = maf[1])
  fathers <- rbinom(n,size = 2,prob = maf[1])
  # make biological parents for the adoptees:
  biological_mothers <- rbinom(n,size = 2,prob = maf[1])
  biological_fathers <- rbinom(n,size = 2,prob = maf[1])
  # make SNP 2 to 100:
  for(i in 2:number_snp){
    mother <- rbinom(n,size = 2,prob = maf[i])
    father <- rbinom(n,size = 2,prob = maf[i])
    biological_mother <- rbinom(n,size = 2,prob = maf[i])
    biological_father <- rbinom(n,size = 2,prob = maf[i])
    mothers <- cbind(mothers,mother) # add SNP to file
    fathers <- cbind(fathers,father)
    biological_mothers <- cbind(biological_mothers,biological_mother) # add SNP to file
    biological_fathers <- cbind(biological_fathers,biological_father)
  }
  # multiply the beta and the SNPs and sum to a perfect PRS:
  mother_g <- true_eff %*% t(mothers)
  father_g<- true_eff %*% t(fathers)
  # multiply the beta and the SNPs and sum to a perfect PRS:
  biological_mother_g <- true_eff %*% t(biological_mothers)
  biological_father_g<- true_eff %*% t(biological_fathers)
  # scale genetic scores:
  mother_g <- scale(t(mother_g))
  father_g <- scale(t(father_g))
  # scale genetic scores:
  biological_mother_g <- scale(t(biological_mother_g))
  biological_father_g <- scale(t(biological_father_g))
  
  ### 2.1.2 Make the parent phenotypes #####
  var_e <- 1 - var_g 
  var_e_parental <- 1 - (var_father + var_mother)
  
  mother_trait_assortative <- sqrt(var_g)*mother_g + sqrt(var_e)*rnorm(n)
  father_trait_assortative <- sqrt(var_g)*father_g + sqrt(var_e)*rnorm(n)
  biological_mother_trait_assortative <- sqrt(var_g)*biological_mother_g + sqrt(var_e)*rnorm(n)
  biological_father_trait_assortative <- sqrt(var_g)*biological_father_g + sqrt(var_e)*rnorm(n)
  
  ### 2.1.3. ASSORTMENT #####
  # phenotype correlation between mates
  r = r.assortment
  vp <- (var(mother_trait_assortative) + var(father_trait_assortative)) / 2
  bvp <- (var(biological_mother_trait_assortative) + var(biological_father_trait_assortative)) / 2
  
  # create mate selection noise, proportional to r (ie. decreases with r)
  #var_mate = (-2 * r * vp + sqrt(2 * r^2 * vp^2 - 2 * r^2 * vp^2 + 4 * vp^2)) / (2 * r)
  var_mate = vp * (-r + 1)  / r
  #bvar_mate = (-2 * r * bvp + sqrt(2 * r^2 * bvp^2 - 2 * r^2 * bvp^2 + 4 * bvp^2)) / (2 * r)
  bvar_mate = bvp * (-r + 1)  / r
  # mother and father are independently rank ordered by phenotype (with noise), then matched according to rank_pm = order(mother_t4 + (as.vector(sqrt(var_mate)) * rnorm(n)))
  rank_pm = order(mother_trait_assortative + (as.vector(sqrt(var_mate)) * rnorm(n)))
  rank_pf = order(father_trait_assortative + (as.vector(sqrt(var_mate)) * rnorm(n)))
  # mother and father are independently rank ordered by phenotype (with noise), then matched according to brank_pm = order(biological_mother_t4 + (as.vector(sqrt(var_mate)) * rnorm(n)))
  brank_pm = order(biological_mother_trait_assortative + (as.vector(sqrt(var_mate)) * rnorm(n)))
  brank_pf = order(biological_father_trait_assortative + (as.vector(sqrt(var_mate)) * rnorm(n)))
  
  # rank parents
  mothers = mothers[rank_pm,]
  fathers = fathers[rank_pf,]
  # rank bio parents
  biological_mothers = biological_mothers[brank_pm,]
  biological_fathers = biological_fathers[brank_pf,]
  # jointly unrank parents
  samp1 <- sample(1:n,n,FALSE)
  mothers = mothers[samp1 ,]
  fathers = fathers[samp1 ,]
  samp2 <- sample(1:n,n,FALSE)
  biological_mothers = biological_mothers[samp2 ,]
  biological_fathers = biological_fathers[samp2 ,]
  
  #Now that parents are assorted we can create the offspring
  #make child & non-transmitted allele:
  ft <- rbinom(n,size = 1,prob = fathers[,1]/2)
  mt <- rbinom(n,size = 1,prob = mothers[,1]/2)
  children <- ft + mt
  ntc <- (fathers[,1] - ft) + (mothers[,1] - mt)
  # make some sibs (for sib design)
  ft <- rbinom(n,size = 1,prob = fathers/2)
  mt <- rbinom(n,size = 1,prob = mothers/2)
  sibs <- ft + mt
  # For the adoptee design make a SNP for the adoptee, from the biological parents:
  bft <- rbinom(n,size = 1,prob = biological_fathers[,1]/2) # draw 1 allele from the father
  bmt <- rbinom(n,size = 1,prob = biological_mothers[,1]/2) # draw 1 allele from the mother
  adoptees <- bft + bmt # make the child SNP
  # make SNP 2 to number_snp:
  for(i in 2:number_snp){
    ft <- rbinom(n,size = 1,prob = fathers[,i]/2) # draw 1 allele from the father
    mt <- rbinom(n,size = 1,prob = mothers[,i]/2) # draw 1 allele from the mother
    child <- ft + mt # make the child SNP
    nt <- (fathers[,i] - ft) + (mothers[,i]- mt) # make the non-transmitted genotype
    children <- cbind(children,child)
    ntc <- cbind(ntc,nt)
    # For the sib design, make a SNP for a second child:
    ft <- rbinom(n,size = 1,prob = fathers[,i]/2)
    mt <- rbinom(n,size = 1,prob = mothers[,i]/2)
    sib <- ft + mt
    sibs <- cbind(sibs,sib) # add the data
    # For the adoptee design make a SNP for the adoptee, from the biological parents:
    bft <- rbinom(n,size = 1,prob = biological_fathers[,i]/2) # draw 1 allele from the father
    bmt <- rbinom(n,size = 1,prob = biological_mothers[,i]/2) # draw 1 allele from the mother
    adoptee <- bft + bmt # make the child SNP
    adoptees <- cbind(adoptees,adoptee)
  }
    
  # Reorder the parental phenotypes before influencing the kids
  mother_trait_assortative <- mother_trait_assortative[rank_pm][samp1]
  father_trait_assortative <- father_trait_assortative[rank_pf][samp1]
  biological_mother_trait_assortative <- biological_mother_trait_assortative[brank_pm][samp2]
  biological_father_trait_assortative <- biological_father_trait_assortative[brank_pf][samp2]
  
  ### 2.1.4 Create true genetic scores ####
  # multiply the beta and the SNPs and sum to a perfect PRS:
  child_g<- true_eff %*% t(children)
  ntc_g <- true_eff %*% t(ntc)
  sib_g <- true_eff %*% t(sibs)
  adoptee_g <- true_eff %*% t(adoptees)
  # scale genetic scores:
  child_g <- scale(t(child_g))
  ntc_g <- scale(t(ntc_g))
  sib_g <- scale(t(sib_g))
  adoptee_g <- scale(t(adoptee_g))
  
  ### 2.1.5 Create GWAS genetic scores ####
  # multiply the beta and the SNPs and sum to a PRS:
  mother_prs_am <- GWAS_eff %*% t(mothers)
  father_prs_am <- GWAS_eff %*% t(fathers)
  child_prs_am <- GWAS_eff %*% t(children)
  ntc_prs_am <- GWAS_eff %*% t(ntc)
  sib_prs_am <- GWAS_eff %*% t(sibs)
  adoptee_prs_am <- GWAS_eff %*% t(adoptees)
  
  # scale genetic scores:
  mother_prs_am <- scale(t(mother_prs_am))
  father_prs_am<- scale(t(father_prs_am))
  child_prs_am <- scale(t(child_prs_am))
  ntc_prs_am <- scale(t(ntc_prs_am))
  sib_prs_am <- scale(t(sib_prs_am))
  adoptee_prs_am <- scale(t(adoptee_prs_am))  
  
  ## 2.2. SIMULATE TRAITS -------
  ### 2.2.1 With assortative mating but no indirect effect -----------------------------------------------------------------------------------------------------------------  
  
  child_trait_assortative <- sqrt(var_g)*child_g + sqrt(var_e)*rnorm(n) 
  sib_trait_assortative <- sqrt(var_g)*sib_g + sqrt(var_e)*rnorm(n)
  adoptee_trait_assortative <- sqrt(var_g)*adoptee_g + sqrt(var_e)*rnorm(n) 
    
  
  ### 2.2.2. With assortative mating and indirect effect -----------------------------------------------------------------------------------------------------------------
  
  mother_trait_assortative_ind <- mother_trait_assortative
  father_trait_assortative_ind <- father_trait_assortative
  biological_mother_trait_assortative_ind <- biological_mother_trait_assortative
  biological_father_trait_assortative_ind <- biological_father_trait_assortative
  
  child_trait_assortative_ind <- sqrt(var_g)*child_g + sqrt(var_mother)*mother_trait_assortative_ind + sqrt(var_father)*father_trait_assortative_ind + sqrt(var_e_parental)*rnorm(n) 
  sib_trait_assortative_ind <- sqrt(var_g)*sib_g + sqrt(var_mother)*mother_trait_assortative_ind + sqrt(var_father)*father_trait_assortative_ind + sqrt(var_e_parental)*rnorm(n)
  adoptee_trait_assortative_ind <- sqrt(var_g)*adoptee_g + sqrt(var_mother)*mother_trait_assortative_ind + sqrt(var_father)*father_trait_assortative_ind + sqrt(var_e_parental)*rnorm(n) 
  
  # 3. Population stratification -------------------------------------------------------------------------------------------------------------------
  ## 3.1 CREATE PGS -----------------------
  
  ### 3.1.1  Make two subpopulations with different MAF ####
  n <- .5* sample_size #(we make 2 pops, each .5n)
  
  # make true MAF's for number_snp SNPs
  maf <- runif(number_snp,.1,.5)
  maf2 <- runif(number_snp,.1,.5)
  
  ### 3.1.2 Run a GWAS to get GWAS SNP effects (with popstrat)#####
  # The GWAS is run in a population containing two subpopulations,
  # with different maf and different phenotype mean
  SNPs1 <- rbinom(n,size = 2,prob = maf[1])
  SNPs2 <- rbinom(n,size = 2,prob = maf2[1])
  for(i in 2:number_snp){ 
    SNP1 <- rbinom(n,size = 2,prob = maf[i])
    SNP2 <- rbinom(n,size = 2,prob = maf2[i])
    
    SNPs1 <- cbind(SNPs1,SNP1)
    SNPs2 <- cbind(SNPs2,SNP2)
  }
  
  SNPs <- rbind(SNPs1,SNPs2)
  SNP_g <-  true_eff %*% t(SNPs)
  
  env_conf <- env_conf
  GWAS_t  <- sqrt(var_g)*SNP_g + sqrt(1-var_g)*as.numeric(scale(c(rnorm(n),rnorm(n,env_conf,1)))) 
  #this creates different noise for both populations, with one population with a mean of 0 and the second population with a mean of 1  
  GWAS_t <- as.vector(GWAS_t)
  GWAS_eff <- lm(GWAS_t ~ SNPs)  
  
  GWAS_eff1 <- as.numeric(GWAS_eff$coefficients)[2:(number_snp + 1)]
  
  ### Shrink the precision of the effect, like in real GWAS:
  GWAS_eff <- sqrt(.2)*GWAS_eff1 + sqrt(.8)*shrink_eff
  m.error2 <- cor(GWAS_eff,GWAS_eff1)
  

  ### 3.1.3 Create the genotypes of the target populations #### 
  # the target population contains the same two subpopulations as the GWAS 
  
  # POP 1
  # make the first "SNP" for n people:
  mothers <-  rbinom(n,size = 2,prob = maf[1])
  fathers <-  rbinom(n,size = 2,prob = maf[1])
  # make child & non-transmitted allele:
  ft <- rbinom(n,size = 1,prob = fathers/2)
  mt <- rbinom(n,size = 1,prob = mothers/2)
  children <- ft + mt
  ntc <- (fathers - ft) + (mothers- mt) 
  
  # make some sibs (for sib design)
  ft <- rbinom(n,size = 1,prob = fathers/2)
  mt <- rbinom(n,size = 1,prob = mothers/2)
  sibs <- ft + mt
  
  # Make some adoptees (i.e drawn fully independently from the rest):
  # make the first "SNP" for n people:
  bmothers <-  rbinom(n,size = 2,prob = maf[1])
  bfathers <-  rbinom(n,size = 2,prob = maf[1])
  # make child & non-transmitted allele:
  bft <- rbinom(n,size = 1,prob = bfathers/2)
  bmt <- rbinom(n,size = 1,prob = bmothers/2)
  adoptees <- ft + mt
  
  # make SNP 2 to number_snp:
  for(i in 2:number_snp){ 
    
    mother <- rbinom(n,size = 2,prob = maf[i])
    father <- rbinom(n,size = 2,prob = maf[i])
    ft <- rbinom(n,size = 1,prob = father/2) # draw 1 allele from the father
    mt <- rbinom(n,size = 1,prob = mother/2) # draw 1 allele from the mother
    child  <- ft + mt # make the child SNP
    nt <- (father - ft) + (mother- mt) # make the non-transmitted genotype 
    
    mothers <- cbind(mothers,mother) # add SNP to file
    fathers <- cbind(fathers,father)
    children <- cbind(children,child)
    ntc <- cbind(ntc,nt)
    
    # For the sib design, make a SNP for a second child:
    ft <- rbinom(n,size = 1,prob = father/2)
    mt <- rbinom(n,size = 1,prob = mother/2)
    sib  <- ft + mt
    
    sibs <- cbind(sibs,sib) # add the data
    
    # For the adoptee design make a SNP for the adoptee:
    
    bmother <-  rbinom(n,size = 2,prob = maf[i])
    bfather <-  rbinom(n,size = 2,prob = maf[i])
    # make child & non-transmitted allele:
    bft <- rbinom(n,size = 1,prob = bfather/2)
    bmt <- rbinom(n,size = 1,prob = bmother/2)
    adoptee <- bft + bmt
    adoptees <- cbind(adoptees,adoptee)
    bmothers <- cbind(bmothers,bmother) # add SNP to file
    bfathers <- cbind(bfathers,bfather)
  }
  
  # POP 2
  # make the first "SNP" for n people:
  mothers.2 <-  rbinom(n,size = 2,prob = maf2[1])
  fathers.2 <-  rbinom(n,size = 2,prob = maf2[1])
  # make child.2 & non-transmitted allele:
  ft <- rbinom(n,size = 1,prob = fathers.2/2)
  mt <- rbinom(n,size = 1,prob = mothers.2/2)
  children.2 <- ft + mt
  ntc.2 <- (fathers.2 - ft) + (mothers.2- mt) 
  
  # make some sibs (for sib design)
  ft <- rbinom(n,size = 1,prob = fathers.2/2)
  mt <- rbinom(n,size = 1,prob = mothers.2/2)
  sibs.2 <- ft + mt
  
  # Make some adoptees (i.e drawn fully independently from the rest):
  # make the first "SNP" for n people:
  bmothers.2 <-  rbinom(n,size = 2,prob = maf2[1])
  bfathers.2 <-  rbinom(n,size = 2,prob = maf2[1])
  # make child.2 & non-transmitted allele:
  bft <- rbinom(n,size = 1,prob = bfathers.2/2)
  bmt <- rbinom(n,size = 1,prob = bmothers.2/2)
  adoptees.2 <- ft + mt
  
  # make SNP 2 to number_snp:
  for(i in 2:number_snp){ 
    
    mother.2 <- rbinom(n,size = 2,prob = maf2[i])
    father.2 <- rbinom(n,size = 2,prob = maf2[i])
    ft <- rbinom(n,size = 1,prob = father.2/2) # draw 1 allele from the father.2
    mt <- rbinom(n,size = 1,prob = mother.2/2) # draw 1 allele from the mother.2
    child.2  <- ft + mt # make the child.2 SNP
    nt.2 <- (father.2 - ft) + (mother.2- mt) # make the non-transmitted genotype 
    
    mothers.2 <- cbind(mothers.2,mother.2) # add SNP to file
    fathers.2 <- cbind(fathers.2,father.2)
    children.2 <- cbind(children.2,child.2)
    ntc.2 <- cbind(ntc.2,nt.2)
    
    # For the sib design, make a SNP for a second child.2:
    ft <- rbinom(n,size = 1,prob = father.2/2)
    mt <- rbinom(n,size = 1,prob = mother.2/2)
    sib.2  <- ft + mt
    
    sibs.2 <- cbind(sibs.2,sib.2) # add the data
    
    # For the adoptee design make a SNP for the adoptee:
    
    bmother.2 <-  rbinom(n,size = 2,prob = maf2[i])
    bfather.2 <-  rbinom(n,size = 2,prob = maf2[i])
    # make child.2 & non-transmitted allele:
    bft <- rbinom(n,size = 1,prob = bfather.2/2)
    bmt <- rbinom(n,size = 1,prob = bmother.2/2)
    adoptee.2 <- bft + bmt
    adoptees.2 <- cbind(adoptees.2,adoptee.2)
    bmothers.2 <- cbind(bmothers.2,bmother.2) # add SNP to file
    bfathers.2 <- cbind(bfathers.2,bfather.2)
  }
  
  # Combine populations
  mothers <- rbind(mothers,mothers.2)
  fathers <- rbind(fathers,fathers.2)
  children <- rbind(children,children.2)
  sibs <- rbind(sibs,sibs.2)
  ntc <- rbind(ntc,ntc.2)
  bmothers <- rbind(bmothers,bmothers.2)
  bfathers <- rbind(bfathers,bfathers.2)
  adoptees <- rbind(adoptees,adoptees.2)
  
  # In our simulation we simulate the simple situation of adoptees being adopted in their own subpopulation
  
  #### 3.1.3.1 OPTIONAL BUT IMPORTANT, possibility of cross ancestry adoption #####
  # Several cases possible (among others): 
  # adoptees are adopted in the same subpopulation
  # adoptees are adopted in the alternative subpopulation 
  # part of subopulation is "moving" and has noise corresponding to the other population
  # 
# 
#   half10min <- sample_size/2 - sample_size*0.10
#   half10max <- sample_size/2 + sample_size*0.10
#   
#   adoptees <- rbind(adoptees[1:half10min,],adoptees[(half10max+1):sample_size,],adoptees[(half10min +1):half10max,])
#   bmothers <- rbind(bmothers[1:half10min,],bmothers[(half10max+1):sample_size,],bmothers[(half10min+1):half10max,])
#   bfathers <- rbind(bfathers[1:half10min,],bfathers[(half10max+1):sample_size,],bfathers[(half10min+1):half10max,])
#   
  
  
  ###3.1.4 Create true genetic score ####
  # multiply the beta and the SNPs and sum to a perfect PRS:
  mother_g <- true_eff %*% t(mothers)
  father_g<- true_eff %*% t(fathers)
  bmother_g <- true_eff %*% t(bmothers)
  bfather_g<- true_eff %*% t(bfathers)
  child_g<- true_eff %*% t(children)
  ntc_g <-  true_eff %*% t(ntc)
  sib_g <- true_eff %*% t(sibs)
  adoptee_g <-  true_eff %*% t(adoptees)
  
  # scale genetic scores:
  mother_g <- scale(t(mother_g))
  father_g <- scale(t(father_g))
  bmother_g <- scale(t(bmother_g))
  bfather_g <- scale(t(bfather_g))
  child_g <- scale(t(child_g))
  ntc_g <- scale(t(ntc_g))
  sib_g <- scale(t(sib_g))
  adoptee_g <- scale(t(adoptee_g))
  
  ### 3.1.5. Create the GWAS PGS #####
  # multiply the beta and the SNPs and sum to a PRS:
  mother_prs_popstrat <- GWAS_eff %*% t(mothers)
  father_prs_popstrat <- GWAS_eff %*% t(fathers)
  child_prs_popstrat <- GWAS_eff %*% t(children)
  ntc_prs_popstrat <-  GWAS_eff %*% t(ntc)
  sib_prs_popstrat <- GWAS_eff %*% t(sibs)
  adoptee_prs_popstrat <-  GWAS_eff %*% t(adoptees)
  
  # scale genetic scores:
  mother_prs_popstrat <- scale(t(mother_prs_popstrat))
  father_prs_popstrat <- scale(t(father_prs_popstrat))
  child_prs_popstrat <- scale(t(child_prs_popstrat))
  ntc_prs_popstrat <- scale(t(ntc_prs_popstrat))
  sib_prs_popstrat <- scale(t(sib_prs_popstrat))
  adoptee_prs_popstrat <- scale(t(adoptee_prs_popstrat))
  
  
  ## 3.2 SIMULATE TRAITS -------------------------------------------------------------------------------
  ### 3.2.1. With population stratification but no indirect effect -----------------------------------------------------------------------------------------------------------------  
  var_e <- 1 - var_g 
  
  mother_trait_popstrat <- sqrt(var_g)*mother_g + sqrt(var_e)*scale(c(rnorm(n),rnorm(n,env_conf,1)))
  father_trait_popstrat <- sqrt(var_g)*father_g + sqrt(var_e)*scale(c(rnorm(n),rnorm(n,env_conf,1)))
  bmother_trait_popstrat <- sqrt(var_g)*bmother_g + sqrt(var_e)*scale(c(rnorm(n),rnorm(n,env_conf,1)))
  
  child_trait_popstrat   <- sqrt(var_g)*child_g    + sqrt(var_e)*scale(c(rnorm(n),rnorm(n,env_conf,1)))
  sib_trait_popstrat     <- sqrt(var_g)*sib_g      + sqrt(var_e)*scale(c(rnorm(n),rnorm(n,env_conf,1)))
  adoptee_trait_popstrat <- sqrt(var_g)*adoptee_g  + sqrt(var_e)*scale(c(rnorm(n),rnorm(n,env_conf,1)))
  
  
  
  ### 3.2.2. With population stratification and parental indirect effect -----------------------------------------------------------------------------------------------------------------  
  var_e_parental <- 1 - (var_father + var_mother)
  
  child_trait_popstrat_ind   <- sqrt(var_g)*child_g   + sqrt(var_mother)*mother_trait_popstrat + sqrt(var_father)*father_trait_popstrat + sqrt(var_e_parental)*scale(c(rnorm(n),rnorm(n,env_conf,1)))
  sib_trait_popstrat_ind     <- sqrt(var_g)*sib_g     + sqrt(var_mother)*mother_trait_popstrat + sqrt(var_father)*father_trait_popstrat + sqrt(var_e_parental)*scale(c(rnorm(n),rnorm(n,env_conf,1)))
  adoptee_trait_popstrat_ind  <- sqrt(var_g)*adoptee_g + sqrt(var_mother)*mother_trait_popstrat + sqrt(var_father)*father_trait_popstrat + sqrt(var_e_parental)*scale(c(rnorm(n),rnorm(n,env_conf,1)))
  
  # 4. RUN ANALYSES ---------------------------------------------------------------------------------------------------------------------------------------------
  ## 4.1 perform transmitted non-transmitted PRS analysis --------------------------------------------------------------------------
  kong_noind <- summary(lm(child_trait_noind ~ child_prs_noam + ntc_prs_noam))
  kong_ind <- summary(lm(child_trait_ind ~ child_prs_noam + ntc_prs_noam))
  kong_sib <- summary(lm(child_trait_sib ~ child_prs_noam + ntc_prs_noam))
  kong_ind_sib <- summary(lm(child_trait_ind_sib ~ child_prs_noam + ntc_prs_noam))
  kong_ind_prenatal <- summary(lm(child_trait_ind_prenatal ~ child_prs_noam + ntc_prs_noam))
  kong_assortative <- summary(lm(child_trait_assortative ~ child_prs_am + ntc_prs_am))
  kong_assortative_ind <- summary(lm(child_trait_assortative_ind ~ child_prs_am + ntc_prs_am))
  kong_popstrat <- summary(lm(child_trait_popstrat ~ child_prs_popstrat + ntc_prs_popstrat))
  kong_popstrat_ind <- summary(lm(child_trait_popstrat_ind ~ child_prs_popstrat + ntc_prs_popstrat))
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  direct_kong_noind <- kong_noind$coef[2,1] - kong_noind$coef[3,1]
  direct_kong_ind <- kong_ind$coef[2,1] - kong_ind$coef[3,1]
  direct_kong_sib <- kong_sib$coef[2,1] - kong_sib$coef[3,1]
  direct_kong_ind_sib <- kong_ind_sib$coef[2,1] - kong_ind_sib$coef[3,1]
  direct_kong_ind_prenatal <- kong_ind_prenatal$coef[2,1] - kong_ind_prenatal$coef[3,1]
  direct_kong_assortative <- kong_assortative$coef[2,1] - kong_assortative$coef[3,1]
  direct_kong_assortative_ind <- kong_assortative_ind$coef[2,1] - kong_assortative_ind$coef[3,1]
  direct_kong_popstrat <- kong_popstrat$coef[2,1] - kong_popstrat$coef[3,1]
  direct_kong_popstrat_ind <- kong_popstrat_ind$coef[2,1] - kong_popstrat_ind$coef[3,1]
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  indirect_kong_noind <-kong_noind$coef[3,1]
  indirect_kong_ind <- kong_ind$coef[3,1]
  indirect_kong_sib <- kong_sib$coef[3,1]
  indirect_kong_ind_sib <- kong_ind_sib$coef[3,1]
  indirect_kong_ind_prenatal <- kong_ind_prenatal$coef[3,1]
  indirect_kong_assortative<- kong_assortative$coef[3,1]
  indirect_kong_assortative_ind<- kong_assortative_ind$coef[3,1]
  indirect_kong_popstrat<- kong_popstrat$coef[3,1]
  indirect_kong_popstrat_ind<- kong_popstrat_ind$coef[3,1]
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  ## 4.2. perform adoption PGS analysis: ----------------------------------------------------------------------------
  cheesman_noind_adopt <- summary(lm(adoptee_trait_noind ~ adoptee_prs_noam ))
  cheesman_noind_bio <- summary(lm(child_trait_noind ~ child_prs_noam ))
  
  cheesman_ind_adopt <- summary(lm(adoptee_trait_ind ~ adoptee_prs_noam ))
  cheesman_ind_bio <- summary(lm(child_trait_ind ~ child_prs_noam ))
  
  cheesman_sib_adopt <- summary(lm(adoptee_trait_sib ~ adoptee_prs_noam ))
  cheesman_sib_bio <- summary(lm(child_trait_sib ~ child_prs_noam ))
  
  cheesman_ind_sib_adopt <- summary(lm(adoptee_trait_ind_sib ~ adoptee_prs_noam ))
  cheesman_ind_sib_bio <- summary(lm(child_trait_ind_sib ~ child_prs_noam ))
  
  cheesman_ind_prenatal_adopt <- summary(lm(adoptee_trait_ind_prenatal ~ adoptee_prs_noam ))
  cheesman_ind_prenatal_bio <- summary(lm(child_trait_ind_prenatal ~ child_prs_noam ))
  
  cheesman_assortative_adopt <- summary(lm(adoptee_trait_assortative ~ adoptee_prs_am ))
  cheesman_assortative_bio <- summary(lm(child_trait_assortative ~ child_prs_am ))
  
  cheesman_assortative_ind_adopt <- summary(lm(adoptee_trait_assortative_ind ~ adoptee_prs_am ))
  cheesman_assortative_ind_bio <- summary(lm(child_trait_assortative_ind ~ child_prs_am ))
  
  cheesman_popstrat_adopt <- summary(lm(adoptee_trait_popstrat ~ adoptee_prs_popstrat ))
  cheesman_popstrat_bio <- summary(lm(child_trait_popstrat ~ child_prs_popstrat ))
  
  cheesman_popstrat_ind_adopt <- summary(lm(adoptee_trait_popstrat_ind ~ adoptee_prs_popstrat ))
  cheesman_popstrat_ind_bio <- summary(lm(child_trait_popstrat_ind ~ child_prs_popstrat ))
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  direct_cheesman_noind <- cheesman_noind_adopt$coef[2,1]
  direct_cheesman_ind <- cheesman_ind_adopt$coef[2,1]
  direct_cheesman_sib <- cheesman_sib_adopt$coef[2,1]
  direct_cheesman_ind_sib <- cheesman_ind_sib_adopt$coef[2,1]
  direct_cheesman_ind_prenatal<- cheesman_ind_prenatal_adopt$coef[2,1]
  direct_cheesman_assortative<- cheesman_assortative_adopt$coef[2,1]
  direct_cheesman_assortative_ind<- cheesman_assortative_ind_adopt$coef[2,1]
  direct_cheesman_popstrat<- cheesman_popstrat_adopt$coef[2,1]
  direct_cheesman_popstrat_ind<- cheesman_popstrat_ind_adopt$coef[2,1]
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  indirect_cheesman_noind <- cheesman_noind_bio$coef[2,1] - cheesman_noind_adopt$coef[2,1]
  indirect_cheesman_ind <- cheesman_ind_bio$coef[2,1] - cheesman_ind_adopt$coef[2,1]
  indirect_cheesman_sib <- cheesman_sib_bio$coef[2,1] - cheesman_sib_adopt$coef[2,1]
  indirect_cheesman_ind_sib <- cheesman_ind_sib_bio$coef[2,1] - cheesman_ind_sib_adopt$coef[2,1]
  indirect_cheesman_ind_prenatal <- cheesman_ind_prenatal_bio$coef[2,1] - cheesman_ind_prenatal_adopt$coef[2,1]
  indirect_cheesman_assortative <- cheesman_assortative_bio$coef[2,1] - cheesman_assortative_adopt$coef[2,1]
  indirect_cheesman_assortative_ind <- cheesman_assortative_ind_bio$coef[2,1] - cheesman_assortative_ind_adopt$coef[2,1]
  indirect_cheesman_popstrat <- cheesman_popstrat_bio$coef[2,1] - cheesman_popstrat_adopt$coef[2,1]
  indirect_cheesman_popstrat_ind <- cheesman_popstrat_ind_bio$coef[2,1] - cheesman_popstrat_ind_adopt$coef[2,1]
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  ## 4.3 perform siblings PGS analysis: -----------------------------------------------------------------------
  # We run 3 different siblings designs: comparing the within-sibling effect to the between-sibling effect, 
  # to the population effect, and to the total effect (as defined in Selzam et al.)
  
  
  between_noam <- .5*(child_prs_noam + sib_prs_noam)
  within_noam <- child_prs_noam - between_noam
  
  between_am <- .5*(child_prs_am + sib_prs_am)
  within_am <- child_prs_am - between_am
  
  between_popstrat <- .5*(child_prs_popstrat + sib_prs_popstrat)
  within_popstrat <- child_prs_popstrat - between_popstrat
  
  selzam_noind <- summary(lm(child_trait_noind ~ between_noam + within_noam))
  selzam_ind <- summary(lm(child_trait_ind ~ between_noam + within_noam))
  selzam_sib <- summary(lm(child_trait_sib ~ between_noam + within_noam))
  selzam_ind_sib <- summary(lm(child_trait_ind_sib ~ between_noam + within_noam))
  selzam_ind_prenatal <- summary(lm(child_trait_ind_prenatal ~ between_noam + within_noam))
  selzam_assortative <- summary(lm(child_trait_assortative ~ between_am + within_am))
  selzam_assortative_ind <- summary(lm(child_trait_assortative_ind ~ between_am + within_am))
  selzam_popstrat <- summary(lm(child_trait_popstrat ~ between_popstrat + within_popstrat))
  selzam_popstrat_ind <- summary(lm(child_trait_popstrat_ind ~ between_popstrat + within_popstrat))
  
  icc_noind <- cor(child_trait_noind, sib_trait_noind)
  icc_ind <- cor(child_trait_ind, sib_trait_ind)
  icc_sib <- cor(child_trait_sib, sib_trait_sib)
  icc_ind_sib <- cor(child_trait_ind_sib, sib_trait_ind_sib)
  icc_ind_prenatal <- cor(child_trait_ind_prenatal, sib_trait_ind_prenatal)
  icc_assortative <- cor(child_trait_assortative, sib_trait_assortative)
  icc_assortative_ind <- cor(child_trait_assortative_ind, sib_trait_assortative_ind)
  icc_popstrat <- cor(child_trait_popstrat, sib_trait_popstrat)
  icc_popstrat_ind <- cor(child_trait_popstrat_ind, sib_trait_popstrat_ind)
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  direct_selzam_noind <- selzam_noind$coef[3,1]
  direct_selzam_ind <- selzam_ind$coef[3,1]
  direct_selzam_sib <- selzam_sib$coef[3,1]
  direct_selzam_ind_sib <- selzam_ind_sib$coef[3,1]
  direct_selzam_ind_prenatal <- selzam_ind_prenatal$coef[3,1]
  direct_selzam_assortative <- selzam_assortative$coef[3,1]
  direct_selzam_assortative_ind <- selzam_assortative_ind$coef[3,1]
  direct_selzam_popstrat <- selzam_popstrat$coef[3,1]
  direct_selzam_popstrat_ind <- selzam_popstrat_ind$coef[3,1]
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  
  indirect_selzam_noind <- selzam_noind$coef[2,1]-selzam_noind$coef[3,1] 
  indirect_selzam_ind <- selzam_ind$coef[2,1]-selzam_ind$coef[3,1]
  indirect_selzam_sib <- selzam_sib$coef[2,1]-selzam_sib$coef[3,1]
  indirect_selzam_ind_sib <- selzam_ind_sib$coef[2,1]-selzam_ind_sib$coef[3,1]
  indirect_selzam_ind_prenatal <- selzam_ind_prenatal$coef[2,1]-selzam_ind_prenatal$coef[3,1]
  indirect_selzam_assortative <- selzam_assortative$coef[2,1]-selzam_assortative$coef[3,1]
  indirect_selzam_assortative_ind <- selzam_assortative_ind$coef[2,1]-selzam_assortative_ind$coef[3,1]
  indirect_selzam_popstrat <- selzam_popstrat$coef[2,1]-selzam_popstrat$coef[3,1]
  indirect_selzam_popstrat_ind <- selzam_popstrat_ind$coef[2,1]-selzam_popstrat_ind$coef[3,1]
  
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  
  total_noind <- (direct_selzam_noind*(1-icc_noind))+ selzam_noind$coef[2,1]*icc_noind
  total_ind <- (direct_selzam_ind*(1-icc_ind))+ selzam_ind$coef[2,1]*icc_ind
  total_sib <- (direct_selzam_sib*(1-icc_sib))+ selzam_sib$coef[2,1]*icc_sib
  total_ind_sib <- (direct_selzam_ind_sib*(1-icc_ind_sib))+ selzam_ind_sib$coef[2,1]*icc_ind_sib
  total_ind_prenatal <- (direct_selzam_ind_prenatal*(1-icc_ind_prenatal))+ selzam_ind_prenatal$coef[2,1]*icc_ind_prenatal
  total_assortative <- (direct_selzam_assortative*(1-icc_assortative))+ selzam_assortative$coef[2,1]*icc_assortative
  total_assortative_ind <- (direct_selzam_assortative_ind*(1-icc_assortative_ind))+ selzam_assortative_ind$coef[2,1]*icc_assortative_ind
  total_popstrat <- (direct_selzam_popstrat*(1-icc_popstrat))+ selzam_popstrat$coef[2,1]*icc_popstrat
  total_popstrat_ind <- (direct_selzam_popstrat_ind*(1-icc_popstrat_ind))+ selzam_popstrat_ind$coef[2,1]*icc_popstrat_ind
  
  indirect_selzam_noind_total <- total_noind - direct_selzam_noind
  indirect_selzam_ind_total <- total_ind - direct_selzam_ind
  indirect_selzam_sib_total <- total_sib - direct_selzam_sib
  indirect_selzam_ind_sib_total <- total_ind_sib - direct_selzam_ind_sib
  indirect_selzam_ind_prenatal_total <- total_ind_prenatal - direct_selzam_ind_sib
  indirect_selzam_assortative_total <- total_assortative - direct_selzam_assortative
  indirect_selzam_assortative_ind_total <- total_assortative_ind - direct_selzam_assortative_ind
  indirect_selzam_popstrat_total <- total_popstrat - direct_selzam_popstrat
  indirect_selzam_popstrat_ind_total <- total_popstrat_ind - direct_selzam_popstrat_ind
  
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  pop_noind <- summary(lm(child_trait_noind ~ child_prs_noam))
  pop_ind <- summary(lm(child_trait_ind~ child_prs_noam))
  pop_sib <- summary(lm(child_trait_sib ~ child_prs_noam))
  pop_ind_sib <- summary(lm(child_trait_ind_sib ~ child_prs_noam))
  pop_ind_prenatal <- summary(lm(child_trait_ind_prenatal ~ child_prs_noam))
  pop_assortative <- summary(lm(child_trait_assortative ~ child_prs_am))
  pop_assortative_ind <- summary(lm(child_trait_assortative_ind ~ child_prs_am))
  pop_popstrat <- summary(lm(child_trait_popstrat ~ child_prs_popstrat))
  pop_popstrat_ind <- summary(lm(child_trait_popstrat_ind ~ child_prs_popstrat))
  
  
  indirect_selzam_noind_pop <- pop_noind$coef[2,1]-selzam_noind$coef[3,1] 
  indirect_selzam_ind_pop <- pop_ind$coef[2,1]-selzam_ind$coef[3,1]
  indirect_selzam_sib_pop <- pop_sib$coef[2,1]-selzam_sib$coef[3,1]
  indirect_selzam_ind_sib_pop <- pop_ind_sib$coef[2,1]-selzam_ind_sib$coef[3,1]
  indirect_selzam_ind_prenatal_pop <- pop_ind_prenatal$coef[2,1]-selzam_ind_prenatal$coef[3,1]
  indirect_selzam_assortative_pop <- pop_assortative$coef[2,1]-selzam_assortative$coef[3,1]
  indirect_selzam_assortative_ind_pop <- pop_assortative_ind$coef[2,1]-selzam_assortative_ind$coef[3,1]
  indirect_selzam_popstrat_pop <- pop_popstrat$coef[2,1]-selzam_popstrat$coef[3,1]
  indirect_selzam_popstrat_ind_pop <- pop_popstrat_ind$coef[2,1]-selzam_popstrat_ind$coef[3,1]
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  ## 4.4 perform trios analysis --------------------------------------------------------------------------------------
  # Suggested by the reviewer, similar to what was suggested by Conley et al. 
  conley_noind <- summary(lm(child_trait_noind ~ child_prs_noam + father_prs_noam + mother_prs_noam))
  conley_ind <- summary(lm(child_trait_ind ~ child_prs_noam + father_prs_noam + mother_prs_noam))
  conley_sib <- summary(lm(child_trait_sib ~ child_prs_noam + father_prs_noam + mother_prs_noam))
  conley_ind_sib <- summary(lm(child_trait_ind_sib ~ child_prs_noam + father_prs_noam + mother_prs_noam))
  conley_ind_prenatal <- summary(lm(child_trait_ind_prenatal ~ child_prs_noam + father_prs_noam + mother_prs_noam))
  conley_assortative <- summary(lm(child_trait_assortative ~ child_prs_am + father_prs_am + mother_prs_am))
  conley_assortative_ind <- summary(lm(child_trait_assortative_ind ~ child_prs_am + father_prs_am + mother_prs_am))
  conley_popstrat <- summary(lm(child_trait_popstrat ~ child_prs_popstrat + father_prs_popstrat + mother_prs_popstrat))
  conley_popstrat_ind <- summary(lm(child_trait_popstrat_ind ~ child_prs_popstrat + father_prs_popstrat + mother_prs_popstrat))
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  direct_conley_noind <- conley_noind$coef[2,1]
  direct_conley_ind <- conley_ind$coef[2,1]
  direct_conley_sib <- conley_sib$coef[2,1]
  direct_conley_ind_sib <- conley_ind_sib$coef[2,1]
  direct_conley_ind_prenatal <- conley_ind_prenatal$coef[2,1]
  direct_conley_assortative <- conley_assortative$coef[2,1]
  direct_conley_assortative_ind <- conley_assortative_ind$coef[2,1]
  direct_conley_popstrat <- conley_popstrat$coef[2,1]
  direct_conley_popstrat_ind <- conley_popstrat_ind$coef[2,1]
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  indirect_conley_noind <- .5*(conley_noind$coef[3,1] + conley_noind$coef[4,1])
  indirect_conley_ind <- .5*(conley_ind$coef[3,1] + conley_ind$coef[4,1])
  indirect_conley_sib <- .5*(conley_sib$coef[3,1] + conley_sib$coef[4,1])
  indirect_conley_ind_sib <- .5*(conley_ind_sib$coef[3,1] + conley_ind_sib$coef[4,1])
  indirect_conley_ind_prenatal <- .5*(conley_ind_prenatal$coef[3,1] + conley_ind_prenatal$coef[4,1])
  indirect_conley_assortative <- .5*(conley_assortative$coef[3,1] + conley_assortative$coef[4,1])
  indirect_conley_assortative_ind <- .5*(conley_assortative_ind$coef[3,1] + conley_assortative_ind$coef[4,1])
  indirect_conley_popstrat <- .5*(conley_popstrat$coef[3,1] + conley_popstrat$coef[4,1])
  indirect_conley_popstrat_ind <- .5*(conley_popstrat_ind$coef[3,1] + conley_popstrat_ind$coef[4,1])
  
  ## 4.5 Combine, compare and save ---------------------------------------------------------------------------------------------------------------------------------------------
  direct_effects_noind <- c(m.error*sqrt(var_g),
                         direct_kong_noind, 
                         direct_cheesman_noind,
                         direct_selzam_noind,
                         direct_selzam_noind,
                         direct_selzam_noind,
                         direct_conley_noind
  )
  direct_effects_ind <- c(m.error*sqrt(var_g),
                         direct_kong_ind, 
                         direct_cheesman_ind,
                         direct_selzam_ind,
                         direct_selzam_ind,
                         direct_selzam_ind,
                         direct_conley_ind
  )
  direct_effects_sib <- c(m.error*sqrt(var_g),
                           direct_kong_sib, 
                           direct_cheesman_sib,
                           direct_selzam_sib,
                           direct_selzam_sib,
                           direct_selzam_sib,
                           direct_conley_sib
  )
  direct_effects_ind_sib <- c(m.error*sqrt(var_g),
                            direct_kong_ind_sib, 
                            direct_cheesman_ind_sib,
                            direct_selzam_ind_sib,
                            direct_selzam_ind_sib,
                            direct_selzam_ind_sib,
                            direct_conley_ind_sib
  )
  direct_effects_ind_prenatal <- c(m.error*sqrt(var_g),
                                direct_kong_ind_prenatal, 
                                direct_cheesman_ind_prenatal,
                                direct_selzam_ind_prenatal,
                                direct_selzam_ind_prenatal,
                                direct_selzam_ind_prenatal,
                                direct_conley_ind_prenatal
  )
  direct_effects_assortative <- c(m.error*sqrt(var_g),
                                      direct_kong_assortative, 
                                      direct_cheesman_assortative,
                                      direct_selzam_assortative,
                                      direct_selzam_assortative,
                                      direct_selzam_assortative,
                                      direct_conley_assortative
  ) 
  
  direct_effects_assortative_ind <- c(m.error*sqrt(var_g),
                          direct_kong_assortative_ind, 
                          direct_cheesman_assortative_ind,
                          direct_selzam_assortative_ind,
                          direct_selzam_assortative_ind,
                          direct_selzam_assortative_ind,
                          direct_conley_assortative_ind
  ) 
  
  direct_effects_popstrat <- c(m.error2*sqrt(var_g),
                                  direct_kong_popstrat, 
                                  direct_cheesman_popstrat,
                                  direct_selzam_popstrat,
                                  direct_selzam_popstrat,
                                  direct_selzam_popstrat,
                                  direct_conley_popstrat
  ) 
  
  direct_effects_popstrat_ind <- c(m.error2*sqrt(var_g),
                                      direct_kong_popstrat_ind, 
                                      direct_cheesman_popstrat_ind,
                                      direct_selzam_popstrat_ind,
                                      direct_selzam_popstrat_ind,
                                      direct_selzam_popstrat_ind,
                                      direct_conley_popstrat_ind
  ) 
  
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  var_par <- (var_mother + var_father)/2
  indirect_effects_noind <- c(m.error*sqrt(0),
                           indirect_kong_noind, 
                           indirect_cheesman_noind,
                           indirect_selzam_noind,
                           indirect_selzam_noind_total,
                           indirect_selzam_noind_pop,
                           indirect_conley_noind
  )
  
  indirect_effects_ind <- c(m.error*sqrt(var_par)*sqrt(var_g),
                           indirect_kong_ind, 
                           indirect_cheesman_ind,
                           indirect_selzam_ind,
                           indirect_selzam_ind_total,
                           indirect_selzam_ind_pop,
                           indirect_conley_ind
  )
  
  
  indirect_effects_sib <- c(m.error*sqrt(0),
                             indirect_kong_sib, 
                             indirect_cheesman_sib,
                             indirect_selzam_sib,
                             indirect_selzam_sib_total,
                             indirect_selzam_sib_pop,
                             indirect_conley_sib
  )
  
  indirect_effects_ind_sib <- c(m.error*sqrt(var_par)*sqrt(var_g),
                              indirect_kong_ind_sib, 
                              indirect_cheesman_ind_sib,
                              indirect_selzam_ind_sib,
                              indirect_selzam_ind_sib_total,
                              indirect_selzam_ind_sib_pop,
                              indirect_conley_ind_sib
  )
  
  indirect_effects_ind_prenatal <- c(m.error*sqrt(var_par)*sqrt(var_g),
                                  indirect_kong_ind_prenatal, 
                                  indirect_cheesman_ind_prenatal,
                                  indirect_selzam_ind_prenatal,
                                  indirect_selzam_ind_prenatal_total,
                                  indirect_selzam_ind_prenatal_pop,
                                  indirect_conley_ind_prenatal
  )
  
  indirect_effects_assortative <- c(m.error*sqrt(0),
                                        indirect_kong_assortative, 
                                        indirect_cheesman_assortative,
                                        indirect_selzam_assortative,
                                        indirect_selzam_assortative_total,
                                        indirect_selzam_assortative_pop,
                                        indirect_conley_assortative
  )
  
  
  indirect_effects_assortative_ind <- c(m.error*sqrt(var_par)*sqrt(var_g),
                            indirect_kong_assortative_ind, 
                            indirect_cheesman_assortative_ind,
                            indirect_selzam_assortative_ind,
                            indirect_selzam_assortative_ind_total,
                            indirect_selzam_assortative_ind_pop,
                            indirect_conley_assortative_ind
  )
  
  indirect_effects_popstrat <- c(m.error2*sqrt(0),
                                    indirect_kong_popstrat, 
                                    indirect_cheesman_popstrat,
                                    indirect_selzam_popstrat,
                                    indirect_selzam_popstrat_total,
                                    indirect_selzam_popstrat_pop,
                                    indirect_conley_popstrat
  )
  
  
  indirect_effects_popstrat_ind <- c(m.error2*sqrt(var_par)*sqrt(var_g),
                                        indirect_kong_popstrat_ind, 
                                        indirect_cheesman_popstrat_ind,
                                        indirect_selzam_popstrat_ind,
                                        indirect_selzam_popstrat_ind_total,
                                        indirect_selzam_popstrat_ind_pop,
                                        indirect_conley_popstrat_ind
  )
  
  
  ## ---------------------------------------------------------------------------------------------------------------------------------------------
  tab_noind <-  cbind(round(direct_effects_noind,3),round(indirect_effects_noind,3))
  
  rownames(tab_noind) <- c("Truth", "Non-transmitted", "Adoption",
                        "Sibling", "Sibling - total","Sibling - population", "Reviewer")
  colnames(tab_noind) <-  c("Direct effect", "Parental indirect effect")
  
  tab_ind <-  cbind(round(direct_effects_ind,3),round(indirect_effects_ind,3))
  
  rownames(tab_ind) <- c("Truth", "Non-transmitted", "Adoption",
                        "Sibling", "Sibling - total","Sibling - population", "Reviewer")
  colnames(tab_ind) <-  c("Direct effect", "Parental indirect effect")  
  
  tab_sib <-  cbind(round(direct_effects_sib,3),round(indirect_effects_sib,3))
  
  rownames(tab_sib) <- c("Truth", "Non-transmitted", "Adoption",
                          "Sibling", "Sibling - total","Sibling - population", "Reviewer")
  colnames(tab_sib) <-  c("Direct effect", "Parental indirect effect")
  
  tab_ind_sib <-  cbind(round(direct_effects_ind_sib,3),round(indirect_effects_ind_sib,3))
  
  rownames(tab_ind_sib) <- c("Truth", "Non-transmitted", "Adoption",
                           "Sibling", "Sibling - total","Sibling - population", "Reviewer")
  colnames(tab_ind_sib) <-  c("Direct effect", "Parental indirect effect")
  
  tab_ind_prenatal <-  cbind(round(direct_effects_ind_prenatal,3),round(indirect_effects_ind_prenatal,3))
  
  rownames(tab_ind_prenatal) <- c("Truth", "Non-transmitted", "Adoption",
                               "Sibling", "Sibling - total","Sibling - population", "Reviewer")
  colnames(tab_ind_prenatal) <-  c("Direct effect", "Parental indirect effect")
  
  tab_assortative <-  cbind(round(direct_effects_assortative,3),round(indirect_effects_assortative,3))
  
  rownames(tab_assortative) <- c("Truth", "Non-transmitted", "Adoption",
                                     "Sibling", "Sibling - total","Sibling - population", "Reviewer")
  colnames(tab_assortative) <-  c("Direct effect", "Parental indirect effect")  
  
  
  tab_assortative_ind <-  cbind(round(direct_effects_assortative_ind,3),round(indirect_effects_assortative_ind,3))
  
  rownames(tab_assortative_ind) <- c("Truth", "Non-transmitted", "Adoption",
                         "Sibling", "Sibling - total","Sibling - population", "Reviewer")
  colnames(tab_assortative_ind) <-  c("Direct effect", "Parental indirect effect")  
  
  tab_popstrat <-  cbind(round(direct_effects_popstrat,3),round(indirect_effects_popstrat,3))
  
  rownames(tab_popstrat) <- c("Truth", "Non-transmitted", "Adoption",
                                 "Sibling", "Sibling - total","Sibling - population", "Reviewer")
  colnames(tab_popstrat) <-  c("Direct effect", "Parental indirect effect")  
  
  
  tab_popstrat_ind <-  cbind(round(direct_effects_popstrat_ind,3),round(indirect_effects_popstrat_ind,3))
  
  rownames(tab_popstrat_ind) <- c("Truth", "Non-transmitted", "Adoption",
                                     "Sibling", "Sibling - total","Sibling - population", "Reviewer")
  colnames(tab_popstrat_ind) <-  c("Direct effect", "Parental indirect effect")  
  
  # Print results ---------------------------------------------------------------------------------------------------------------------------------------------
  
  myFunction=function(){
    tab_noind=tab_noind
    tab_ind=tab_ind
    tab_sib=tab_sib
    tab_ind_sib=tab_ind_sib
    tab_ind_prenatal=tab_ind_prenatal
    tab_assortative=tab_assortative
    tab_assortative_ind=tab_assortative_ind
    tab_popstrat=tab_popstrat
    tab_popstrat_ind=tab_popstrat_ind
    sapply(ls(),function(x)get(x),simplify=F,USE.NAMES=T)
  }
  
  myResults=myFunction()
  
  myResults
} 

# Run simulation ##########################################################################

library(tidyverse)
library(data.table)

len <- 9 #length(all_sim) number of table in the output of run_sim
all_sim_total_direct  <- vector(mode = "list", length = len)
all_sim_total_indirect  <- vector(mode = "list", length = len)
for (sim in 1:100){ 
  all_sim <- run_sim(20000, 100, 0.5, 0.2, 0.2, 0.15, 0.1, 0.5, 3) 
  #print(length(all_sim))
  for (table in 1:length(all_sim)){
    all_sim_total_direct[[table]] <- cbind(all_sim_total_direct[[table]], all_sim[[table]][,1] )
    all_sim_total_indirect[[table]] <- cbind(all_sim_total_indirect[[table]], all_sim[[table]][,2] )
  }
}

all_simulations <- NULL
#order of the table was re-order by alphabetic order so: tab_assortative, tab_assortative_ind, tab_ind, tab_ind_prenatal, tab_ind_sib, tab_noind, tab_popstrat, tab_popstrat_ind, tab_sib

name_simulation<- c("AM", "AM_ind", "ind", "prenatal", "ind_sib", "noind", "popstrat", "popstrat_ind", "sib")
for (simulation in 1:len){
  direct <- as.data.frame(t(all_sim_total_direct[[simulation]]))
  direct$effect <- "direct"
  indirect <- as.data.frame(t(all_sim_total_indirect[[simulation]]))
  indirect$effect <- "indirect"
  total <- rbind(direct, indirect)
  total$simulation <- name_simulation[[simulation]]
  all_simulations <- rbind(all_simulations, total)
}


save(all_simulations, file="all_simulations_100_20000_biggereffects_popstrat_210616.Rda")
load("all_simulations_100_20000_biggereffects_popstrat_210616.Rda")


# Figures #############################################################################
library(reshape2)
library(ggplot2)
library(tidyverse)
head(all_simulations)

names(all_simulations)[names(all_simulations) == "Reviewer"] <- "Trios"
total_results_long<- melt(data          = all_simulations,
                          id.vars       = c("effect", "simulation"),
                          measure.vars  = c("Truth", "Non-transmitted", 
                                            "Adoption",  "Sibling",  
                                            "Sibling - total", "Sibling - population",
                                            "Trios"),
                          variable.name = "design",
                          value.name    = "value")





#1. Figure Compare different Sibling designs ######

sib <- total_results_long[which(total_results_long$design == "Truth" | 
                                   total_results_long$design == "Sibling" | 
                                   #total_results_long$design == "Sibling - total"| 
                                  total_results_long$design == "Sibling - population"),]

sib <- sib[which(sib$simulation == "ind" | sib$simulation == "noind"),]
sib$simulation <- factor(sib$simulation, levels = c("noind", "ind"))

# Change value Truth 

test <- sib %>%
  group_by(simulation, effect, design) %>%
  summarise(mean = mean(value))

test
sib$value[which(sib$design == "Truth" & sib$effect == "direct")] <- 0.311
sib$value[sib$design == "Truth" & sib$effect == "indirect" & sib$simulation == "noind"] <- 0
sib$value[sib$design == "Truth" & sib$effect == "indirect" & sib$simulation == "ind"] <- 0.139
sib$design <- as.character(sib$design)
sib$design[which(sib$design == "Sibling")] <- "Sibling - between"

sibfig <- sib[which(sib$design != "Truth"),] 
sibfig_truth <- test[which(test$design == "Truth"),] 
sibfig$effect <- factor(sibfig$effect, levels = c("direct", "indirect"),
                        labels = c("Direct effect", "Parental indirect effect")
)
sibfig_truth$effect <- factor(sibfig_truth$effect, levels = c("direct", "indirect"),
                        labels = c("Direct effect", "Parental indirect effect")
)


ggplot(data=sibfig, aes(x=simulation, y=value, fill=design))+ 
  #geom_violin()+ 
  #geom_boxplot(width=0.1)+
  geom_boxplot(outlier.shape=NA)+ 
  facet_wrap(~ effect)+ 
  theme_minimal(base_size = 22)+
  theme(panel.grid=element_blank())+ 
  ylab(element_blank()) + 
  xlab("Effect in simulated data")+
  ggtitle("Estimated direct and parental indirect genetic effects")+ #, subtitle= "100 * N = 20000" genetic effects"
  geom_segment(data=sibfig_truth, # https://stackoverflow.com/questions/45617136/combine-ggplot-facet-wrap-with-geom-segment-to-draw-mean-line-in-scatterplot
                aes(x=as.numeric(simulation) - 0.45,y=mean, xend=as.numeric(simulation) +0.45, yend=mean, color="red"),
                inherit.aes=FALSE,
                linetype= "longdash", #blank", "solid", "dashed", "dotted", "dotdash", "longdash", and "twodash".
                size= 1.2) + 
  scale_x_discrete(labels=c("noind" = "No indirect\neffect", "ind" = "Parental\nindirect\neffect"))+
  #theme(axis.text.x = element_text(angle = 45, hjust=1))+ #rotate labels
  labs(fill="Designs", color="Simulated") + 
  scale_color_manual(name=element_blank(), values= "red", labels= "True parental\nindirect effect") + 
    guides(
      color = guide_legend(order = 0),
      fill = guide_legend(order = 1)
    )
  #theme(legend.title = element_blank())

#2. Supp Fig.  Figure Compare core designs for all conditions ######

all <- total_results_long[which(total_results_long$design == "Truth" | 
                                  total_results_long$design == "Sibling - population"| 
                                  total_results_long$design == "Adoption" | 
                                  total_results_long$design == "Non-transmitted"| 
                                  total_results_long$design == "Trios"
                                  ),]
all$simulation <- factor(all$simulation, levels = c("noind", "ind", "prenatal", "sib", "ind_sib", "AM", "AM_ind", "popstrat", "popstrat_ind"))
all$design <- factor(all$design, levels = c("Truth", "Sibling - population", "Adoption", "Non-transmitted", "Trios"))
all$effect <- factor(all$effect, levels = c("direct", "indirect"),
                     labels = c("Direct effect", "Parental indirect effect"))

# Change value Truth 

test <- all %>%
  group_by(simulation, effect, design) %>%
  summarise(mean = mean(value))

test[test$design == "Truth",] 
#CHECK VALUES
all$value[which(all$design == "Truth" & all$effect == "direct")] <- 0.311
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "noind"] <-0
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "ind"] <- 0.139
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "prenatal"] <- 0.139
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "sib"] <- 0
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "ind_sib"] <- 0.139
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "AM_ind"] <- 0.139
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "AM"] <- 0
all$value[which(all$design == "Truth" & all$effect == "direct" & all$simulation == "popstrat_ind")] <- 0.231
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "popstrat_ind"] <- 0.103
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "popstrat"] <- 0

allfig <- all[which(all$design != "Truth"),] 
allfig_truth <- test[which(test$design == "Truth"),] 


ggplot(data=allfig, aes(x=simulation, y=value, fill=design))+ 
  #geom_violin()+ 
  #geom_boxplot(width=0.1)+
  geom_boxplot(outlier.shape=NA)+ 
  facet_wrap(~ effect)+ 
  theme_minimal(base_size = 22)+
  theme(panel.grid=element_blank())+ 
  ylab(element_blank()) + 
  xlab("Effect in simulated data")+
  ggtitle("Estimated direct and parental indirect genetic effects")+ #, subtitle= "100 * N = 20000" genetic effects"
  # scale_x_discrete(labels=c("noind" = "No indirect\neffect", "ind" = "Parental\nindirect\neffect", 
  #                           "prenatal" = "Pre\nand post-\nnatal\nparental\nindirect\neffect", 
  #                           "sib"= "Sibling\nindirect\neffect", "ind_sib"= "Sibling\nand\nparental\nindirect\neffect", 
  #                           "AM"= "Assortative\nmating", "AM_ind"= "Assortative\nmating\nand\nparental\nindirect\neffect"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ #rotate labels
  scale_x_discrete(labels=c("noind" = "No indirect effect", "ind" = "Parental\nindirect effect", 
                            "prenatal" = "Pre and post-natal\nparental indirect effect", 
                            "sib"= "Sibling\nindirect effect", "ind_sib"= "Sibling and\nparental indirect effect", 
                            "AM"= "Assortative mating","AM_ind"= "Assortative mating and\nparental indirect effect", 
                            "popstrat"= "Population stratification", "popstrat_ind"= "Population stratification and\nparental indirect effect"))+
  labs(color="Simulated", fill="Designs" )+ 
  geom_segment(data=allfig_truth, # https://stackoverflow.com/questions/45617136/combine-ggplot-facet-wrap-with-geom-segment-to-draw-mean-line-in-scatterplot
             aes(x=as.numeric(simulation) - 0.45,y=mean, xend=as.numeric(simulation) +0.45, yend=mean, color="red"),
             inherit.aes=FALSE,
             linetype= "longdash", #blank", "solid", "dashed", "dotted", "dotdash", "longdash", and "twodash".
             size= 1.2) + 
  scale_color_manual(name=element_blank(), values= "red", labels= "True parental\nindirect effect") + 
  guides(
    color = guide_legend(order = 0),
    fill = guide_legend(order = 1)
  )



#3. Figure 3 short ###### 

all <- total_results_long[which(total_results_long$design == "Truth" | 
                                  total_results_long$design == "Sibling - population"| 
                                  total_results_long$design == "Adoption" | 
                                  total_results_long$design == "Non-transmitted"
),]
all$simulation <- factor(all$simulation, levels = c("noind", "ind", "prenatal", "sib", "ind_sib", "AM", "AM_ind","popstrat", "popstrat_ind"))
all$design <- factor(all$design, levels = c("Truth", "Sibling - population", "Adoption", "Non-transmitted"), labels = c("Truth", "Sibling", "Adoption", "Non-transmitted"))
all$effect <- factor(all$effect, levels = c("direct", "indirect"),
                     labels = c("Direct effect", "Parental indirect effect"))

# Change value Truth 

test <- all %>%
  group_by(simulation, effect, design) %>%
  summarise(mean = mean(value))

test[test$design == "Truth",]
all$value[which(all$design == "Truth" & all$effect == "direct")] <- 0.311
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "noind"] <- 0
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "ind"] <- 0.139
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "prenatal"] <- 0.139
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "sib"] <- 0
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "ind_sib"] <- 0.139
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "AM_ind"] <- 0.139
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "AM"] <- 0
all$value[which(all$design == "Truth" & all$effect == "direct" & all$simulation == "popstrat_ind")] <- 0.231
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "popstrat_ind"] <- 0.103
all$value[all$design == "Truth" & all$effect == "indirect" & all$simulation == "popstrat"] <- 0

allfig <- all[which(all$design != "Truth"),] 
allfig_truth <- test[which(test$design == "Truth"),] 

allfig <- allfig[which(allfig$effect != "Direct effect"),] 
allfig_truth <- allfig_truth[which(allfig_truth$effect != "Direct effect"),] 


p <- ggplot(data=allfig, aes(x=simulation, y=value, fill=design))+ 
  #geom_violin()+ 
  #geom_boxplot(width=0.1)+
  geom_boxplot(outlier.shape=NA)+
  #facet_wrap(~ effect)+ 
  theme_minimal(base_size = 18)+
  theme(panel.grid=element_blank())+ 
  ylab("Estimated parental indirect genetic effects") + 
  #coord_cartesian(ylim=c(-0.2, 0.5))+
  xlab("Components and biases included in simulated data")+
  #ggtitle("Estimated parental indirect genetic effects")+ #, subtitle= "100 * N = 20000" genetic effects"
  # scale_x_discrete(labels=c("noind" = "No indirect\neffect", "ind" = "Parental\nindirect\neffect", 
  #                           "prenatal" = "Pre\nand post-\nnatal\nparental\nindirect\neffect", 
  #                           "sib"= "Sibling\nindirect\neffect", "ind_sib"= "Sibling\nand\nparental\nindirect\neffect", 
  #                           "AM"= "Assortative\nmating", "AM_ind"= "Assortative\nmating\nand\nparental\nindirect\neffect"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+ #rotate labels
  scale_x_discrete(labels=c("noind" = "1. No indirect effect", "ind" = "2. Parental indirect effect", 
                            "prenatal" = "3. Pre and post-natal\nparental indirect effect", 
                            "sib"= "4. Sibling indirect effect", "ind_sib"= "5. Sibling and\nparental indirect effect", 
                            "AM"= "6. Assortative mating", "AM_ind"= "7. Assortative mating and\nparental indirect effect",
                            "popstrat"= "8. Population stratification", "popstrat_ind"= "9. Population stratification and\nparental indirect effect"))+
  labs(color="Simulated", fill="Designs" )+ 
  geom_segment(data=allfig_truth, # https://stackoverflow.com/questions/45617136/combine-ggplot-facet-wrap-with-geom-segment-to-draw-mean-line-in-scatterplot
               aes(x=as.numeric(simulation) - 0.45,y=mean, xend=as.numeric(simulation) +0.45, yend=mean, color="red"),
               inherit.aes=FALSE,
               linetype= "longdash", #blank", "solid", "dashed", "dotted", "dotdash", "longdash", and "twodash".
               size= 1.2) + 
  scale_color_manual(name=element_blank(), values= "red", labels= "True postnatal\nparental indirect effect") + 
  guides(
    color = guide_legend(order = 0),
    fill = guide_legend(order = 1)
  )
#theme(legend.title = element_blank())

p


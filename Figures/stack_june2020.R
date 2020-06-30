#rosa
#non cog stack plots june 2020
setwd("~/Dropbox/ROSA/Perline_cognoncog/figures")
library(dplyr)
library(ggplot2)
library(reshape2)
library(readxl)
library(tidyr)

#1. plot % of total effect of Cog and Noncog polygenic score  accounted for by direct and indirect effects. 
#WTIHIN EACH SAMPLE, ACROSS OUTCOMES AND METHODS
#NB Within-cohort results were aggregated across outcome measures/ages and methods.
stack_1_rosa <- read_excel("~/Dropbox/ROSA/Perline_cognoncog/figures/stack_1_rosa.xlsx")
stack_1_rosa$PRSxEffect <- paste(stack_1_rosa$PRS, stack_1_rosa$Effect)
#make variables into ordered factors
stack_1_rosa$Sample <- factor(stack_1_rosa$Sample, levels=c("UKB", "TEDS","NTR"))
stack_1_rosa$Effect <- factor(stack_1_rosa$Effect, levels=c( "Indirect","Direct"))

ggplot(stack_1_rosa,aes(y=Value, x=Sample, fill=factor(PRSxEffect))) + 
  theme_bw()+
  #coord_flip(ylim = c(0,1))+
  geom_bar(stat = "identity")+
  #theme(legend.position="none")+
  theme(legend.title=element_blank())+
  theme(axis.ticks.y=element_blank(),
        axis.title.y=element_blank())+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  facet_grid(~PRS) + 
  theme(strip.text.x = element_text(size = 16))+
  scale_fill_manual(values=c("blue", "light blue","orange", "yellow"))
  
#2. plot % of total effect of Cog and Noncog polygenic score  accounted for by direct and indirect effects. 
#ACROSS METHODS, WITHIN SAMPLE/OUTCOME
stack_2_rosa <- read_excel("~/Dropbox/ROSA/Perline_cognoncog/figures/stack_2_rosa.xlsx")
stack_2_rosa$PRSxEffect <- paste(stack_2_rosa$PRS, stack_2_rosa$Effect)
stack_2_rosa$Effect <- factor(stack_2_rosa$Effect, levels=c( "Indirect","Direct"))
stack_2_rosa$Method <- factor(stack_2_rosa$Method, levels=c( "Sibling","Trio","Adoption"))

ggplot(stack_2_rosa,aes(y=Value, x=Method, fill=factor(PRSxEffect))) + 
  theme_bw()+
  #coord_flip(ylim = c(0,1))+
  geom_bar(stat = "identity")+
  #theme(legend.position="none")+
  theme(legend.title=element_blank())+
  theme(axis.ticks.y=element_blank(),
        axis.title.y=element_blank())+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  facet_grid(~PRS) + 
  theme(strip.text.x = element_text(size = 16))+
  scale_fill_manual(values=c("blue", "light blue","orange", "yellow"))


#2. plot % of total effect of Cog and Noncog polygenic score  accounted for by direct and indirect effects. 
#ACROSS METHODS, WITHIN SAMPLE/OUTCOME
stack_3_rosa <- read_excel("~/Dropbox/ROSA/Perline_cognoncog/figures/stack_3_rosa.xlsx")
stack_3_rosa$PRSxEffect <- paste(stack_3_rosa$PRS, stack_3_rosa$Effect)
stack_3_rosa$Effect <- factor(stack_3_rosa$Effect, levels=c( "Indirect","Direct"))
stack_3_rosa$Outcome <- factor(stack_3_rosa$Outcome, levels=c( "UKB EA","TEDS 12","TEDS 16","NTR CITO","NTR EA"))

ggplot(stack_3_rosa,aes(y=Value, x=Outcome, fill=factor(PRSxEffect))) + 
  theme_bw()+
  #coord_flip(ylim = c(0,1))+
  geom_bar(stat = "identity")+
  #theme(legend.position="none")+
  theme(legend.title=element_blank())+
  theme(axis.ticks.y=element_blank(),
        axis.title.y=element_blank())+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_text(size = 6.5)) +
  theme(axis.text.y = element_text(size = 16)) +
  facet_grid(~PRS) + 
  theme(strip.text.x = element_text(size = 16))+
  scale_fill_manual(values=c("blue", "light blue","orange", "yellow"))





# rosa cheesman & Perline Demange 
# july 2020
# plotting noncog indirect genetic effect results



# Libraries #############################

library(dplyr)
library(ggplot2)
library(reshape2)
library(psych)
library(readr)
library(stringr)


# 1. (Fig 2.A) make a stacked bar plot for overall metaanalysed results ###########################
all <- read_csv("../Meta-analysis/meta_all_20210519.csv")
all$Effect <- factor(all$Effect, levels=c( "Direct","Indirect"))
all$PGS <- factor(all$PGS, levels=c("NonCog","Cog"))
all$PGSxEffect <- paste(all$PGS, all$Effect)
all$PGSxEffect <- factor(all$PGSxEffect, levels=c("NonCog Direct", "Cog Direct", "NonCog Indirect","Cog Indirect"))

all$blob<-ifelse(all$PGSxEffect=="NonCog Direct",1,
                 ifelse(all$PGSxEffect=="NonCog Indirect",2,
                        ifelse(all$PGSxEffect=="Cog Direct",3,
                               ifelse(all$PGSxEffect=="Cog Indirect",4,NA))))

all <- all[order(all$blob),]


#for SEs, need the direct effect to be cumulation of direct + indirect
#all$est3<-c( (0.13302074+0.08679965),0.08679965,(0.15872786+0.10410311),0.10410311)
all$est3<-c( (all$estimate[1]+all$estimate[2]),all$estimate[2],(all$estimate[3]+all$estimate[4]),all$estimate[4])



# Function to increase the vertical spacing between legends keys here but doesn't work
# no possiblity to increase spacing in vertical legends 
# https://stackoverflow.com/questions/11366964/is-there-a-way-to-change-the-spacing-between-legend-items-in-ggplot2

tiff("Fig2A_20210615.tiff")
plot_A <- ggplot(all,aes(y=estimate, x=PGS, fill=factor(PGSxEffect))) + 
  theme_bw(base_size=14)+
  #coord_flip(ylim = c(0,1))+
  geom_bar(stat = "identity", width=0.9)+
  theme(legend.title=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, margin=margin(t=0,r=2, b=0, l=0, unit="lines")),
        legend.key.size = unit(3, "lines"), 
        legend.text = element_text(size=14, margin = margin(t=1, b=1, unit="lines")))+
  #scale_fill_manual(values=c( "orange","yellow","blue","light blue"))   +
  scale_fill_manual(values=c(  "#ff9900","#0033cc","#ffcc66", "#3399ff"))   +
  #ggtitle("A. Meta-analytic results") +
  labs(y = "Effect of polygenic score on educational outcome", x = " ") +
  geom_errorbar(aes(ymin = est3-1.96*se, ymax = est3+1.96*se), width = 0.06) +
  scale_y_continuous(breaks = seq(0, 0.30, by = 0.05))
plot_A
dev.off()

# # 2. stacked bar plot for results per method metanalaysis results #########
# method <- read_csv("../Meta-analysis/meta_method_20210504.csv")
# 
# #add adoption results (not meta analysed as not multuple cohorts or outcomes)
# adopt <- read.table("../UKB/Adoptees/summary_mean_CI_adoption_UKB_20200529.csv")
# 
# adoption<- as.data.frame(cbind(c(replicate(2, "NonCog"),replicate(2, "Cog")),
#                                   rep(c('Direct', 'Indirect'), 2),
#                                   rep(c('adopt'), 4),
#                                c(adopt$direct_NonCog[1],adopt$indirect_NonCog[1],adopt$direct_Cog[1],adopt$indirect_Cog[1]),#noncog direct indirect, cog same
#                                c(adopt$direct_NonCog[4],adopt$indirect_NonCog[4],adopt$direct_Cog[4],adopt$indirect_Cog[4]),
#                                NA,NA,NA,NA))
# colnames(adoption)=colnames(method)
# method<-rbind(method,adoption)
# 
# method$Effect <- factor(method$Effect, levels=c( "Direct","Indirect"))
# method$PGS <- factor(method$PGS, levels=c("NonCog","Cog"))
# #method$Sample <- factor(method$Sample, levels=c( "sib","adopt","trio"))
# method$PGSxEffect <- paste(method$PGS, method$Effect)
# method$PGSxEffect <- factor(method$PGSxEffect, levels=c("NonCog Direct", "NonCog Indirect", "Cog Direct","Cog Indirect"))
# method$estimate<-as.numeric(method$estimate)
# #method$est3<-c( (0.13302074+0.08679965),0.08679965,(0.15872786+0.10410311),0.10410311)
# 
# ggplot(method,aes(y=estimate, x=factor(Sample), fill=factor(PGSxEffect))) + 
#   theme_bw()+
#   geom_bar(stat = "identity")+
#   theme(legend.title=element_blank())+
#   theme(axis.ticks.y=element_blank(),
#         axis.title.y=element_blank())+
#   theme(axis.title.x=element_blank())+
#   theme(axis.text.y = element_text(size = 10)) +
#   theme(strip.text.x = element_text(size = 14))+
#   facet_grid(~PGS) + 
#   theme(strip.background = element_blank(),strip.text.y = element_blank())+
#   scale_fill_manual(values=c( "orange","yellow","blue","light blue"))   +
#  # geom_errorbar(aes(ymin = est3-se, ymax = est3+se), width = .05) +
#   scale_y_continuous(breaks = seq(0, 0.40, by = 0.05))
# 
# 
# 
# ## 3.1 Compare methods with all data + metanaalytic results #########
# 
# all <- read.table("../dataforfigures_20210503.csv")
# 
# split2 <- str_split_fixed(all$Methods, "_", 2)
# colnames(split2) <- c("Method", "Sample")
# head(split2)
# split2 <- print.data.frame(as.data.frame(split2), quote=FALSE)
# 
# all <- cbind(all, split2)
# all$Samples <- paste(all$Sample, all$pheno, sep="_")
# head(all)
# 
# 
# plotdata=all[,c("original", "leftCI", "rightCI","Type","PRS", "Method","Samples")]
# plotdata <- plotdata[which(plotdata$Type == "indirect"| plotdata$Type == "direct"),]
# head(plotdata)
# str(plotdata)
# 
# 
# data<-plotdata
# meta  <- read.table("meta_method_forfig.csv", header=T, sep=",")
# data <- rbind(data, meta)
# head(data)
# 
# 
# data$Type<- factor(data$Type, levels=c("direct","indirect"))
# data$Method <- factor(data$Method, levels=c("Siblings", "Adoption", "Trios"))
# data$Samples <- factor(data$Samples, levels=c("UKB_EA", "TEDS_GCSE", "TEDS_12yo", "NTR_EA", "NTR_CITO", "Meta"))
# 
# 
# #subset cog so easier color
# datacog <- data[data$PRS == "Cog",]
# datanoncog <- data[data$PRS == "NonCog",]
# 
# plotcog <- ggplot(datacog, aes(y=original, x=Samples)) +
#   geom_point(aes(x=Samples, y=original, shape= Samples, color = Type), size=5,position = position_dodge(.5))+#,position = position_dodge(4))+
#   geom_errorbar(aes(x=Samples,ymin=rightCI, ymax=leftCI, color=Type),  width=.9,position = position_dodge(.5)) +#, position=position_dodge(4)) + 
#   theme_bw()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         panel.background=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   labs(y = "Estimated effect of polygenic score on educational outcome", x = " ") +
#   facet_grid(~Method+PRS, scales="free", space="free") +
#   scale_color_manual(values=c("blue", "light blue")) +
#   scale_shape_manual(values=c(21, 24, 25, 22, 23, 18))
# 
# plotnoncog <- ggplot(datanoncog, aes(y=original, x=Samples)) +
#   geom_point(aes(x=Samples, y=original, shape= Samples, color = Type), size=5,position = position_dodge(.5))+#,position = position_dodge(4))+
#   geom_errorbar(aes(x=Samples,ymin=rightCI, ymax=leftCI, color=Type),  width=.9,position = position_dodge(.5)) +#, position=position_dodge(4)) + 
#   theme_bw()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(), 
#         panel.background=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   labs(y = "Estimated effect of polygenic score on educational outcome", x = " ") +
#   facet_grid(~Method+PRS, scales="free", space="free") +
#   scale_color_manual(values=c("orange", "yellow")) +
#   scale_shape_manual(values=c(21, 24, 25, 22, 23, 18))
# 
# png("plotcog.png")
# plotcog
# dev.off()
# 
# 
# png("plotnoncog.png")
# plotnoncog
# dev.off()
# 
# 
# 
# # 3.2 Compare method meta-analysis all results ######################################
# 
# method <- read_csv("../Meta-analysis/meta_method_20210503.csv")
# 
# all <- read.table("../dataforfigures_20210503.csv", stringsAsFactors = F)
# split2 <- str_split_fixed(all$Methods, "_", 2)
# colnames(split2) <- c("Method", "Sample")
# head(split2)
# split2 <- print.data.frame(as.data.frame(split2), quote=FALSE)
# 
# all <- cbind(all, split2)
# all$Samples <- paste(all$Sample, all$pheno, sep="_")
# head(all)
# 
# adoption <- all[all$Method == "Adoption", ]
# adoption <- adoption[adoption$Type == "direct" | adoption$Type == "indirect", ]
# adoption$Type[adoption$Type == "direct"] <- "Direct"
# adoption$Type[adoption$Type == "indirect"] <- "Indirect"
# adoption <- adoption[, c("PRS", "Type", "Method", "original", "leftCI", "rightCI")]
# adoption
# method <- method[, c("PGS", "Effect", "Sample", "estimate", "ci.lb", "ci.ub")]
# method
# colnames(adoption) <- c("PGS", "Effect", "Method", "estimate", "ci.lb", "ci.ub")
# colnames(method) <- colnames(adoption)
# method <- rbind(method, adoption)
# 
# method$Method[method$Method == "sib"] <- "Siblings"
# method$Method[method$Method == "trio"] <- "Non-transmitted"
# method$Method <- factor(method$Method, c("Siblings", "Adoption", "Non-transmitted"))
#                         
# plot <- 
#   ggplot(method, aes(y=estimate, x=Effect)) +
#   geom_point(aes(x=Effect, y=estimate, color = interaction(Effect, PGS)), size=5,position = position_dodge(.5))+#,position = position_dodge(4))+
#   geom_errorbar(aes(x=Effect,ymin=ci.lb, ymax=ci.ub, color = interaction(Effect, PGS)),  width=.9,position = position_dodge(.5)) +#, position=position_dodge(4)) + 
#   theme_bw()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         panel.background=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   labs(y = "Estimated effect of polygenic score on educational outcome", x = " ") +
#   facet_grid(~Method, scales="free", space="free") +
#   scale_color_manual(values=c("blue",  "orange", "light blue","yellow")) +
#   ggtitle("Per Method, including all outcomes") +
#   scale_shape_manual(values=c(21, 24, 25, 22, 23, 18))
# plot


# (Fig 2.C) 3.3 Compare method only including EA ##############

method <- read_csv("../Meta-analysis/meta_methodEA_20210519.csv")

all <- read.table("../dataforfigures_20210519.csv", stringsAsFactors = F, header=T)
split2 <- str_split_fixed(all$Methods, "_", 2)
colnames(split2) <- c("Method", "Sample")
head(split2)
split2 <- print.data.frame(as.data.frame(split2), quote=FALSE)
all <- cbind(all, split2)
all$Samples <- paste(all$Sample, all$pheno, sep="_")
head(all)

adoption <- all[all$Method == "Adoption", ]
adoption <- adoption[adoption$Type == "direct" | adoption$Type == "indirect", ]
adoption$Type[adoption$Type == "direct"] <- "Direct"
adoption$Type[adoption$Type == "indirect"] <- "Indirect"
adoption <- adoption[, c("PRS", "Type", "Method", "original", "leftCI", "rightCI")]
adoption

trios <- all[all$Method == "Trios", ]
trios <- trios[trios$pheno == "EA", ]
trios <- trios[trios$Type == "direct" | trios$Type == "indirect", ]
trios$Type[trios$Type == "direct"] <- "Direct"
trios$Type[trios$Type == "indirect"] <- "Indirect"
trios <- trios[, c("PRS", "Type", "Method", "original", "leftCI", "rightCI")]
trios


method$Sample <- "Siblings"
method <- method[, c("PGS", "Effect", "Sample", "estimate", "ci.lb", "ci.ub")]

colnames(adoption) <- c("PGS", "Effect", "Method", "estimate", "ci.lb", "ci.ub")
colnames(method) <- colnames(adoption)
colnames(trios) <- colnames(adoption)
method <- rbind(method, adoption, trios)

method$Method[method$Method == "Trios"] <- "Non-transmitted"
method$Method <- factor(method$Method, c("Siblings", "Adoption", "Non-transmitted"))
method$PGS <- factor(method$PGS, c("NonCog", "Cog"))

tiff("Fig2C_20210615.tiff")
plot_C <- 
  ggplot(method, aes(y=estimate, x=Effect)) +
  geom_point(aes(x=Effect, y=estimate, color = interaction(Effect, PGS)), size=7,position = position_dodge(.6))+#,position = position_dodge(4))+
  geom_errorbar(aes(x=Effect,ymin=ci.lb, ymax=ci.ub, color = interaction(Effect, PGS)),  width=.6,size= 1, position = position_dodge(.6)) +#, position=position_dodge(4)) + 
  theme_bw(base_size = 14)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size = 12),
        #axis.ticks.x=element_blank(),
        axis.text.y=element_text(size = 12),
        strip.text.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14, margin=margin(t=0,r=2, b=0, l=0, unit="lines")),
        strip.background=element_rect(fill="white"),
        panel.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position= "none")+
  labs(y = "Effect of polygenic score on educational outcome", x = " ") +
  facet_grid(~Method, scales="free", space="free") +
  scale_color_manual(values=c( "#ff9900","#ffcc66", "#0033cc", "#3399ff")) +
  #ggtitle("C. Educational attainment by method") +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=1.5, linetype="dotted")+
  scale_shape_manual(values=c(21, 24, 25, 22, 23, 18)) 
plot_C
dev.off()

# # 3.4 Compare method only including EA not meta ##############

all <- read.table("../dataforfigures_20210519.csv", stringsAsFactors = FALSE, header=T)
split2 <- str_split_fixed(all$Methods, "_", 2)
colnames(split2) <- c("Method", "Sample")
head(split2)
split2 <- print.data.frame(as.data.frame(split2), quote=FALSE)
all <- cbind(all, split2)
all$Samples <- paste(all$Sample, all$pheno, sep="_")
all$Method <- as.character(all$Method)
all$Sample <- as.character(all$Sample)
str(all)

adoption <- all[all$Method == "Adoption", ]
adoption <- adoption[adoption$Type == "direct" | adoption$Type == "indirect", ]
adoption$Type[adoption$Type == "direct"] <- "Direct"
adoption$Type[adoption$Type == "indirect"] <- "Indirect"
adoption <- adoption[, c("PRS", "Type", "Method", "Sample", "original", "leftCI", "rightCI")]
adoption

trios <- all[all$Method == "Trios", ]
trios <- trios[trios$pheno == "EA", ]
trios <- trios[trios$Type == "direct" | trios$Type == "indirect", ]
trios$Type[trios$Type == "direct"] <- "Direct"
trios$Type[trios$Type == "indirect"] <- "Indirect"
trios <- trios[, c("PRS", "Type", "Method", "Sample", "original", "leftCI", "rightCI")]
trios

sib <- all[all$Method == "Siblings", ]
sib <- sib[sib$pheno == "EA", ]
sib <- sib[sib$Type == "direct" | sib$Type == "indirect", ]
sib$Type[sib$Type == "direct"] <- "Direct"
sib$Type[sib$Type == "indirect"] <- "Indirect"
sib <- sib[, c("PRS", "Type", "Method", "Sample", "original", "leftCI", "rightCI")]
sib

colnames(adoption) <- c("PGS", "Effect", "Method", "Sample", "estimate", "ci.lb", "ci.ub")
colnames(sib) <- colnames(adoption)
colnames(trios) <- colnames(adoption)
method <- rbind(sib, adoption, trios)

method$Method[method$Method == "Trios"] <- "Non-transmitted"
method$Method <- factor(method$Method, c("Siblings", "Adoption", "Non-transmitted"))
method$PGS <- factor(method$PGS, c("NonCog", "Cog"))

plot <-
  ggplot(method, aes(y=estimate, x=Effect)) +
  geom_point(aes(x=Effect, y=estimate, color = interaction(Effect, PGS)), size=5,position = position_dodge(.5))+#,position = position_dodge(4))+
  geom_errorbar(aes(x=Effect,ymin=ci.lb, ymax=ci.ub, color = interaction(Effect, PGS)),  width=.6,position = position_dodge(.5)) +#, position=position_dodge(4)) +
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(y = "Estimated effect of polygenic score on educational outcome", x = " ") +
  facet_grid(~Sample+Method, scales="free", space="free") +
  scale_color_manual(values=c( "#ff9900","#ffcc66", "#0033cc", "#3399ff")) +
  ggtitle("Per Method, Educational Attainment") +
  scale_shape_manual(values=c(21, 24, 25, 22, 23, 18))
plot
# 
# 
# # 4.2 Compare child outcome all results  only siblings method ###################
# 
# all <- read.table("../dataforfigures_20210503.csv", stringsAsFactors = FALSE)
# split2 <- str_split_fixed(all$Methods, "_", 2)
# colnames(split2) <- c("Method", "Sample")
# head(split2)
# split2 <- print.data.frame(as.data.frame(split2), quote=FALSE)
# all <- cbind(all, split2)
# all$Samples <- paste(all$Sample, all$pheno, sep="_")
# all$Method <- as.character(all$Method)
# all$Sample <- as.character(all$Sample)
# str(all)
# 
# child <- all[!all$pheno == "EA", ]
# child <- child[child$Type == "direct" | child$Type == "indirect", ]
# child <- child[child$Method == "Siblings",]
# child$PRS <- factor(child$PRS, c("NonCog", "Cog"))
# plot <- 
#   ggplot(child, aes(y=original, x=Samples)) +
#   geom_point(aes(x=Samples, y=original, color = interaction(Type, PRS)), size=5,position = position_dodge(.5))+#,position = position_dodge(4))+
#   geom_errorbar(aes(x=Samples,ymin=leftCI, ymax=rightCI, color = interaction(Type, PRS)),  width=.9,position = position_dodge(.5)) +#, position=position_dodge(4)) + 
#   theme_bw()+
#   theme(axis.title.x=element_blank(),
#         panel.background=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   labs(y = "Estimated effect of polygenic score on educational outcome", x = " ") +
#   facet_grid(~Sample, scales="free", space="free") +
#   scale_color_manual(values=c("#ff9900", "#ffcc00", "#0033cc",  "#3399ff")) +
#   ggtitle("Child outcome in sibling method") +
#   scale_shape_manual(values=c(21, 24, 25, 22, 23, 18))
# plot
# 
# # 4.3 Compare child outcome 12yo  only siblings method ###################
# 
# all <- read.table("../dataforfigures_20210503.csv", stringsAsFactors = FALSE)
# split2 <- str_split_fixed(all$Methods, "_", 2)
# colnames(split2) <- c("Method", "Sample")
# head(split2)
# split2 <- print.data.frame(as.data.frame(split2), quote=FALSE)
# all <- cbind(all, split2)
# all$Samples <- paste(all$Sample, all$pheno, sep="_")
# all$Method <- as.character(all$Method)
# all$Sample <- as.character(all$Sample)
# str(all)
# 
# child <- all[!all$pheno == "EA" & !all$pheno == "GCSE" , ]
# child <- child[child$Type == "direct" | child$Type == "indirect", ]
# child <- child[child$Method == "Siblings",]
# child$PRS <- factor(child$PRS, c("NonCog", "Cog"))
# plot <- 
#   ggplot(child, aes(y=original, x=Type)) +
#   geom_point(aes(x=Type, y=original, color = interaction(Type, PRS)), size=5,position = position_dodge(.5))+#,position = position_dodge(4))+
#   geom_errorbar(aes(x=Type,ymin=leftCI, ymax=rightCI, color = interaction(Type, PRS)),  width=.9,position = position_dodge(.5)) +#, position=position_dodge(4)) + 
#   theme_bw()+
#   labs(y = "Estimated effect of polygenic score on educational outcome", x = " ") +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         panel.background=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   facet_grid(~Sample, scales="free", space="free") +
#   scale_color_manual(values=c("#ff9900", "#ffcc00", "#0033cc",  "#3399ff")) +
#   ggtitle("Child outcome in sibling method") +
#   geom_hline(yintercept=0, linetype="dashed")+
#   scale_shape_manual(values=c(21, 24, 25, 22, 23, 18))
# plot

# (Fig 2.B) 4.4 Compare 12yo and EA  only siblings method ###################

all <- read.table("../dataforfigures_20210519.csv", stringsAsFactors = FALSE, header=T)
split2 <- str_split_fixed(all$Methods, "_", 2)
colnames(split2) <- c("Method", "Sample")
head(split2)
split2 <- print.data.frame(as.data.frame(split2), quote=FALSE)
all <- cbind(all, split2)
all$Samples <- paste(all$Sample, all$pheno, sep="_")
all$Method <- as.character(all$Method)
all$Sample <- as.character(all$Sample)
str(all)

child <- all[!all$pheno == "GCSE" , ]
child <- child[child$Type == "direct" | child$Type == "indirect", ]
child <- child[child$Method == "Siblings",]
child$Type[child$Type == "direct"] <- "Direct"
child$Type[child$Type == "indirect"] <- "Indirect"
child$Type <- factor(child$Type, c("Direct", "Indirect"))
child$PRS <- factor(child$PRS, c("NonCog", "Cog"))
child$pheno[child$pheno == "CITO"] <- "12yo"
child$pheno[child$pheno == "12yo"] <- "Age 12 achievement"
child$pheno[child$pheno == "EA"] <- "Adult attainment"
child$pheno <- factor(child$pheno, c("Age 12 achievement", "Adult attainment"))
# child$Sample[child$Sample == "TEDS"] <- "UK"
# child$Sample[child$Sample == "UKB"] <- "UK"


library(grid)
library(gtable)

#https://stackoverflow.com/questions/40732543/seeking-workaround-for-gtable-add-grob-code-broken-by-ggplot-2-2-0
OverlappingStripLabels = function(plot) {
  
  # Get the ggplot grob
  g = ggplotGrob(plot)
  
  ### Collect some information about the strips from the plot
  # Get a list of strips
  strip = lapply(grep("strip-t", g$layout$name), function(x) {g$grobs[[x]]})
  
  
  # Number of strips
  NumberOfStrips = sum(grepl(pattern = "strip-t", g$layout$name))
  
  # Number of rows
  NumberOfRows = length(strip[[1]])
  
  # Panel spacing and its unit
  plot_theme <- function(p) {
    plyr::defaults(p$theme, theme_get())
  }
  PanelSpacing = plot_theme(plot)$panel.spacing
  unit = attr(PanelSpacing, "unit")
  
  # Map the boundaries of the new strips
  Nlabel = vector("list", NumberOfRows)
  map = vector("list", NumberOfRows)
  for(i in 1:NumberOfRows) {
    
    for(j in 1:NumberOfStrips) {
      Nlabel[[i]][j] = getGrob(grid.force(strip[[j]][i]), gPath("GRID.text"), grep = TRUE)$label
    }
    
    map[[i]][1] = TRUE
    for(j in 2:NumberOfStrips) {
      map[[i]][j] = Nlabel[[i]][j] != Nlabel[[i]][j-1]
    }
  }
  
  
  
  ## Construct gtable to contain the new strip
  # Set the widths of the strips, based on widths of the panels and PanelSpacing
  panel = subset(g$layout, grepl("panel", g$layout$name), l, drop = TRUE)                       
  StripWidth = list()
  for(i in seq_along(panel)) StripWidth[[i]] = unit.c(g$width[panel[i]], PanelSpacing)
  
  newStrip  = gtable(widths = unit.c(unit(unlist(StripWidth), c("null", unit)))[-2*NumberOfStrips], 
                     heights = strip[[1]]$heights)
  
  
  ## Populate the gtable  
  seqLeft = list()
  for(i in 1:NumberOfRows) {  
    Left = which(map[[i]] == TRUE)
    seqLeft[[i]] = if((i-1) < 1) 2*Left - 1 else sort(unique(c(seqLeft[[i-1]], 2*Left - 1))) 
    seqRight = c(seqLeft[[i]][-1] -2, (2*NumberOfStrips-1))
    newStrip = gtable_add_grob(newStrip, lapply(strip[(seqLeft[[i]]+1)/2], `[`, i), t = i, l = seqLeft[[i]], r = seqRight)
  }
  
  #remove second row 
  #newStrip <- newStrip[1]
  
  ## Put the strip into the plot
  # Get the locations of the original strips
  pos = subset(g$layout, grepl("strip-t", g$layout$name), t:r)
  
  ## Use these to position the new strip
  pgNew = gtable_add_grob(g, newStrip, t = unique(pos$t), l = min(pos$l), r = max(pos$r))
  
  return(pgNew)
}



tiff("Fig2B_20210615.tiff")
plot_B <- 
  ggplot(child, aes(y=original, x=Type)) +
  geom_point(aes(x=Type, y=original, color = interaction(Type, PRS)), size=7,
             position = position_dodge(.6))+#,position = position_dodge(4))+
  geom_errorbar(aes(x=Type,ymin=leftCI, ymax=rightCI, color = interaction(Type, PRS)),
                width=.6,size=1, position = position_dodge(.6)) + 
  theme_bw()+
  labs(y = "Effect of polygenic score on educational outcome", x = " ") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size = 12),
        #axis.ticks.x=element_blank(),
        axis.text.y=element_text(size = 12),
        strip.text.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14, margin=margin(t=0,r=2, b=0, l=0, unit="lines")),
        legend.position = "none",
        panel.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(c(0, 0.5, 0) , "lines"),
        strip.background=element_rect(fill="white"))+
  facet_grid(~pheno+Sample, scales="free", space="free") +
  scale_color_manual(values=c("#ff9900", "#ffcc66", "#0033cc",  "#3399ff")) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_vline(xintercept=1.5, linetype="dotted")+
  #ggtitle("B. Sibling method by cohort") +
  scale_shape_manual(values=c(21, 24, 25, 22, 23, 18))
plot_B
grid.newpage()
grid.draw(OverlappingStripLabels(plot_B))
dev.off()

# tiff("Fig2B_20210503_testse.tiff")
# plot_B <- 
#   ggplot(child, aes(y=original, x=Type)) +
#   geom_point(aes(x=Type, y=original, color = interaction(Type, PRS)), size=7,
#              position = position_dodge(.6))+#,position = position_dodge(4))+
#   geom_errorbar(aes(x=Type,ymin=original - 1.96*se, ymax=original + 1.96*se, color = interaction(Type, PRS)),
#                 width=.6,size=1, position = position_dodge(.6)) + 
#   theme_bw()+
#   labs(y = "Estimated effect of polygenic score on educational outcome", x = " ") +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_text(size = 12),
#         #axis.ticks.x=element_blank(),
#         axis.text.y=element_text(size = 12),
#         strip.text.x = element_text(size = 14), 
#         axis.title.y = element_text(size = 14, margin=margin(t=0,r=2, b=0, l=0, unit="lines")),
#         legend.position = "none",
#         panel.background=element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.spacing.x = unit(c(0, 0.5, 0) , "lines"),
#         strip.background=element_rect(fill="white"))+
#   facet_grid(~pheno+Sample, scales="free", space="free") +
#   scale_color_manual(values=c("#ff9900", "#ffcc66", "#0033cc",  "#3399ff")) +
#   #ggtitle("12yo and Adult outcomes, Sibling method") +
#   geom_hline(yintercept=0, linetype="dashed")+
#   geom_vline(xintercept=1.5, linetype="dotted")+
#   scale_shape_manual(values=c(21, 24, 25, 22, 23, 18))
# plot_B
# grid.newpage()
# grid.draw(OverlappingStripLabels(plot_B))
# dev.off()




# (Supp Fig.4) 5. All ratios ###################

all <- read.table("../dataforfigures_20210519.csv", stringsAsFactors = F, header=T)
split2 <- str_split_fixed(all$Methods, "_", 2)
colnames(split2) <- c("Method", "Sample")
head(split2)
split2 <- print.data.frame(as.data.frame(split2), quote=FALSE)
all <- cbind(all, split2)
all$Samples <- paste(all$Sample, all$pheno, sep="_")
all$Method <- as.character(all$Method)
all$Sample <- as.character(all$Sample)
head(all)

all$Type
ratio <- all[all$Type == "ratio_tot",]
ratio$Method[ratio$Method == "Trios"] <- "Non-transmitted"
ratio$Method <- factor(ratio$Method, c("Siblings", "Adoption", "Non-transmitted"))

# save as"Supp_Fig4_20210519.tiff"
plot <- 
  ggplot(ratio, aes(y=original, x=Samples)) +
  geom_point(aes(x=Samples, y=original, color = PRS), size=5,position = position_dodge(.5))+#,position = position_dodge(4))+
  geom_errorbar(aes(x=Samples,ymin=leftCI, ymax=rightCI, color = PRS),  width=.9,position = position_dodge(.5)) +#, position=position_dodge(4)) + 
  theme_bw()+
  theme(axis.title.x=element_blank(),
        panel.background=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(y = "Estimated effect of polygenic score on educational outcome") +
  facet_grid(~Method, scales="free", space="free") +
  scale_color_manual(values=c("blue",  "orange")) +
  geom_hline(yintercept=0.5)+
  ggtitle("Ratio indirect/population") 
plot


# (Supp Fig 3) 6. All total ###################

all <- read.table("../dataforfigures_20210519.csv", stringsAsFactors = F, header=T)
split2 <- str_split_fixed(all$Methods, "_", 2)
colnames(split2) <- c("Method", "Sample")
head(split2)
split2 <- print.data.frame(as.data.frame(split2), quote=FALSE)
all <- cbind(all, split2)
all$Samples <- paste(all$Sample, all$pheno, sep="_")
all$Method <- as.character(all$Method)
all$Sample <- as.character(all$Sample)
head(all)

ratio <- all[all$Type == "total",]
ratio$Method[ratio$Method == "Trios"] <- "Non-transmitted"
ratio$Method <- factor(ratio$Method, c("Siblings", "Adoption", "Non-transmitted"))

#save as "Supp_Fig3_20210519.tiff"
plot <- 
  ggplot(ratio, aes(y=original, x=Samples)) +
  geom_point(aes(x=Samples, y=original, color = PRS), size=5,position = position_dodge(.5))+#,position = position_dodge(4))+
  geom_errorbar(aes(x=Samples,ymin=leftCI, ymax=rightCI, color = PRS),  width=.7,position = position_dodge(.5)) +#, position=position_dodge(4)) + 
  theme_bw()+
  theme(axis.title.x=element_blank(),
        panel.background=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(y = "Estimated effect of polygenic score on educational outcome") +
  facet_grid(~Method, scales="free", space="free") +
  scale_color_manual(values=c("blue",  "orange")) +
  ggtitle("Population genetic effects") 
plot


# 7. Figure with 1, 3.3 and 4.4 -----------
#https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
library(gridExtra)
grid.arrange(plot_A, plot_B, plot_C)

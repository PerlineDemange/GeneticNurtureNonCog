# Compare adoptees and non-adoptees birth weight and birth place 
# Code by R. Cheesmand & P. Demange 
# K clustering based on https://uc-r.github.io/kmeans_clustering


# 1. Get birth weight and birth place data with Python #####

#ID field birth weight : 20022
#ID field birth place coordinate east: 130
#ID field birth place coordinate north: 129

module load pre2019
module load python/3.4.2

echo "Field,Description" > birth.in.csv
echo "20022,Birthweight" >> birth.in.csv
echo "130,BirthCoordEast" >> birth.in.csv
echo "129,BirthCoordNorth" >> birth.in.csv

python3 phenotypes/extract_ukbb_variables.py \
phenotypes/ukb30545.csv \
birth.in.csv \
birth.out.csv 

# 2. Clean birth weight and birth place data in R #####
module unload python
module load 2019
module load R/3.5.1-foss-2019b

R
# Set- up 
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(psych)
library(knitr)
library(data.table)
library(ggplot2)
# https://uc-r.github.io/kmeans_clustering


setwd("/home/pdemange/UKB/Adoption/")
birth_weight_place <- read.csv("../birth.out.csv",header=T,sep = ",")


birth_weight_place<-fread('birth_weight_place.txt',h=T)

colnames(birth_weight_place)<-c("ID",
                                'Place.of.birth.in.UK.north.co.ordinate.0.0',
                                'Place.of.birth.in.UK.north.co.ordinate.1.0',
                                'Place.of.birth.in.UK.north.co.ordinate.2.0',                                
                                'Place.of.birth.in.UK.east.co.ordinate.0.0',
                                'Place.of.birth.in.UK.east.co.ordinate.1.0',
                                'Place.of.birth.in.UK.east.co.ordinate.2.0',
                                'Birth.weight.0.0','Birth.weight.1.0','Birth.weight.2.0')


# remove -1 values for coords
birth_weight_place[birth_weight_place<0] <- NA

# combine instances --> 3 single variables.
birth_weight_place$Birth_weight <- rowMeans(birth_weight_place[,c('Birth.weight.0.0',
                                                                  'Birth.weight.1.0',
                                                                  'Birth.weight.2.0')],
                                            na.rm=TRUE)
birth_weight_place$Birth_place_east <- rowMeans(birth_weight_place[,c('Place.of.birth.in.UK.east.co.ordinate.0.0',
                                                                      'Place.of.birth.in.UK.east.co.ordinate.1.0',
                                                                      'Place.of.birth.in.UK.east.co.ordinate.2.0')], 
                                                na.rm=TRUE)
birth_weight_place$Birth_place_north <- rowMeans(birth_weight_place[,c('Place.of.birth.in.UK.north.co.ordinate.0.0',
                                                                       'Place.of.birth.in.UK.north.co.ordinate.1.0',
                                                                       'Place.of.birth.in.UK.north.co.ordinate.2.0')],
                                                 na.rm=TRUE)
birth_weight_place<-as.data.frame(birth_weight_place)
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
birth_weight_place[is.nan(birth_weight_place)] <- NA

# 3. Compare adoptees and non-adoptees  #####
## *3.1 Merge with adopted and non-adopted IDs.-----
adop<- fread("../EA/adoptees.csv", data.table=F, header=T) 
control <- fread("../EA/nonadopted_6.5k.csv", data.table=F, header=T) 

adop_birth<-merge(adop, birth_weight_place, by.x ='IID',by.y="ID")
control_birth<-merge(control, birth_weight_place, by.x ='IID',by.y="ID")

## *3.2 Compare birth weight-----

var.test(adop_birth$Birth_weight, control_birth$Birth_weight) #p=e-09
t.test(adop_birth$Birth_weight, control_birth$Birth_weight, var.equal = FALSE)
# Welch Two Sample t-test
# 
# data:  adop_birth$Birth_weight and control_birth$Birth_weight
# t = -10.316, df = 3126.6, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.2576322 -0.1753415
# sample estimates:
#   mean of x mean of y
# 3.118329  3.334815

## *3.3 Compare birth place-----
birth_place <- rbind(adop_birth, control_birth)
write.table(birth_place, file="adoptees_birth-weight_place.csv", row.names=F, quote=F) 

#rest of the analyses done on local machine as LISA is not functionning 
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(psych)
library(knitr)
library(data.table)
library(ggplot2)

setwd("C:/Users/Admin/Documents/GeneticNurtureNonCog/UKB/Adoptees") 
birth_place <- fread("adoptees_birth-weight_place.csv")

coords<-birth_place[,c(4,16:17)]
head(coords)
coords<-na.omit(coords)
coords$Adoptee<-ifelse(coords$Adopted== 1,'Yes',
                       ifelse(coords$Adopted==0,'No',NA))
coords$Birth_place_east_sc<-scale(coords$Birth_place_east)
coords$Birth_place_north_sc<-scale(coords$Birth_place_north)


ggplot(coords, aes(x=Birth_place_east_sc, y=Birth_place_north_sc, colour = Adoptee)) + 
  geom_point(alpha = 0.5 )+
  theme_classic()


# **3.3.1 K-means clustering --------
# cluster=kmeans(coords[,c('Birth_place_east_sc','Birth_place_north_sc')],2)
# plot(cluster$centers)

set.seed(123)

df_a<-subset(coords, Adoptee=='Yes', select=c('Birth_place_east_sc','Birth_place_north_sc'))
df_na<-subset(coords, Adoptee=='No', select=c('Birth_place_east_sc','Birth_place_north_sc'))

# define clusters such that the total intra-cluster variation 
# (known as total within-cluster variation or total within-cluster sum of square) is minimized
# For each k, calculate the total within-cluster sum of square (wss)
# Plot the curve of wss according to the number of clusters k.
# The location of a bend (knee) in the plot is generally considered as an indicator of the appropriate number of clusters.

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(df_a, k, nstart = 10 )$tot.withinss
}
# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15
# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares",
     main="Geographical clustering of adoptees' birth places")

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(df_na, k, nstart = 10 )$tot.withinss
}
# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15
# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)
plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares",
     main="Geographical clustering of non-adoptees' birth places")

# **3.3.2 Visualize clusters --------

# cluster everyone together
df<-subset(coords, select=c('Birth_place_east_sc','Birth_place_north_sc'))

final <- kmeans(df, 4, nstart = 25)
fviz_cluster(final, data = df)
print(final)

# adopted 
final_a <- kmeans(df_a, 4, nstart = 25)
print(final_a)

fviz_cluster(final_a, data = df_a)

# non-adopted
final_na <- kmeans(df_na, 4, nstart = 25)
print(final_na)
fviz_cluster(final_na, data = df_na)



# **3.3.3 Compare numbers in each region --------

aa<-as.data.frame(table(final_a$cluster))
naa<-as.data.frame(table(final_na$cluster))
naa$Region=c('north','midland','wales','south')
aa$Region=c('wales','midland','south','north')
aa$Var1=NULL
naa$Var1=NULL

colnames(aa)[1]<-"N_adopted"
colnames(naa)[1]<-"N_nonadopted"
ns<-merge(aa,naa, by = 'Region')

ns$N_adopted_percent <- ns$N_adopted / sum(ns$N_adopted)*100
ns$N_nonadopted_percent <- ns$N_nonadopted / sum(ns$N_nonadopted)*100
kable(ns)

# |Region  | N_adopted| N_nonadopted| N_adopted_percent| N_nonadopted_percent|
#   |:-------|---------:|------------:|-----------------:|--------------------:|
#   |midland |       825|         3463|         13.975944|            57.116939|
#   |north   |       517|          590|          8.758259|             9.731156|
#   |south   |      1537|          746|         26.037608|            12.304140|
#   |wales   |      3024|         1264|         51.228189|            20.847765|

# **3.3.4 Compare cluster means --------
na<-final_na$centers
rownames(na)<-c('north','midland','wales','south')
colnames(na)<-c('Nonadopted_Birth_place_east','Nonadopted_Birth_place_north')
na<-na[c(3:2,4,1),]
ad<-final_a$centers
rownames(ad)<-c('wales','midland','south','north')
colnames(ad)<-c('Adopted_Birth_place_east','Adopted_Birth_place_north')

all<-cbind(ad, na)
all<-all[,c(1,3,2,4)]
kable(round(all,3))


# |        | Adopted_Birth_place_east| Nonadopted_Birth_place_east| Adopted_Birth_place_north| Nonadopted_Birth_place_north|
#   |:-------|------------------------:|---------------------------:|-------------------------:|----------------------------:|
#   |wales   |                   -0.081|                       1.318|                     0.339|                       -1.007|
#   |midland |                   -0.947|                      -0.073|                    -1.098|                        0.384|
#   |south   |                    1.321|                      -0.998|                    -1.019|                       -1.050|
#   |north   |                   -1.509|                      -1.510|                     1.934|                        1.993|
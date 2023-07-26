# this script will introduce you to occupancy modeling using package 'unmarked'
# The objectives of this exercise are:
# 1) to introduce you to simulating occupancy data
# 2) to set up occupancy data for package 'unmarked'
# 3) to run simple (single-season) occupancy models and evaluate outputs

install.packages(c("unmarked", "MuMIn"))

library(unmarked)
library(MuMIn)


# OCCUPANCY SIMULATIONS ---------------------------------------------------

# 1) let's simulate occupancy data
# Objective: understand the structure of single-season occupancy data, and evaluate patterns in simulated data 

## First, set up empty vector of 'true' occupancy (STATE PROCESS)
n_sites=100
n_vis = 10

tocc <- rep(0, n_sites)
tocc <- as.data.frame(tocc)
names(tocc) <- c("T1")
tocc
 
### set up an empty observed occupancy matrix (OBSERVATION PROCESS)
obsocc <- matrix(rep(0, n_sites*n_vis), ncol=(n_vis),
                 dimnames=list(paste("site", 1:n_sites, sep=""), 
                               paste("v", 1:n_vis, sep="")))
obsocc <- as.data.frame(obsocc)
obsocc

## Now let's play with various values for psi (occupancy) and p (detection)
psi = 0.7
p = 0.8

# fill the observed occupancy matrix
for(i in 1:n_sites) {
  tocc[i,1] = rbinom(1, 1, psi)                      # Simulate site occupancy at t1
  obsocc[i,1:n_vis] = rbinom(n_vis, 1, tocc[i,1]*p)   # Simulate observed P/A at t1 using detection prob. p
}

# examine the simulated probability of occupancy
sim_psi = sum(tocc)/n_sites
sim_psi

# examine the true and observed occupancy 
tocc
sum(tocc)

obsocc
sum(obsocc)

# calculate the sum of successful detections per site 
obsocc$sums <- rowSums(obsocc)
obsocc

# calculate the number of sites where the species was actually observed
# the difference between the true occupancy and number of sites where species
# was detected gives us the observation error

# which sites have 0 detections are saved in a vector, 
# then calculate the length of that particular vector
# subtract this value from the number of sites!!!
n_sites - length(which(obsocc$sums == 0))


Group <- c(rep("G",50), rep("B",50))

dataUMF <-  unmarkedFrameOccu(
  y = obsocc,
  siteCovs = data.frame(Group)
  #obsCovs = list()
)

plot(dataUMF)
str(dataUMF)

# run an occupancy model with 2 submodels: 
# 1st ~ is for detection, 2nd ~ is for occupancy

model.Null = occu(~1 ~1, data=dataUMF)      
summary(model.Null)

model.1 = occu(~1 ~Group, data=dataUMF)
summary(model.1)

oc = predict(model.1, type="state")
oc

det = predict(model.1, type="det")
det

### Let's simulate different occupancy rates at two groups of sites and evaluate whether the model is able to 
# identify the differences. We will simulate another 100 sites with psi = 0.3 and p=0.3 

tocc1 <- rep(0, n_sites)
tocc1 <- as.data.frame(tocc1)
names(tocc1) <- c("T1")
tocc1

### set up an empty observed occupancy matrix (OBSERVATION PROCESS)
obsocc1 <- matrix(rep(0, n_sites*n_vis), ncol=(n_vis),
                 dimnames=list(paste("site", 1:n_sites, sep=""), paste("v", 1:n_vis, sep="")))
obsocc1 <- as.data.frame(obsocc1)
obsocc1

## Now let's play with various values for psi (occupancy) and p (detection)
psi1 = 0.3
p1 = 0.3

# fill the observed occupancy matrix
for(i in 1:n_sites) {
  tocc1[i,1] = rbinom(1, 1, psi1)                      # Simulate site occupancy at t1
  obsocc1[i,1:n_vis] = rbinom(n_vis, 1, tocc1[i,1]*p1)   # Simulate observed P/A at t1 using detection prob. p
}

tocc1
sum(tocc1)

obsocc1
sum(obsocc1)

# put together the 2 datasets: tocc and tocc1, obsocc and obsocc1

tocc2 = rbind(tocc, tocc1)
tocc2

obsocc2 = rbind(obsocc, obsocc1)
obsocc2

# create new object for unmarked and run models
Group2 <- c(rep("G",100), rep("B",100))

dataUMF2 <-  unmarkedFrameOccu(
  y = obsocc2,
  siteCovs = data.frame(Group2),
  #obsCovs = list()
)

plot(dataUMF2)
str(dataUMF2)

# run an occupancy model with 2 submodels: 
# 1st ~ is for detection, 2nd ~ is for occupancy

model.Null = occu(~1 ~1, data=dataUMF2)      
summary(model.Null)
plogis(0.0802)
plogis(0.494)

model.1 = occu(~1 ~Group2-1, data=dataUMF2)
summary(model.1)
plogis(-0.663)
plogis(0.848)

model.2 = occu(~Group2-1 ~Group2-1, data=dataUMF2)
summary(model.2)
plogis(-0.624)
plogis(0.847)


# create model selection table
modlist <- list(Null = model.Null, psiGroup = model.1, psiGroup.pGroup = model.2)

# create model selection table
selection = model.sel(modlist)
selection

psi.pred = predict(model.2, type="state")
psi.pred

plogis(coef(model.2))



###### ANALYSIS OF OCCUPANCY DATASET ########
# let's analyze an occupancy dataset of mink frog (Lithobates septentrionalis) 
# occurrence from the Adirondack mountains

# This document describes the occupancy analysis of Mink Frogs (Lithobates septentrionalis) 
# in the Adirondack Mountains, New York State

# The specific objectives of this document are:
#   To document the data management and occupancy analysis process;
#   To conduct a simple occupancy analysis in R using the unmarked package; and
#   To explore occupancy predictions and prepare the lay the fondation for multi-year (dynamic) occupancy models

setwd('C:/Users/viorel/Dropbox/OhioU/teaching/PopEco_Fall2021/occupancy')

mink <- read.csv("minkfrogs_occ.csv", header=T)

names(mink)
str(mink)

attach(mink) # we need to 'attach' the file for unmarked to work; usually it is not recommended to use this function

library(unmarked)
library(MuMIn)

# calculate the naive occupancy (not adjusted for detection probability)
# divide the number where species was detected to thenumber of sites
minkfrog_naive <- sum(apply(mink[,3:6], 1, sum) > 0) / nrow(mink)
minkfrog_naive

# Set up an R object that can be read by the 'unmarked' package
# the object contains 3 types of data:
#   1) detection/nondetection data (the observations grouped into a matrix "y")
#   2) site-level covariates; these can be any variables that do not vary with survey occasion,
#             such as landscape characteristics, habitat types, site characteristics etc.
#             - these variables are used to model OCCUPANCY; they can also be used to model DETECTION (not shown here)
#   3) observation-level coveriates; these can be any variables recorded at the time of survey
#             that may influence detection
#             - these variables are used to model DETECTION

minkUMF <- with(mink, {
  unmarkedFrameOccu(
    
    y = cbind(LISE1, LISE2, LISE3, LISE4),
    
    siteCovs = data.frame(JulyTemp, DO, DistMFbreed, Beaver),
    
    obsCovs = list(Date = cbind(Julian1, Julian2, Julian3, Julian4),
                   Wind = cbind(Wind1, Wind2, Wind3, Wind4),
                   Sky = cbind(Sky1, Sky2, Sky3, Sky4)))
})

plot(minkUMF)
str(minkUMF)

# first we need to evaluate the best predictors for detection
# therefore, we run a model that contains some site-related variables (for modeling occuopancy) 

mink.d1 <- occu(~Date + Wind + Sky ~JulyTemp + DO + DistMFbreed, data=minkUMF)
#summary(mink.d1)

mink.d2 = occu(~Date + Sky ~JulyTemp + DO + DistMFbreed, data=minkUMF)
summary(mink.d2)

mink.d3 = occu(~Date  ~JulyTemp + DO + DistMFbreed, data=minkUMF)
#summary(mink.d3)

mink.d4 = occu(~Sky  ~JulyTemp + DO + DistMFbreed, data=minkUMF)
#summary(mink.d4)

mink.d5 = occu(~Date + Wind  ~JulyTemp + DO + DistMFbreed, data=minkUMF)
#summary(mink.d5)

# create a list for the models and given them names 
modlist_d <- list(DateWindSky = mink.d1, DateSky = mink.d2, Date = mink.d3, 
                   Sky = mink.d4, DateWind = mink.d5)

# create model selection table
selection_d = model.sel(modlist_d)
selection_d


# BEST models is mink.d1: variables Date and Sky and Wind are best for detection


# we now use these 3 variables to build several models that explore occupancy predictors

# climate predictors 
mink.o1 = occu(~Date + Sky + Wind ~JulyTemp, data=minkUMF)

# landscape predictors 
mink.o2 = occu(~Date + Sky + Wind ~DistMFbreed + factor(Beaver), data=minkUMF)

# pond predictors 
mink.o3 = occu(~Date + Sky + Wind ~DO, data=minkUMF)
summary(mink.o3)

# climate and pond predictors
mink.o4 = occu(~Date + Sky + Wind ~JulyTemp + DO, data=minkUMF)
summary(mink.o4)

# landscape and pond predictors
mink.o5 = occu(~Date + Sky + Wind ~DistMFbreed + factor(Beaver) + DO, data=minkUMF)

# landscape and climate predictors
mink.o6 = occu(~Date + Sky + Wind ~DistMFbreed + factor(Beaver) + JulyTemp, data=minkUMF)

modlist_o <- list(JTemp = mink.o1, MFBreedBeaver = mink.o2, DO = mink.o3, JTempDO = mink.o4,
                  MFBreedBeaverDO = mink.o5, FBreedBeaverJTemp = mink.o6)

# create model selection table
selection_o = model.sel(modlist_o)
selection_o

# model averaging: instead of choosing the BEST model, 
# draw inferences on occupancy and detection
# based on all models, weighed by their ACIc weight
# use all models whose cumulative AICc weight <=0.95
# This function returns the averaged model coefficients 
# (RECENT WORK CLAIMS THIS IS WRONG, SO WE WILL NOT DWELL ON THEM) 
mod.avg = model.avg(selection_o, cumsum(weight) <= .95)
mod.avg


# detection probability estimates by occasions using model averaging (THIS IS CORRECT)
detpred = predict(mod.avg, type="det")
detpred
hist(detpred$fit)
mean(detpred$fit)

# occupancy predictions by site using moedl averaging (THIS IS CORRECT)
psipred = predict(mod.avg, type="state")
psipred
hist(psipred$fit)

# get the probability of occupancy CORRECTED for imperfect detection
# also known as Proportion of Areas Occupied (PAO) = mean of 'psipred'
PAO = mean(psipred$fit)
PAO

# difference between naive and occupancy predictions
PAO-minkfrog_naive

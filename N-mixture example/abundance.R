
### This code analyzes data on longicorn beetle abundance collected
### in the Iron Gates Natural Park, Romania, and published in 
### Brodie, B. S., V. D. Popescu, R. Iosif, et al. 2019. Non-lethal monitoring of 
### longicorn beetle communities using generic pheromone lures and occupancy models. 
### Ecological Indicators 101:330-340.

setwd('...')

beetles <- read.csv("cerambycidae.csv")
names(beetles)
attach(beetles)
str(beetles)

# for Ruben to understand, run abundance models for species: Xarv, Pdet, Ptest, Mneb, Aclav, Cmut

library(unmarked)
#library(emdbook)
library(MuMIn)
options(scipen=999)


### setup 'unmarked' dataframe using site-specific and survey-specific covariates ###
# survey-specific covariates are: JulianDay (day of year) and Precipitation
# site-specific covariates are Site, Canopy cover, Ground cover, Tree density,
# Number of tree species, Dominant tree species, Average tree diameter,
# Number of stumps, and Volume of downed wood

# bind the survey-specific covariates and standardize them
JulianDay = matrix(cbind(Julian1, Julian2, Julian3, Julian4,
                         Julian5, Julian6, Julian7, Julian8,
                         Julian9, Julian10, Julian11, Julian12))
JulianDay = matrix(scale(JulianDay), nrow=50, byrow=F)

Precipitation = matrix(cbind(Precip1, Precip2, Precip3, Precip4,
                             Precip5, Precip6, Precip7, Precip8,
                             Precip9, Precip10, Precip11, Precip12))
Precipitation = matrix(scale(Precipitation), nrow=50, byrow=F)

# create an unmarked objects for theone of the focal species, Plagionotus detritus

PdetUMF <- with(beetles, {
  unmarkedFramePCount(
    y = cbind(Pdet1, Pdet2, Pdet3, Pdet4, Pdet5, Pdet6, Pdet7, Pdet8, Pdet9, Pdet10, Pdet11, Pdet12),
    siteCovs = data.frame(Site, Lures, scale(Canopy), scale(GroundCover), scale(TreeDens), 
                          scale(NoTreeSp), Dominant.1, scale(AvgDiameter), scale(NoStumps), scale(VolDeadWood1)),
    obsCovs = list(Julian = JulianDay,
                   Precip = Precipitation))
  
})

plot(PdetUMF)


### BULD AND RUN DETECTION MODELS (using a full model for occupancy)
# run several models for detection, as we did for the Mink Frog case study


m0 <- pcount(~1 ~Site+scale.Canopy.+scale.TreeDens.+scale.AvgDiameter.+
               scale.NoStumps.+scale.VolDeadWood1., mixture="P", K=200, data=PdetUMF)

m1 <- pcount(~Julian ~Site+scale.Canopy.+scale.TreeDens.+scale.AvgDiameter.+
               scale.NoStumps.+scale.VolDeadWood1., mixture="P", K=200, data=PdetUMF)

m2 <- pcount(~Precip ~Site+scale.Canopy.+scale.TreeDens.+scale.AvgDiameter.+
               scale.NoStumps.+scale.VolDeadWood1., mixture="P", K=200, data=PdetUMF)

m3 <- pcount(~Julian+Precip ~Site+scale.Canopy.+scale.TreeDens.+scale.AvgDiameter.+
               scale.NoStumps.+scale.VolDeadWood1., mixture="P", K=200, data=PdetUMF)


### TO DO IN CLASS: BUILD 3 OTHER MODELS FOR DETECTION ###

modelList <- fitList(Null = m0, "mod1" = m1, "mod2" = m2, "mod3" = m3)
modSel_Pdet <- modSel(modelList, nullmod = 'Null')
modSel_Pdet



### BUILD AND RUN MODELS FOR ABUNDANCE ###

Pdet_Null <- pcount(~Julian+Precip  ~1,
                  mixture="P", K=200, data=PdetUMF)

Pdet_m0 <- pcount(~Julian+Precip  ~scale(NoStumps)+scale(VolDeadWood1),
                  mixture="P", K=200, data=PdetUMF) 

Pdet_m1 <- pcount(~Julian+Precip  ~scale(GroundCover)+scale(Canopy),
                  mixture="P", K=200, data=PdetUMF)
summary(Pdet_m1)

### TO DO IN CLASS: BUILD SEVERAL OTHER MODELS FOR MODLEING ABUNDANCE ###

# rank models by AIC for PDET
modelList_Pdet <- fitList(Null = Pdet_Null, m0 = Pdet_m0, m1 = Pdet_m1)

modSel_Pdet <- modSel(modelList_Pdet, nullmod = 'Null')
modSel_Pdet


# extract coeficients for best model
summary(Pdet_m1)
coef(Pdet_m1)

modelaveragedmodels <- model.avg(modSel_Pdet)

# estimate model averaged individuals per survey #
Pdet_abund = predict(Pdet_m1, type="state")

#avg abu for Eselnita
Pdet_abund_Eselnita <- Pdet_abund[1:25,]
avg_Eselnita <- mean(Pdet_abund_Eselnita[,1])
se_Eselnita <- sqrt(sum((1/25^2)*Pdet_abund_Eselnita[,2]^2))
cv_Eselnita <- se_Eselnita/avg_Eselnita
ci.lo_Eselnita <- avg_Eselnita - se_Eselnita * qnorm(0.975) # 95% CI
ci.hi_Eselnita <- avg_Eselnita + se_Eselnita * qnorm(0.975)
avg_Eselnita
se_Eselnita
cv_Eselnita
ci.lo_Eselnita
ci.hi_Eselnita

#avg abu for Mala
Pdet_abund_Mala <- Pdet_abund[26:50,]
avg_Mala <- mean(Pdet_abund_Mala[,1])
se_Mala <- sqrt(sum((1/25^2)*Pdet_abund_Mala[,2]^2))
cv_Mala <- se_Mala/avg_Mala
ci.lo_Mala <- avg_Mala - se_Mala * qnorm(0.975) # 95% CI
ci.hi_Mala <- avg_Mala + se_Mala * qnorm(0.975)
avg_Mala
se_Mala
cv_Mala
ci.lo_Mala
ci.hi_Mala

# estimate model averaged probability of detection #
Pdet_det = predict(modelList_Pdet, type="det")

#avg det
avg <- mean(Pdet_det[,1])
se <- sqrt(sum((1/(nrow(Pdet_det))^2)*Pdet_det[,2]^2))
cv <- se/avg
ci.lo <- avg - se * qnorm(0.975) # 95% CI
ci.hi <- avg + se * qnorm(0.975)
avg
se
cv
ci.lo
ci.hi

# END MODELS

library(unmarked)
library(MuMIn)

options(digits = 4) 

# can change to your Working Directory (computer path)
setwd("~/Desktop/Popescu Lab/THESIS/CLEANED EXCEL")
setwd("C:\\Users\\Viorel\\Dropbox\\OhioU\\projects\\Amphibians\\WAYNE occupancy\\Jan 16th")

sites <- read.csv("Nov 15 Site Data.csv", header=T)
str(sites)
summary(sites)
frogs <- read.csv("Nov 8.csv", header=T)
str(frogs)
frogs$Spp <- as.factor(frogs$Spp)
str(frogs)
summary(frogs)
levels(frogs$Spp)

AT <- subset(frogs, Spp == "American Toad")
str(AT)

# extract observations; values >1 will be automatically converted to 1's when we crete the 
# object for unmarked
yMat <- AT[,4:59]
yMat

# extract observation covariates (used for detection)
det.covs.TEMP <- AT[,172:227]
str(det.covs.TEMP)
det.covs.TEMP <- as.numeric(det.covs.TEMP)

det.covs.SKY <- AT[,228:283]
str(det.covs.SKY)
summary(det.covs.SKY)
levels(det.covs.SKY)

det.covs.TIME <- AT[,60:115]
str(det.covs.TIME)
summary(det.covs.TIME)
which(is.na(det.covs.TIME))
#need to change?

det.covs.MONTH <- AT[,352:407]
str(det.covs.MONTH)


# read in the site data
site.covariates <- read.csv("Site Covariates.csv", header=T)
str(site.covariates)

# scale site covariates (land cover etc.)
s.100.dev <- scale(site.covariates$P100V21V22)
s.100.devhigh <- scale(site.covariates$P100V23V24)
s.100.for <- scale(site.covariates$P100V31V41V42V43)
s.100.opn <- scale(site.covariates$P100V52V71)
s.100.agr <- scale(site.covariates$P100V81V82)
s.100.wet <- scale(site.covariates$P100V90V95)
s.250.dev <- scale(site.covariates$P250V21V22)
s.250.devhigh <- scale(site.covariates$P250V23V24)
s.250.for <- scale(site.covariates$P250V31V41V42V43)
s.250.opn <- scale(site.covariates$P250V52V71)
s.250.agr <- scale(site.covariates$P250V81V82)
s.250.wet <- scale(site.covariates$P250V90V95)
s.500.dev <- scale(site.covariates$P500V21V22)
s.500.devhigh <- scale(site.covariates$P500V23V24)
s.500.for <- scale(site.covariates$P500V31V41V42V43)
s.500.opn <- scale(site.covariates$P500V52V71)
s.500.agr <- scale(site.covariates$P500V81V82)
s.500.wet <- scale(site.covariates$P500V90V95)
s.beaver <- scale(site.covariates$Beaver)
s.dw <- scale(site.covariates$DistNearWetland)
s.dr <- scale(site.covariates$DistNearRoad)
s.dd <- scale(site.covariates$DistNearDevelop)
s.num <- scale(site.covariates$WSNumber)


# site covariates data.frame
site.covs <- data.frame(s.100.dev, s.100.devhigh, s.100.for, s.100.agr, s.100.opn, s.100.wet, s.250.dev, s.250.devhigh, s.250.for, s.250.agr, s.250.opn, s.250.wet, s.500.dev, s.500.devhigh, s.500.for, s.500.agr, s.500.opn, s.500.wet, s.dw, s.beaver, s.dr, s.dw, s.num)


# prepare the unmarked multi-year colext data.frame
nyear <- 14 #14 years of data
nsites <- 30 # 30 sites

# create yearly site covariates for modeling extinction/colonization
years <- matrix(rep(seq(1:(nyear-1)), nsites), ncol=nyear-1, byrow=T,
                dimnames=list(paste("site", 1:nsites, sep=""), 
                              paste("Yr", rep(1:(nyear-1)), sep="")))
# we need to add a column that will never be used 
yrs.Junk <- matrix(rep("Junk", nsites*1), ncol=1)
years <- cbind(years,yrs.Junk)
colnames(years) <- paste("Yr", rep(1:nyear), sep="")
years <- as.data.frame(years)
years <- data.frame(lapply(years, as.factor))
years

# prepare the matrix for yearly covariates for modeling colonization and extinction
# we have to transpose this matrix to be read in unmarked
yearly_covs = data.frame(YEARS = matrix(t(years)))
yearly_covs


# construct the unmarked object
oUMF <- unmarkedMultFrame(
  y = yMat,                
  siteCovs = site.covs,   # site covariates for modeling occupancy in Year1
  obsCovs = list(SKY = det.covs.SKY, TIME = det.covs.TIME, MONTH = det.covs.MONTH, TEMP = det.covs.TEMP), 
    #TEMP = , TIME = , MONTH = ), # Temp, Time of day, Month, Sky
  yearlySiteCovs = yearly_covs,   # yearly site covariates for modeling extinction and colonization
  numPrimary=nyear)

summary(oUMF)
plot(oUMF)

# model estimates psi (occupancy), gamma (colonization), epsilon (extinction), p (detection)
Null.fit <- colext(psiformula = ~1, gammaformula = ~1, 
                   epsilonformula = ~1, pformula = ~1, oUMF)
summary(Null.fit)
#!!sky covs change to 1-5 and check spaces again 




#step 1 is identify best detection model 
nd <- data.frame(YEARS=c('1','2','3','4','5','6','7','8',
                         '9','10','11','12','13'))

model.1 <- colext(psiformula = ~s.100.for + s.dw + s.dr+ s.100.dev +s.100.wet,
                  gammaformula = ~1, 
                  epsilonformula = ~1,
                  pformula = ~as.factor(SKY) + TEMP + MONTH + TIME, oUMF)
summary(model.1)
predict(model.1, type = "ext", newdata = nd)

model.2 <- colext(psiformula = ~s.100.for + s.dw + s.dr+ s.100.dev +s.100.wet,
                  gammaformula = ~1, 
                  epsilonformula = ~1,
                  pformula = ~as.factor(SKY) + TEMP + MONTH, oUMF)
summary(model.2)
predict(model.2, type = "ext", newdata = nd)

model.3 <- colext(psiformula = ~s.100.for + s.dw + s.dr+ s.100.dev +s.100.wet,
                  gammaformula = ~1, 
                  epsilonformula = ~1,
                  pformula = ~as.factor(SKY) + TIME + MONTH, oUMF)
summary(model.3)
predict(model.3, type = "ext", newdata = nd)

model.4 <- colext(psiformula = ~s.100.for + s.dw + s.dr+ s.100.dev +s.100.wet,
                  gammaformula = ~1, 
                  epsilonformula = ~1,
                  pformula = ~as.factor(SKY) + TEMP, oUMF)
summary(model.4)
predict(model.4, type = "ext", newdata = nd)

model.5 <- colext(psiformula = ~s.100.for + s.dw + s.dr+ s.100.dev +s.100.wet,
                  gammaformula = ~1, 
                  epsilonformula = ~1,
                  pformula = ~as.factor(SKY) + MONTH, oUMF)
summary(model.5)
predict(model.5, type = "ext", newdata = nd)

model.6 <- colext(psiformula = ~s.100.for + s.dw + s.dr+ s.100.dev +s.100.wet,
                  gammaformula = ~1, 
                  epsilonformula = ~1,
                  pformula = ~as.factor(SKY) + TIME, oUMF)
summary(model.6)
predict(model.6, type = "ext", newdata = nd)

#if error is it tilde in fromt of TEMP
model.7 <- colext(psiformula = ~s.100.for + s.dw + s.dr+ s.100.dev +s.100.wet,
                  gammaformula = ~1, 
                  epsilonformula = ~1,
                  pformula = ~TEMP + MONTH + TIME, oUMF)
summary(model.7)
predict(model.7, type = "ext", newdata = nd)

model.8 <- colext(psiformula = ~s.100.for + s.dw + s.dr+ s.100.dev +s.100.wet,
                  gammaformula = ~1, 
                  epsilonformula = ~1,
                  pformula = ~TEMP + MONTH, oUMF)
summary(model.8)
predict(model.8, type = "ext", newdata = nd)
predict(model.8, type = "col", newdata = nd)

model.9 <- colext(psiformula = ~s.100.for + s.dw + s.dr+ s.100.dev +s.100.wet,
                  gammaformula = ~1, 
                  epsilonformula = ~1,
                  pformula = ~TEMP + TIME, oUMF)
summary(model.9)
predict(model.9, type = "ext", newdata = nd)

model.10 <- colext(psiformula = ~s.100.for + s.dw + s.dr+ s.100.dev +s.100.wet,
                  gammaformula = ~1, 
                  epsilonformula = ~1,
                  pformula = ~MONTH + TIME, oUMF)
summary(model.10)
predict(model.10, type = "ext", newdata = nd)

model.11 <- colext(psiformula = ~s.100.for + s.dw + s.dr+ s.100.dev +s.100.wet,
                  gammaformula = ~1, 
                  epsilonformula = ~1,
                  pformula = ~as.factor(SKY) + TEMP + TIME, oUMF)
summary(model.11)
predict(model.11, type = "ext", newdata = nd)


# Model selection

my_models <- list(m1 = model.1, m2 = model.2, m3 = model.3, m4 = model.4, m5 = model.5, m6 = model.6, m7 = model.7, m8 = model.8, m9 = model.9, m10 = model.10, m11 = model.11)

modsel <- model.sel(my_models)
modsel

install.packages("Hmisc")
library(Hmisc)

#run once to find coorelations and determine what cannot run
rcorr(as.matrix(site.covs), type="spearman")

#AT Models 
#epsilon ~1 for 15 models below
model.LC100 <- colext(psiformula = ~s.100.dev + s.100.for + s.100.agr + s.100.opn + s.100.devhigh,
                   gammaformula = ~1, 
                   epsilonformula = ~1,
                   pformula = ~TEMP + MONTH, oUMF)
summary(model.LC100)
predict(model.LC100, type = "ext", newdata = nd)

model.LC250 <- colext(psiformula = ~s.250.dev + s.250.for + s.250.opn + s.250.devhigh,
                      gammaformula = ~1, 
                      epsilonformula = ~1,
                      pformula = ~TEMP + MONTH, oUMF)
summary(model.LC250)
predict(model.LC250, type = "ext", newdata = nd)

model.LC500 <- colext(psiformula = ~s.500.dev + s.500.opn + s.500.devhigh,
                      gammaformula = ~1, 
                      epsilonformula = ~1,
                      pformula = ~TEMP + MONTH, oUMF)
summary(model.LC500)
predict(model.LC500, type = "ext", newdata = nd)

model.LS1 <- colext(psiformula = ~s.dw + s.dr,
                      gammaformula = ~1, 
                      epsilonformula = ~1,
                      pformula = ~TEMP + MONTH, oUMF)
summary(model.LS1)
predict(model.LS1, type = "ext", newdata = nd)

model.LS1 <- colext(psiformula = ~s.dw + s.dr,
                    gammaformula = ~1, 
                    epsilonformula = ~1,
                    pformula = ~TEMP + MONTH, oUMF)
summary(model.LS1)
predict(model.LS1, type = "ext", newdata = nd)

model.LS2 <- colext(psiformula = ~s.dw,
                    gammaformula = ~1, 
                    epsilonformula = ~1,
                    pformula = ~TEMP + MONTH, oUMF)
summary(model.LS2)
predict(model.LS2, type = "ext", newdata = nd)

model.LS3 <- colext(psiformula = ~s.dw + s.dd,
                    gammaformula = ~1, 
                    epsilonformula = ~1,
                    pformula = ~TEMP + MONTH, oUMF)
summary(model.LS3)
predict(model.LS3, type = "ext", newdata = nd)

model.LS4 <- colext(psiformula = ~s.dw + s.dr + s.dd,
                    gammaformula = ~1, 
                    epsilonformula = ~1,
                    pformula = ~TEMP + MONTH, oUMF)
summary(model.LS4)
predict(model.LS4, type = "ext", newdata = nd)

#LSModel selection lines 280-283 ran already

#AT.LS <- list(m.1 = model.LS1, m.2 = model.LS2, m.3 = model.LS3, m.4 = model.LS4)

#modsel <- model.sel(AT.LS)
#modsel

#LS 2 and LS1 are best models of 4 LS ones, below combine with LC for 6 more models

model.B101 <- colext(psiformula = ~s.100.dev + s.100.for + s.100.agr + s.100.opn + s.100.devhigh + s.dr + s.dw,
                      gammaformula = ~1, 
                      epsilonformula = ~1,
                      pformula = ~TEMP + MONTH, oUMF)
summary(model.B101)
predict(model.B101, type = "ext", newdata = nd)

model.B251 <- colext(psiformula = ~s.250.dev + s.250.for + s.250.opn + s.250.devhigh + s.dr + s.dw,
                      gammaformula = ~1, 
                      epsilonformula = ~1,
                      pformula = ~TEMP + MONTH, oUMF)
summary(model.B251)
predict(model.B251, type = "ext", newdata = nd)

model.B501 <- colext(psiformula = ~s.500.dev + s.500.opn + s.500.devhigh + s.dr + s.dw,
                      gammaformula = ~1, 
                      epsilonformula = ~1,
                      pformula = ~TEMP + MONTH, oUMF)
summary(model.B501)
predict(model.B501, type = "ext", newdata = nd)

model.B102 <- colext(psiformula = ~s.100.dev + s.100.for + s.100.agr + s.100.opn + s.100.devhigh + s.dw,
                     gammaformula = ~1, 
                     epsilonformula = ~1,
                     pformula = ~TEMP + MONTH, oUMF)
summary(model.B102)
predict(model.B102, type = "ext", newdata = nd)

model.B252 <- colext(psiformula = ~s.250.dev + s.250.for + s.250.opn + s.250.devhigh + s.dw,
                     gammaformula = ~1, 
                     epsilonformula = ~1,
                     pformula = ~TEMP + MONTH, oUMF)
summary(model.B252)
predict(model.B252, type = "ext", newdata = nd)

model.B502 <- colext(psiformula = ~s.500.dev + s.500.opn + s.500.devhigh + s.dw,
                     gammaformula = ~1, 
                     epsilonformula = ~1,
                     pformula = ~TEMP + MONTH, oUMF)
summary(model.B502)
predict(model.B502, type = "ext", newdata = nd)

#epsilon ~year for 15 models below
model.LC100.y <- colext(psiformula = ~s.100.dev + s.100.for + s.100.agr + s.100.opn + s.100.devhigh,
                      gammaformula = ~1, 
                      epsilonformula = ~YEARS,
                      pformula = ~TEMP + MONTH, oUMF)
summary(model.LC100.y)
predict(model.LC100.y, type = "ext", newdata = nd)

model.LC250.y <- colext(psiformula = ~s.250.dev + s.250.for + s.250.opn + s.250.devhigh,
                      gammaformula = ~1, 
                      epsilonformula = ~YEARS,
                      pformula = ~TEMP + MONTH, oUMF)
summary(model.LC250.y)
predict(model.LC250.y, type = "ext", newdata = nd)

model.LC500.y <- colext(psiformula = ~s.500.dev + s.500.opn + s.500.devhigh,
                      gammaformula = ~1, 
                      epsilonformula = ~YEARS,
                      pformula = ~TEMP + MONTH, oUMF)
summary(model.LC500.y)
predict(model.LC500.y, type = "ext", newdata = nd)

model.LS1.y <- colext(psiformula = ~s.dw + s.dr,
                    gammaformula = ~1, 
                    epsilonformula = ~YEARS,
                    pformula = ~TEMP + MONTH, oUMF)
summary(model.LS1.y)
predict(model.LS1.y, type = "ext", newdata = nd)

model.LS1.y <- colext(psiformula = ~s.dw + s.dr,
                    gammaformula = ~1, 
                    epsilonformula = ~YEARS,
                    pformula = ~TEMP + MONTH, oUMF)
summary(model.LS1.y)
predict(model.LS1.y, type = "ext", newdata = nd)

model.LS2.y <- colext(psiformula = ~s.dw,
                    gammaformula = ~1, 
                    epsilonformula = ~YEARS,
                    pformula = ~TEMP + MONTH, oUMF)
summary(model.LS2.y)
predict(model.LS2.y, type = "ext", newdata = nd)

model.LS3.y <- colext(psiformula = ~s.dw + s.dd,
                    gammaformula = ~1, 
                    epsilonformula = ~YEARS,
                    pformula = ~TEMP + MONTH, oUMF)
summary(model.LS3.y)
predict(model.LS3.y, type = "ext", newdata = nd)

model.LS4.y <- colext(psiformula = ~s.dw + s.dr + s.dd,
                    gammaformula = ~1, 
                    epsilonformula = ~YEARS,
                    pformula = ~TEMP + MONTH, oUMF)
summary(model.LS4.y)
predict(model.LS4.y, type = "ext", newdata = nd)

model.B101.y <- colext(psiformula = ~s.100.dev + s.100.for + s.100.agr + s.100.opn + s.100.devhigh + s.dr + s.dw,
                     gammaformula = ~1, 
                     epsilonformula = ~YEARS,
                     pformula = ~TEMP + MONTH, oUMF)
summary(model.B101.y)
predict(model.B101.y, type = "ext", newdata = nd)

model.B251.y <- colext(psiformula = ~s.250.dev + s.250.for + s.250.opn + s.250.devhigh + s.dr + s.dw,
                     gammaformula = ~1, 
                     epsilonformula = ~YEARS,
                     pformula = ~TEMP + MONTH, oUMF)
summary(model.B251.y)
predict(model.B251.y, type = "ext", newdata = nd)

model.B501.y <- colext(psiformula = ~s.500.dev + s.500.opn + s.500.devhigh + s.dr + s.dw,
                     gammaformula = ~1, 
                     epsilonformula = ~YEARS,
                     pformula = ~TEMP + MONTH, oUMF)
summary(model.B501.y)
predict(model.B501.y, type = "ext", newdata = nd)

model.B102.y <- colext(psiformula = ~s.100.dev + s.100.for + s.100.agr + s.100.opn + s.100.devhigh + s.dw,
                     gammaformula = ~1, 
                     epsilonformula = ~YEARS,
                     pformula = ~TEMP + MONTH, oUMF)
summary(model.B102.y)
predict(model.B102.y, type = "ext", newdata = nd)

model.B252.y <- colext(psiformula = ~s.250.dev + s.250.for + s.250.opn + s.250.devhigh + s.dw,
                     gammaformula = ~1, 
                     epsilonformula = ~YEARS,
                     pformula = ~TEMP + MONTH, oUMF)
summary(model.B252.y)
predict(model.B252.y, type = "ext", newdata = nd)

model.B502.y <- colext(psiformula = ~s.500.dev + s.500.opn + s.500.devhigh + s.dw,
                     gammaformula = ~1, 
                     epsilonformula = ~YEARS,
                     pformula = ~TEMP + MONTH, oUMF)
summary(model.B502.y)
predict(model.B502.y, type = "ext", newdata = nd)

#AT Model Selection
model.list <- list(mLC1 = model.LC100, mLC2 = model.LC250, mLC3 = model.LC500, mLS1 = model.LS1, 
                   mLS2 = model.LS2, mLS3 = model.LS3, mLS4 = model.LS4, mB101 = model.B101, 
                   mB251 = model.B251, mB501 = model.B501, mB102 = model.B102, mB252 = model.B252, 
                   mB502 = model.B502, mLC1y = model.LC100.y, mLC2y = model.LC250.y, 
                   mLC3y = model.LC500.y, mLS1y = model.LS1.y, mLS2y = model.LS2.y, mLS3y = model.LS3.y, 
                   mLS4y = model.LS4.y, mB101y = model.B101.y, mB251y = model.B251.y, 
                   mB501y = model.B501.y, mB102y = model.B102.y, mB252y = model.B252.y, 
                   mB502y = model.B502.y)

modsel <- model.sel(model.list)
modsel

###############################
### New code from Viorel 17 Jan
###############################

# Assess goodness-of-fit
parboot(model.LS2)
plot(model.LS2)

# predicted occupancy in year 1
E.psi1 <- predict(model.LS2, type = "psi")
E.psi1
mean_psi1 <- mean(E.psi1[,1])
mean_psi1

# predict occupancy probabilities for years 2 through 14
E.psi <- projected(model.LS2)
E.psi

#E.psi.smoothed <- smoothed(model.LS2)
#E.psi.smoothed

## Find bootstrap standard errors for smoothed trajectory
model.LS2 <- nonparboot(model.LS2, B = 100)  # This takes a while!
SEs <- model.LS2@smoothed.mean.bsse

# PREDICTED OCCUPANCY AND STANDARD ERRORS!

psi.SEs <- data.frame(cbind(E.psi[2,], SEs[1,]))
psi.SEs
colnames(psi.SEs) <- c("OccProb", "StdErr")
psi.SEs
psi.SEs$Year <- c("05","06","07","08","09","10","11","12","13","14","15","16","17","18")
psi.SEs

# plot predicted occupancy
tiff('AmericanToad_occupancy.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')

ggplot(psi.SEs, aes(x=Year, y=OccProb)) + 
  geom_errorbar(aes(ymin=OccProb-StdErr, ymax=OccProb+StdErr), width=.1) +
  #expand_limits(y=0) +                        # Expand y range
  #scale_y_continuous(breaks=0:1 * 0.1) +
  geom_line() +
  geom_point() +
  ylab("Occupancy probability") +
  ylim(0,1) +
  xlab("Year") +
  theme_bw()

dev.off()

# Find the estimates of number of sites occupied in each year
nsites <- 30
round(E.psi*nsites)

# compare predicted occupancy in YEAR 1 (2005) to the naive occupancy in YEAR 1
# examine the simulated naive occupancy
#naive_occ <- sum(apply(AT[,4:7], 1, sum) > 0) / nrow(AT)

# predict extinction probabilities
E.ext <- predict(model.LS2, type = "ext", newdata = nd)
E.ext

# predict colonization probabilities
E.col <- predict(model.LS2, type = "col", newdata = nd)
E.col

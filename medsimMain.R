# SAMPLE CODE FOR GENERATING SIMULATED DATA FOR
# MEDIATION ANALYSIS WITH SURVIVAL OUTCOME AND TIME-UPDATED MEDIATOR
#######################################################
####################################################
## ----------------------------------------------------------------------------
##  Code to accompany chapter 8 of the thesis by Kamaryn Tanner:
##    "Dynamic Prediction, Mediation and Communication for Survival
##      Outcomes, with Applications to Cystic Fibrosis"
##  
##  Contact for code: kamaryn.tanner1@lshtm.ac.uk
## ----------------------------------------------------------------------------
##  Date        :       08-July-2021
##  R Version   :       4.0.2                                                           
##  
##  Other files :       Scenario1.R contains parameter values for simulation     
##                      medsimSupportFxns.R contains helper functions called in this file
##  
##  R packages needed: MASS, plyr, dplyr
##                      
## ------------------------------------------------------------------
#####################################################################

## ---- Load R libraries and source files ---- ##

library(MASS)
library(plyr)
library(dplyr)
source("medsimSupportFxns.R")

## ---- Load parameter settings to control the simulation ---- ##

source("Scenario1.R")

## ---- Work with sub-scenarios ---- ##
## 3 subscenarios: Both DE and IE, No DE, and No IE
subScenList <- c("Both", "NoDE", "NoIE")
truthOut <- list()

## Set the seed only once
set.seed(35711)


## ----
##  First generate the truth -- a massive number from which we calculate the DE and IE.
##  Then, separately simulate data to be run through methods of Vansteelandt et al. and 
##  Aalen et al. for study
## ----

## Number individuals in the truth dataset
numIndivT <- 5000000

## Generate data for four cases:
##  A=1 and M(A=1)
##  A=0 and M(A=0)
##  A=1 but M set to level as if A=0
##  A=0 but M set to level as if A=1
for(subScen in subScenList){
  cat(subScen)
  cat("\n")
  datS <- list()
  datS[[1]] <- data.frame(id=seq(1, numIndivT, by=1), A=1, AforM=1, Y=0, Event=0)
  datS[[2]] <- data.frame(id=seq(1, numIndivT, by=1), A=0, AforM=0, Y=0, Event=0)
  datS[[3]] <- data.frame(id=seq(1, numIndivT, by=1), A=1, AforM=0, Y=0, Event=0)
  datS[[4]] <- data.frame(id=seq(1, numIndivT, by=1), A=0, AforM=1, Y=0, Event=0)
  
  ## ---- Generate binary Z0 that affects A, Y, M  (the same for all 4 cases)
  Z0.common <- rbinom(n=numIndivT, size=1, prob=0.5)
  for(i in 1:4){ datS[[i]]$Z0 <- Z0.common}
  
  ## ---- Get random intercept and slope for mediator 
  randomM <- getRandMed( muM0=muM0, muM1=muM1, Mcorr=Mcorr, sdM0=sdM0, sdM1=sdM1, nDat=numIndivT)
  ## ---- Generate the mediator values for our 4 potential outcomes
  ##  first get baseline mediator measurement, common across 4 cases 
  M0.common <- rnorm(n=numIndivT, mean=randomM[,1]+betaZ0*Z0.common, sd=1)
  for(i in 1:4){ datS[[i]]$M0 <- M0.common}
  ##  then get the rest of the mediator values  
  datS <- lapply(datS, genMTruth, subScen = subScen)
  rm(randomM)
  rm(M0.common)
  rm(Z0.common)
  
  ## ---- Generate event times -- 
  censor_time <- visitTimes[length(visitTimes)]
  datS <- lapply(datS, genY, subScen=subScen)
  # Inspect histograms of event times for the 4 cases
  par(mfrow=c(2,2))
  histY <- function(x){ hist(x$Y)}
  lapply( datS, FUN=histY)
  
  ## To accommodate method of Aalen requirement that analysis is conditional 
  ## on surviving to the first mediator measurement
  datS <- lapply(datS, FUN=removeBeforeM1)

  ## Pick time points to evaluate survival probability as percentile of event occurrence 
  evalTimes <- quantile(datS[[1]]$Y[datS[[1]]$Event==1], probs=c(0.2, 0.5, 0.8), names=FALSE)

  ## ---- Calculate the true estimand values ---- ##
  survProbs <- lapply(datS, FUN=getSurvProb, times=evalTimes)
  
  ##SA1M1/SA0M0
  TrueTE <- survProbs[[1]]/survProbs[[2]]
  ##SA1M0/SA0M0
  TrueDE <- survProbs[[3]]/survProbs[[2]]
  ##SA1M1/SA1M0
  TrueIE <- survProbs[[1]]/survProbs[[3]]
  ##(SA1M1-SA1M0)/(SA1M1-SA0M0)
  TrueMP <- (survProbs[[1]]-survProbs[[3]])/(survProbs[[1]]-survProbs[[2]])

  ## View the survival curves
  par(mfrow=c(1,1))
  plot(c(visitTimes[2],evalTimes), c(1, survProbs[[1]]), type="l", ylim=c(0,1.05), col=2, lwd=2,
       ylab="Survivor function", xlab="Time", 
       main=paste0(scenName, " truth: " , subScen, " sub-scenario" ),
       xlim=c(0, max(visitTimes)), yaxt="none")
  axis(2, seq(0, 1.0, 0.25))
  lines(c(visitTimes[2],evalTimes), c(1, survProbs[[2]]), col=3, lwd=2)
  lines(c(visitTimes[2],evalTimes), c(1, survProbs[[3]]), col=1, lwd=2, lty=2)
  legend("topright", 
         legend=c(as.expression(bquote("S"[(A1)(M1)])), 
                  as.expression(bquote("S"[(A0)(M0)])), 
                  as.expression(bquote("S"[(A1)(M0)]))), 
         col=c(2,3,1), lty=c(1,1,2))
  
  truthOut[[subScen]] <- c(subScen=subScen, survProbs=list(survProbs), 
                           data.frame(evalTimes))
  
}

## save file
saveRDS(truthOut, file=paste0("Truth-", scenName, ".rds"))
rm(datS)
    


################################################################
## Part 2 is the generation of data for the mediation analysis methods

numIndivPerDataset <- 2200  #remember, some will be lost b/c of not surviving to 1st mediator measurement
numDatasets <- 1000
numIndiv <- numIndivPerDataset*numDatasets 
##Again, loop over the three sub-scenarios
for(subScen in subScenList){
  cat(subScen)
  cat("\n")
  datV <- data.frame(id=seq(1, numIndiv, by=1), A=rbinom(n=numIndiv,1,p=0.5), Y=0, Event=0)
  
  randomM <- getRandMed( muM0=muM0, muM1=muM1, Mcorr=Mcorr, sdM0=sdM0, sdM1=sdM1, nDat=numIndiv)
  
  ## ---- Generate a binary Z0
  datV <- genZ0(x=datV)
  ## ---- Generate the mediator values
  datV <- genMdata(x=datV, subScen=subScen)
  ## ---- Generate event times
  datV <- genY(x=datV, subScen=subScen)
  
  par(mfrow=c(1,1))
  hist(datV$Y)
  
  ## To accommodate method of Aalen requirement that analysis is conditional 
  ## on surviving to the first mediator measurement
  datV <- datV[datV$Y>=visitTimes[2],]

  ## Save data to file by sub-scenario
  #  Data is stored in one big data frame so separate code is needed to divide this up
  #  into individual simulated datasets.
  out <- list(name=paste0(scenName,subScen), visitTimes=visitTimes, evalTimes1 = truthOut[[subScen]]$evalTimes, simDat=datV)
  saveRDS(out, file=paste0("SimDat-",scenName, subScen, ".rds"))

}



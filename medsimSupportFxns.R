#######################################################
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
##  This file contains the helper functions called from medSimMain.R
## ------------------------------------------------------------------
#####################################################################



## ---- Mediator 

# First generate a random intercept and a random time slope for each
# observation via a multivariate normal distribution with mean intercept and slope
getRandMed <- function( muM0, muM1, Mcorr, sdM0, sdM1, nDat){
  
  covarM <- Mcorr * sdM0 * sdM1
  SigmaM <- matrix(c(sdM0^2, covarM, covarM, sdM1^2),2,2)
  randomM <- mvrnorm(n=nDat, mu=c(muM0, muM1), SigmaM)
  return(randomM)
}

# Generate mediator values for truth data
genMTruth <- function(x, subScen="Both"){ 
  if(subScen == "NoIE"){ ##case when no IE because A doesn't affect M
    betaAj. <- betaAjNoIE
  }else{
    betaAj. <- betaAj
  }
  # we have M0; start with M1
  cn <- paste0("M", visitTimes)
  for(i in 2:(length(cn)-1)){
    x[[cn[i]]] <- rnorm(n=nrow(x), mean=randomM[,1]+betaZ0*x$Z0+ 
                          randomM[,2]*visitTimes[i] +betaAj.[i]*x$AforM,sd=1)
  }
  return( x)
}

## Generate survival times 
## Code adapted from Keogh et al. (2021) "Simulating longitudinal data from marginal 
##      structural models using the additive hazard model"  Biometrical Journal
genY <- function(x, subScen="Both"){
  countNeg <- 0   ## additive hazards models run the risk of having a negative hazard
  T.obs <- rep(NA,nrow(x))
  x$Event <- rep(NA,nrow(x))
  if(subScen=="NoDE"){ ##change for case when no DE because A doesn't affect Y
    alphaA. <- alphaANoDE
  }else{
    alphaA. <- alphaA
  }
  for(v in 1:(length(visitTimes)-1)){
    timeInt <- visitTimes[v+1] - visitTimes[v]
    colName <- paste0("M", visitTimes[v])
    u.t <- runif(nrow(x),0,1)
    haz <- hazFac*(alpha0 + alphaA.*x$A + alpha_Mj[v]*x[[colName]] + alpha_Z0*x$Z0 )
    # Provide feedback about whether many hazards are negative
    countNeg <- length(which(haz<=0))
    # Uncomment the below line to get feedback about whether negative hazards are a problem
    #print(paste0("visit  ", v, " countNeg= ", countNeg))
    new.t <- -log(u.t)/haz
    #don't let hazard go negative 
    T.obs <- ifelse(is.na(T.obs) & new.t<timeInt & haz>0, visitTimes[v]+new.t,T.obs)
  }
  x$Event <- ifelse(is.na(T.obs),0,1)
  x$Y <- ifelse(is.na(T.obs),censor_time,T.obs)
  return(x)
}

## Removes individuals who did not survive to first mediator measurement
removeBeforeM1<-function(x){
  x <- x[x$Y>=visitTimes[2],]
  return(x)
}

##For calculating the truth
getSurvProb <- function(x, times){
  sp <- rep(0.0, length(times))
  for(i in 1:length(times)){
    sp[i] <- length(which(x$Y > times[i]))
  }
  sp <- sp/nrow(x)
  return(sp)
}

genZ0 <- function(x){ 
  x$Z0 <- if_else(x$A==1, rbinom(n=nrow(x),1,p=muZ0A1), rbinom(n=nrow(x),1,p=muZ0A0))
  return( x)
}

genMdata <- function(x, subScen="Both"){ 
  betaAj. <- betaAj
  if(subScen == "NoIE"){ ##change for case when no IE because A doesn't affect M
    betaAj. <- betaAjNoIE  
  }
  cn <- paste0("M", visitTimes)
  for(i in 1:(length(cn)-1)){
    x[[cn[i]]] <- rnorm(n=nrow(x), mean=randomM[,1]+betaZ0*x$Z0+ 
                          randomM[,2]*visitTimes[i] +betaAj.[i]*x$A,sd=1)
  }
  return( x)
}


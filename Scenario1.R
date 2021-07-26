##############################
## Parameters for the Minimal scenario
##
## This is the base scenario 
## 17 September 2020 -- M time-varying
## ## B0 affects M, Y and A
## 
##
## No negative hazards
##############################

## Name of this scenario
scenName <- "Scenario1"

## Visit times 
## Note: the last number represents the end of study when all are censored
visitTimes <- c(0, 1, 2, 3, 4) 


## ---- 1.0 Parameters for the binary baseline covariate

muZ0A1 <- 0.6  ##This is the prob of a 1 when A==1
muZ0A0 <- 0.4  ##This is the prob of a 1 when A==0


## ---- 2.0 Parameters for the mediator

## random intercept and slope control
muM0<-   2.9   # Mean random intercept 
muM1<-   0     # Mean random slope 
Mcorr<- -0.3   #Correlation between the two
sdM0<-   0.5
sdM1<-   0.1 

## mediator model looks like:
## M(t) <- draw from N(mean=rndInt+rndSlp*visitTime + betaAj(t)*A + betaZ0*Z0,sd=1)
betaAj <- c(0, 2, 2, 2)  # Effect of the exposure on the mediator at t=0,1,2,3
betaAjNoIE <- c(0,0,0,0) # For no IE, exposure doesn't affect mediator
betaZ0 <- 0.5            # Effect of Z0 on mediator

## ---- 3.0 Parameters for the survival time

## additive hazards model looks like:
## haz <- alpha0 + alphaA*A + alpha_Mj(t)*M(t) + alphaZ0*Z0
alpha0 <- 0.75    # Baseline hazard 
alphaA <- 0.35    # Effect of the exposure on the hazard
alphaANoDE <- 0.0 # For no DE, exposure doesn't affect hazard of event
alpha_Mj <- c(0.3,0.3,0.3,0.3)  # Effect of M on the hazard at t=0,1,2,3
alpha_Z0 <- 0.35  # Effect of Z0 on hazard

## Param to control absolute magnitude of hazard
##   (avoids needing to change all parameters to adjust)
hazFac <- 0.35




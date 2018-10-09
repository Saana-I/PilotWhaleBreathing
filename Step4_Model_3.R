######### Set this indicator to FALSE, if you want to run code without re-fitting 
######### the models (save to individual files, e.g., Step4_Model_3_fit_(1).Rd)

  fitModels <- FALSE

## R-tools
####################################################

  library(bbmle)
  
  ### Function to calculate coefficient of determination (CoD)
  myRsq <- function(obs, pred, type="CoD") {
    n <- length(obs)
    SStot <- sum((obs-mean(obs))^2)
    SSreg <- sum((pred-mean(obs))^2)
    SSres <- sum((obs-pred)^2)
    if(type=="CoD") {return(1-SSres/SStot)}
    if(type=="EV") {return((SSreg/n)/(SStot/n))}
  }
  
  ### Functions to do logit and inverse logit transforms
  ilogit <- function(x) {return(exp(x)/(1+exp(x)))} 
  logit <- function(p) {return(log(p/(1-p)))}
  
  ### Function to calculate expected values for instantaneous breathing rate
  E.ibr <-function(Ydata, Xdata, e0, e1, e2, e3, e4, e5, DR0, DR1, DR2, DR3, BMR0, BMR1, BMR2, BMR3, FC0, FC1, FC2, FC3, r, imax, retAll=F){ #input params can be any subset of beta1, beta2, beta3, beta4. if a parameter is not input, or is NA, the function assumes that the corresponding term should be omitted from the model.
    
    # Input check
    e0 <- ifelse(missing(e0), NA, e0)
    e1 <- ifelse(missing(e1), NA, e1)
    e2 <- ifelse(missing(e2), NA, e2)
    e3 <- ifelse(missing(e3), NA, e3)
    e4 <- ifelse(missing(e4), NA, e4)
    e5 <- ifelse(missing(e5), NA, e5)
    
    BMR0 <- ifelse(missing(BMR0), NA, BMR0)
    BMR1 <- ifelse(missing(BMR1), NA, BMR1)
    BMR2 <- ifelse(missing(BMR2), NA, BMR2)
    BMR3 <- ifelse(missing(BMR3), NA, BMR3)
    
    DR0 <- ifelse(missing(DR0), NA, DR0)
    DR1 <- ifelse(missing(DR1), NA, DR1)
    DR2 <- ifelse(missing(DR2), NA, DR2)
    DR3 <- ifelse(missing(DR3), NA, DR3)
    
    FC0 <- ifelse(missing(FC0), NA, FC0)
    FC1 <- ifelse(missing(FC1), NA, FC1)
    FC2 <- ifelse(missing(FC2), NA, FC2)
    FC3 <- ifelse(missing(FC3), NA, FC3)
    
    r <- ifelse(missing(r), NA, r)
    imax <- ifelse(missing(imax), NA, imax)
    
    # Basal cost of diving (excluding fluking cost and exposures)
    BMR <- exp(BMR0 + BMR1*Xdata$ind.calf + BMR2*Xdata$ind.small + BMR3*Xdata$ind.large)
    
    # Cost of fluking during diving (function of body size)
    FC <- exp(FC0 + FC1*Xdata$ind.calf + FC2*Xdata$ind.small + FC3*Xdata$ind.large)
    
    # Baseline dive rate (excluding effects of dive history)
    DR <- exp(DR0 + DR1*Xdata$ind.calf + DR2*Xdata$ind.small + DR3*Xdata$ind.large)
    
    # Factor increase/decrease in dive rate/breathing rate during exposures
    Ef <- exp(e0*Xdata$SIL + e1*Xdata$PB_BBN + e2*Xdata$MLFAS + e3*Xdata$PB_SON + e4*Xdata$PB_KWF + e5*Xdata$PB_KWM)
    
    # Number of breaths required to recover from the dive
    nb <- Ef*BMR*Xdata$dive.dur + FC*Xdata$dive.fluken
    
    # Baseline breathing/ dive rate
    beta0 <- Ef*DR
    
    eb <- pmax(nb-1,0)           # Remaining dive cost after first breath
    beta1 <- pmin(eb*r,imax)     # Initial breathing rate
    tr <- eb/beta1               # Time to recovery
    
    nr <- dim(Xdata)[1]
    nc <- dim(Ydata)[1]
    
    Tv <- matrix(Ydata$ST, nr, nc, byrow=T)-matrix(Xdata$exp.t, nr, nc) # Effect time (when effect starts)
    B1 <- matrix(beta1, nr, nc)
    Tr <- matrix(tr, nr, nc) 
    E <- B1*exp(-Tv/Tr)
    E <- E*(Tv>0)
    
    res <- beta0 + apply(E, 2, sum, na.rm=T)
    
    if(retAll) {
      return(list(beta0=beta0, nb=nb, beta1=beta1, tr=tr, Tv=Tv, E=E))
    } else {
      return(res)
    }
  }
  
  ### Function to calculate minus log-likelihood of a model given the data
  LL <-function(k_gamma, e0=0, e1=0, e2=0, e3=0, e4=0, e5=0, DR0,DR1=0,DR2=0,DR3=0, BMR0,BMR1=0,BMR2=0,BMR3=0, FC0,FC1=0,FC2=0,FC3=0, r, imax){
    
    # convert from working to natural
    parlist <- list(k_gamma=k_gamma, 
                    e0=e0,e1=e1,e2=e2,e3=e3,e4=e4,e5=e5,
                    DR0=DR0,DR1=DR1,DR2=DR2,DR3=DR3,
                    BMR0=BMR0,BMR1=BMR1,BMR2=BMR2,BMR3=BMR3,
                    FC0=FC0, FC1=FC1, FC2=FC2, FC3=FC3,
                    r=r, imax=imax)
    nparams <- pw2pn(parlist)
    
    whales <- unique(divedata$ind)
    D <- rep(NA, dim(Ydata)[1])
    for(w in 1:length(whales)) {
      wBool <- Xdata$ind==whales[w]
      D[wBool]<-E.ibr(Ydata=Ydata[wBool,], Xdata=Xdata[wBool,], 
                      e0=nparams$e0, e1=nparams$e1, e2=nparams$e2, e3=nparams$e3, e4=nparams$e4, e5=nparams$e5,
                      DR0=nparams$DR0, DR1=nparams$DR1, DR2=nparams$DR2, DR3=nparams$DR3,  
                      BMR0=nparams$BMR0, BMR1=nparams$BMR1, BMR2=nparams$BMR2, BMR3=nparams$BMR3, 
                      FC0=nparams$FC0, FC1=nparams$FC1, FC2=nparams$FC2, FC3=nparams$FC3, 
                      r=nparams$r, imax=nparams$imax)
    }
    
    D <- D[Xdata$dive.filter]
    ibr <- Ydata$ibr[Xdata$dive.filter]
    
    k_gamma[k_gamma==0] <- 1e-9 #to prevent log(0), /0
    D[D==0] <- 1e-9#to prevent log(0), /0
    s<-k_gamma/D
    a<-(D^2)/k_gamma
    a[a==0] <- a[a==0]+1e-9 #to prevent "value out of range in lgamma"
    a[ a<0 & is.integer(a)] <- a[ a<0 & is.integer(a)]+1e-9 #to prevent "value out of range in lgamma"
    LLike <- sum((ibr/s)-(a-1)*log(ibr)+lgamma(a)+a*log(s))
    
    return(LLike)
  }
  
  ### Function to get predicted values from the model
  PredD <- function(params, Ydata, Xdata, transform=T)
  {
    
    if(transform) {nparams <- pw2pn(params)} else {nparams <- params}
    
    if(length(nparams$e0)==0) {nparams$e0 <- 0}
    if(length(nparams$e1)==0) {nparams$e1 <- 0}
    if(length(nparams$e2)==0) {nparams$e2 <- 0}
    if(length(nparams$e3)==0) {nparams$e3 <- 0}
    if(length(nparams$e4)==0) {nparams$e4 <- 0}
    if(length(nparams$e5)==0) {nparams$e5 <- 0}
    
    whales <- unique(Xdata$ind)
    D <- rep(NA, dim(Ydata)[1])
    for(w in 1:length(whales)) {
      wBool <- Xdata$ind==whales[w]
      D[wBool]<-E.ibr(Ydata=Ydata[wBool,], Xdata=Xdata[wBool,], 
                      e0=nparams$e0, e1=nparams$e1, e2=nparams$e2, e3=nparams$e3, e4=nparams$e4, e5=nparams$e5,
                      DR0=nparams$DR0, DR1=nparams$DR1, DR2=nparams$DR2, DR3=nparams$DR3, 
                      BMR0=nparams$BMR0, BMR1=nparams$BMR1, BMR2=nparams$BMR2, BMR3=nparams$BMR3, 
                      FC0=nparams$FC0, FC1=nparams$FC1,FC2=nparams$FC2, FC3=nparams$FC3, 
                      r=nparams$r, imax=nparams$imax)
    }
    
    return(D)
  }
  
  ### Function to calculate AIC of the resulting model fit
  getAIC <- function(FittedModel) {return(2*length(FittedModel@coef)-2*FittedModel@min)}
  
  ### Function to transform natural parameters to working space
  pn2pw <- function(nparams, do=c("k_gamma", "r", "r0", "r1", "imax")) {
    if(do[1]=="all") {
      return(lapply(nparams, log))
    } else {
      tempI <- match(do,names(nparams))
      tempI <- tempI[!is.na(tempI)]
      wparams <- nparams
      for(j in 1:length(tempI)) {
        wparams[[tempI[j]]] <- log(nparams[[tempI[j]]])
      }
      return(wparams)
    }}
  
  ### Function to transform working parameters to natural space
  pw2pn <- function(wparams, do=c("k_gamma", "r", "r0", "r1", "imax")) {
    if(do[1]=="all") {
      return(lapply(wparams, exp))
    } else {
      tempI <- match(do,names(wparams))
      tempI <- tempI[!is.na(tempI)]
      nparams <- wparams
      for(j in 1:length(tempI)) {
        nparams[[tempI[j]]] <- exp(wparams[[tempI[j]]])
      }
      return(nparams)
    }}
  
  ### Function to extract results
  getResults <- function(mymodel, Xdata, returnAll=F) {
    
    wparams <- summary(mymodel)@coef[,"Estimate"]
    nparams <- pw2pn(as.list(wparams))
    SE <- summary(mymodel)@coef[,"Std. Error"] # sqrt(diag(solve(mymodel@details$hessian)))
    wparams.l <- wparams-qnorm(0.975)*SE
    wparams.u <- wparams+qnorm(0.975)*SE
    nparams.l <- pw2pn(as.list(wparams.l))
    nparams.u <- pw2pn(as.list(wparams.u))
    
    nparams2 <- nparams
    nparams2.l <- nparams.l
    nparams2.u <- nparams.u
    
    if(length(nparams2$e0)==0) {nparams2$e0 <- 0}
    if(length(nparams2$e1)==0) {nparams2$e1 <- 0}
    if(length(nparams2$e2)==0) {nparams2$e2 <- 0}
    if(length(nparams2$e3)==0) {nparams2$e3 <- 0}
    if(length(nparams2$e4)==0) {nparams2$e4 <- 0}
    if(length(nparams2$e5)==0) {nparams2$e5 <- 0}
    
    if(length(nparams2$DR1)==0) {nparams2$DR1 <- 0}
    if(length(nparams2$DR2)==0) {nparams2$DR2 <- 0}
    if(length(nparams2$DR3)==0) {nparams2$DR3 <- 0}
    
    if(length(nparams2$BMR1)==0) {nparams2$BMR1 <- 0}
    if(length(nparams2$BMR2)==0) {nparams2$BMR2 <- 0}
    if(length(nparams2$BMR3)==0) {nparams2$BMR3 <- 0}
    
    if(length(nparams2$FC1)==0) {nparams2$FC1 <- 0}
    if(length(nparams2$FC2)==0) {nparams2$FC2 <- 0}
    if(length(nparams2$FC3)==0) {nparams2$FC3 <- 0}
    
    if(length(nparams2.l$e0)==0) {nparams2.l$e0 <- 0}
    if(length(nparams2.l$e1)==0) {nparams2.l$e1 <- 0}
    if(length(nparams2.l$e2)==0) {nparams2.l$e2 <- 0}
    if(length(nparams2.l$e3)==0) {nparams2.l$e3 <- 0}
    if(length(nparams2.l$e4)==0) {nparams2.l$e4 <- 0}
    if(length(nparams2.l$e5)==0) {nparams2.l$e5 <- 0}
    
    if(length(nparams2.l$DR1)==0) {nparams2.l$DR1 <- 0}
    if(length(nparams2.l$DR2)==0) {nparams2.l$DR2 <- 0}
    if(length(nparams2.l$DR3)==0) {nparams2.l$DR3 <- 0}
    
    if(length(nparams2.l$BMR1)==0) {nparams2.l$BMR1 <- 0}
    if(length(nparams2.l$BMR2)==0) {nparams2.l$BMR2 <- 0}
    if(length(nparams2.l$BMR3)==0) {nparams2.l$BMR3 <- 0}
    
    if(length(nparams2.l$FC1)==0) {nparams2.l$FC1 <- 0}
    if(length(nparams2.l$FC2)==0) {nparams2.l$FC2 <- 0}
    if(length(nparams2.l$FC3)==0) {nparams2.l$FC3 <- 0}
    
    if(length(nparams2.u$e0)==0) {nparams2.u$e0 <- 0}
    if(length(nparams2.u$e1)==0) {nparams2.u$e1 <- 0}
    if(length(nparams2.u$e2)==0) {nparams2.u$e2 <- 0}
    if(length(nparams2.u$e3)==0) {nparams2.u$e3 <- 0}
    if(length(nparams2.u$e4)==0) {nparams2.u$e4 <- 0}
    if(length(nparams2.u$e5)==0) {nparams2.u$e5 <- 0}
    
    if(length(nparams2.u$DR1)==0) {nparams2.u$DR1 <- 0}
    if(length(nparams2.u$DR2)==0) {nparams2.u$DR2 <- 0}
    if(length(nparams2.u$DR3)==0) {nparams2.u$DR3 <- 0}
    
    if(length(nparams2.u$BMR1)==0) {nparams2.u$BMR1 <- 0}
    if(length(nparams2.u$BMR2)==0) {nparams2.u$BMR2 <- 0}
    if(length(nparams2.u$BMR3)==0) {nparams2.u$BMR3 <- 0}
    
    if(length(nparams2.u$FC1)==0) {nparams2.u$FC1 <- 0}
    if(length(nparams2.u$FC2)==0) {nparams2.u$FC2 <- 0}
    if(length(nparams2.u$FC3)==0) {nparams2.u$FC3 <- 0}
    
    BMR <- matrix(NA, length(Xdata$dive.dur), 3)
    FC <- matrix(NA, length(Xdata$dive.dur), 3)
    DR <- matrix(NA, length(Xdata$dive.dur), 3)
    Ef <- matrix(NA, length(Xdata$dive.dur), 3)
    nb <- matrix(NA, length(Xdata$dive.dur), 3)
    beta0 <- matrix(NA, length(Xdata$dive.dur),3)
    
    
    BMR[,1] <- exp(nparams2.l$BMR0 + nparams2.l$BMR1*Xdata$ind.calf + nparams2.l$BMR2*Xdata$ind.small + nparams2.l$BMR3*Xdata$ind.large)
    FC[,1] <- exp(nparams2.l$FC0 + nparams2.l$FC1*Xdata$ind.calf + nparams2.l$FC2*Xdata$ind.small + nparams2.l$FC3*Xdata$ind.large)
    DR[,1] <- exp(nparams2.l$DR0 + nparams2.l$DR1*Xdata$ind.calf + nparams2.l$DR2*Xdata$ind.small + nparams2.l$DR3*Xdata$ind.large)
    Ef[,1] <- exp(nparams2.l$e0*Xdata$SIL + nparams2.l$e1*Xdata$PB_BBN + nparams2.l$e2*Xdata$MLFAS + nparams2.l$e3*Xdata$PB_SON + nparams2.l$e4*Xdata$PB_KWF + nparams2.l$e5*Xdata$PB_KWM)
    nb[,1] <- (Ef[,1]*BMR[,1])*Xdata$dive.dur + FC[,1]*Xdata$dive.fluken
    beta0[,1] <- Ef[,1]*DR[,1]
    
    BMR[,2] <- exp(nparams2$BMR0 + nparams2$BMR1*Xdata$ind.calf + nparams2$BMR2*Xdata$ind.small + nparams2$BMR3*Xdata$ind.large)
    FC[,2] <- exp(nparams2$FC0 + nparams2$FC1*Xdata$ind.calf + nparams2$FC2*Xdata$ind.small + nparams2$FC3*Xdata$ind.large)
    DR[,2] <- exp(nparams2$DR0 + nparams2$DR1*Xdata$ind.calf + nparams2$DR2*Xdata$ind.small + nparams2$DR3*Xdata$ind.large)
    Ef[,2] <- exp(nparams2$e0*Xdata$SIL + nparams2$e1*Xdata$PB_BBN + nparams2$e2*Xdata$MLFAS + nparams2$e3*Xdata$PB_SON + nparams2$e4*Xdata$PB_KWF + nparams2$e5*Xdata$PB_KWM)
    nb[,2] <- (Ef[,2]*BMR[,2])*Xdata$dive.dur + FC[,2]*Xdata$dive.fluken
    beta0[,2] <- Ef[,2]*DR[,2]
    
    BMR[,3] <- exp(nparams2.u$BMR0 + nparams2.u$BMR1*Xdata$ind.calf + nparams2.u$BMR2*Xdata$ind.small + nparams2.u$BMR3*Xdata$ind.large)
    FC[,3] <- exp(nparams2.u$FC0 + nparams2.u$FC1*Xdata$ind.calf + nparams2.u$FC2*Xdata$ind.small + nparams2.u$FC3*Xdata$ind.large)
    DR[,3] <- exp(nparams2.u$DR0 + nparams2.u$DR1*Xdata$ind.calf + nparams2.u$DR2*Xdata$ind.small + nparams2.u$DR3*Xdata$ind.large)
    Ef[,3] <- exp(nparams2.u$e0*Xdata$SIL + nparams2.u$e1*Xdata$PB_BBN + nparams2.u$e2*Xdata$MLFAS + nparams2.u$e3*Xdata$PB_SON + nparams2.u$e4*Xdata$PB_KWF + nparams2.u$e5*Xdata$PB_KWM)
    nb[,3] <- (Ef[,3]*BMR[,3])*Xdata$dive.dur + FC[,3]*Xdata$dive.fluken
    beta0[,3] <- Ef[,3]*DR[,3]
    
    if(returnAll) {
      return(list(
        nparams=nparams2,
        nparams.l=nparams2.l,
        nparams.u=nparams2.u,
        BMR=BMR,
        FC=FC,
        DR=DR,
        Ef=Ef,
        nb=nb,
        beta0=beta0))
    } else {
      return(list(
        nparams=nparams,
        nparams.l=nparams.l,
        nparams.u=nparams.u,
        BMR=BMR,
        FC=FC,
        DR=DR,
        Ef=Ef,
        nb=nb,
        beta0=beta0))
    }
  }
  

## Data
####################################################

  divedata <- read.csv("divedata.csv")

  ### Transform data
  divedata$ind <- as.character(divedata$ind)
  divedata$ind.type <- as.character(divedata$ind.type)
  divedata$ind.col <- as.character(divedata$ind.col)
  divedata$Session <- as.character(divedata$Session)
  divedata$dive.GMTtime <- as.POSIXct(as.character(divedata$dive.GMTtime), tz="GMT")
  
  divedata$dive.rate <- 1/(divedata$dive.dur/60)
  divedata$ind.calf <- as.numeric(divedata$ind.calf)
  divedata$ind.small <- as.numeric(divedata$ind.small)
  divedata$ind.large <- as.numeric(divedata$ind.large)

  ### Exclude exposures and non-focal animals
  divedata$dive.MLFAS <- divedata$dive.MFAS + divedata$dive.LFAS
  bBool <- divedata$dive.SIL+divedata$dive.PB_BBN+divedata$dive.MLFAS+divedata$dive.PB_SON+divedata$dive.PB_KWF+divedata$dive.PB_KWM+divedata$dive.PB_HW
  bBool <- as.numeric(bBool)==0  
  fBool <- divedata$ind!="gm09_137c" & divedata$ind!="gm09_138b" & divedata$ind!="gm13_169b" & divedata$ind!="gm14_180b"
  
  ### Set up data for the model
  
  # Response data
  Ydata <- data.frame(
    ibr=1/(divedata$dive.dur/60), # instantaneous breathing rate
    ST=divedata$dive.sfromtot/60) # start time, seconds from deployment start time
  
  # Explanatory data
  Xdata <- data.frame(
    exp.t=divedata$dive.end/60, # surface time at the end of the dive
    dive.dur=divedata$dive.dur/60, # dive duration
    dive.fluken=divedata$dive.fluken, # number of fluke strokes
    ind=divedata$ind, # deployment id number
    ind.calf=divedata$ind.calf, # association with a calf (0/1)
    ind.small=divedata$ind.small, # small body size (0/1)
    ind.large=divedata$ind.large, # large body size (0/1)
    surf.dur=divedata$surf.dur/60, # any logging/surface period following the dive
    SIL=as.numeric(divedata$dive.SIL>0.5), # no-sonar approach (0/1)
    PB_BBN=as.numeric(divedata$dive.PB_BBN>0.5), # broad-band control playback
    MLFAS=as.numeric(divedata$dive.MLFAS>0.5), # MFAS/LFAS towed sonar approach
    PB_SON=as.numeric(divedata$dive.PB_SON>0.5), # playback of sonar (LFAS) sounds
    PB_KWF=as.numeric(divedata$dive.PB_KWF>0.5), # fish-feeding killer whale playback
    PB_KWM=as.numeric(divedata$dive.PB_KWM>0.5), # mammal-feeding killer whale playback
    dive.filter=divedata$dive.depth < 31 & bBool & fBool) # filter for response data
  
## Model fitting
####################################################

  #### Make a list of model structures - in each case several parameters
  #### are fixed to zero in the model and not estimated
  
  fixeverytime <- c("FC1","FC2","FC3","e0", "e1", "e2", "e3", "e4", "e5")
  
  mlist <- list(
    
    # Diving cost
    
    # 1 - fix nothing
    fixeverytime,
    
    # 2 - fix everything 
    c("BMR1","BMR2", "BMR3", fixeverytime),
    
    # 3 - include exposures + calf effect on DMR
    c("BMR2","BMR3",fixeverytime), 
    # 4 - include exposures + small body size effect on DMR
    c("BMR1","BMR3", fixeverytime), 
    # 5 - include exposures + large body size effect on DMR
    c("BMR1", "BMR2", fixeverytime),
    
    # 6 - include exposures + calf & small body size effect on DMR
    c("BMR3",fixeverytime), 
    # 7 - include exposures + calf & large body size effect on DMR
    c("BMR2", fixeverytime), 
    # 8 - include exposures + small & large body size effect on DMR
    c("BMR1", fixeverytime))
  
  
  if(fitModels) {
    
  for(m in 1:length(mlist)) {
    
  myfix <- mlist[[m]]
  
  MLM_list <- list()
  nparams_list <- list()
  fit_mlik <- 0
  fit_convergence <- 0
  fit_time <- 0
  
  for(j in 1:5) {
    
      if(!is.null(myfix)) {
        fixlist <- as.list(rep(0,length(myfix)))
        names(fixlist) <- myfix
      } else {fixlist <- NULL}
    
      nparams <- list()
      nparams$k_gamma <- rgamma(1,5,1)

      nparams$BMR0 <- rnorm(1, mean=0, sd=1)
      nparams$BMR1 <- 0
      nparams$BMR2 <- 0
      nparams$BMR3 <- 0
      
      nparams$FC0 <- rnorm(1, mean=0, sd=2)
      nparams$FC1 <- 0
      nparams$FC2 <- 0
      nparams$FC3 <- 0
      
      nparams$DR0 <- rnorm(1, mean=0, sd=1)
      nparams$DR1 <- 0
      nparams$DR2 <- 0
      nparams$DR3 <- 0
      
      nparams$r <- rgamma(1,1,1)
      if(any(myfix=="imax")) {nparams$imax <- 100} else {nparams$imax <- runif(1,10,20)}
      
      if(any(myfix=="e0")) {nparams$e0 <- 0} else {nparams$e0 <- rnorm(1, mean=0, sd=1)}
      if(any(myfix=="e1")) {nparams$e1 <- 0} else {nparams$e1 <- rnorm(1, mean=0, sd=1)}
      if(any(myfix=="e2")) {nparams$e2 <- 0} else {nparams$e2 <- rnorm(1, mean=0, sd=1)}
      if(any(myfix=="e3")) {nparams$e3 <- 0} else {nparams$e3 <- rnorm(1, mean=0, sd=1)}
      if(any(myfix=="e4")) {nparams$e4 <- 0} else {nparams$e4 <- rnorm(1, mean=0, sd=1)}
      if(any(myfix=="e5")) {nparams$e5 <- 0} else {nparams$e5 <- rnorm(1, mean=0, sd=1)}
      
      
    nparams_list[[j]] <- nparams
    wparams <- pn2pw(nparams)
    
    temp <- Sys.time()
    MLM_list[[j]] <- try(list(all.full = mle2(LL, start = wparams, 
                                              fixed=fixlist,
                                               method = "L-BFGS-B", control=list(maxit=300))),
                          silent=T)
    temp[2] <- Sys.time()
    fit_time[j] <- diff(temp) # 8.344126 hours
    
    if(class(MLM_list[[j]])!="try-error") {
      fit_mlik[j] <- MLM_list[[j]]$all.full@min
      fit_convergence[j] <- MLM_list[[j]]$all.full@details$convergence
      print(paste("#", j, "convergence: ", fit_convergence[j]==0))
      print(paste(round(MLM_list[[j]]$all.full@min,2)))
    } else {
      fit_mlik[j] <- NA
      fit_convergence[j] <- NA
      print(paste("#", j, "failed"))}
    
    save(MLM_list, nparams_list, fit_time, fit_mlik, fit_convergence, E.ibr,
         file=paste("Step4_Model_3_fit_(",m,").Rd",sep=""))
    
  }}}


## Summarize each model and conduct model selection
####################################################
  
  nm <- length(mlist) # number of models
  
  mtab <- data.frame(ID=1:nm)
  mtab$runs <- NA
  mtab$conv <- NA
  mtab$mlik_min <- NA
  mtab$mlik_max <- NA
  mtab$npars <- NA
  mtab$AIC <- NA
  
  mtab$CoD <- NA
  mtab$EV <- NA
  mtab$corr <- NA
  
  mtab$drop_num <- NA
  mtab$drop_par <- NA
  mtab$drop_pvl <- NA
  
  mtab$k_gamma <- NA
  mtab$r <- NA
  mtab$imax <- NA
  
  mtab$FC0 <- NA
  
  mtab$DR0 <- NA
  mtab$DR1 <- NA
  mtab$DR2 <- NA
  mtab$DR3 <- NA
  
  mtab$BMR0 <- NA
  mtab$BMR1 <- NA
  mtab$BMR2 <- NA
  mtab$BMR3 <- NA
  
  ## Model predictions 
  
  # Set up dataframe for values used for prediction
  preddata <- data.frame(
    ind.type=c("S","M","L","SC","MC","LC"),
    ind.calf=c(0,0,0,1,1,1),
    ind.small=c(1,0,0,1,0,0),
    ind.large=c(0,0,1,0,0,1))
  
  preddata$SIL <- 0
  preddata$PB_BBN <- 0
  preddata$MLFAS <- 0
  preddata$PB_SON <- 0
  preddata$PB_KWF <- 0
  preddata$PB_KWM <- 0
  
  preddata$dive.dur <- 0
  preddata$dive.fluken <- 0
  
  preddata1 <- preddata
  preddata1$dive.dur <- 1
  
  preddata2 <- preddata
  preddata2$dive.fluken <- 1
  
  ## Initialize arrays to store predictions
  DC <- array(NA, c(length(mtab$ID), 6, 3)) # DC cost - lower, mean, upper
  FC <- array(NA, c(length(mtab$ID), 6, 3)) # FC cost - lower, mean upper
  TR <- array(NA, c(length(mtab$ID), 6, 3)) # Time to recovery - lower, mean upper
  NB <- array(NA, c(length(mtab$ID), 6, 4)) # 25%, 75% and max dive cost
  NR <- matrix(NA, length(mtab$ID), length(Ydata$ibr)) # Non-recovery periods (boolean)
  NDR <- array(NA, c(length(mtab$ID), 6, 3)) # 2.5%, 50%, and 97.5% quantiles non-recovery rate


  for(m in 1:nm) {
    
    load(paste("Step4_Model_3_fit_(",m,").Rd",sep=""))
    temp <- fit_mlik
    temp[is.na(fit_convergence) | fit_convergence==1] <- max(fit_mlik)
    MLM <- MLM_list[[which.min(temp)]]
    print(paste("Model", m))
    print(paste("Convergence",MLM$all.full@details$convergence))
    print(MLM$all.full@details$message)
    
    mi <- mtab$ID==m
    mtab$runs[mi] <- length(fit_mlik)
    mtab$conv[mi] <- sum(fit_convergence==0, na.rm=T)/length(fit_mlik)
    mtab$mlik_max[mi] <- max(fit_mlik, na.rm=T)
    mtab$mlik_min[mi] <- min(fit_mlik, na.rm=T)
    mtab$npars[mi] <- length(MLM$all.full@coef)
    mtab$AIC[mi] <- AIC(MLM$all.full) # 2*length(MLM$all.full@coef)+2*MLM$all.full@details$value
    
    tObj <- getResults(MLM$all.full, Xdata, returnAll=T)
    nparams<-tObj$nparams 
    D <- PredD(nparams, Ydata, Xdata, transform=F)
    
    nb <- tObj$nb[,2]
    
    mtab$CoD[mi] <- myRsq(Ydata$ibr, D, type="CoD")
    mtab$EV[mi] <- myRsq(Ydata$ibr, D, type="EV")
    mtab$corr[mi] <- cor(Ydata$ibr,D) 
    
    plot(Ydata$ibr, D, main=mtab$ID[m], xlab="Observed", ylab="Predicted", col=Xdata$dive.filter)
    abline(0,1,col="red")
    abline(v=nparams$imax, col="blue")
    abline(h=nparams$imax, col="blue")
    
    # Realized diving costs
    NB[mi,,1] <- tapply(nb[nb>1], divedata$ind.type[nb>1], quantile, 0.025)
    NB[mi,,2] <- tapply(nb[nb>1], divedata$ind.type[nb>1], quantile, 0.5)
    NB[mi,,3] <- tapply(nb[nb>1], divedata$ind.type[nb>1], quantile, 0.975)
    NB[mi,,4] <- tapply(nb[nb>1], divedata$ind.type[nb>1], max, na.rm=T)
    
    # Non-dive recovery rates
    NR[mi,] <- nb<1 & round(D,1)==round(tObj$beta0[,2],1)
    NDR[mi,,1] <- tapply(divedata$dive.rate[NR[mi,]], divedata$ind.type[NR[mi,]], quantile, 0.025)
    NDR[mi,,2] <- tapply(divedata$dive.rate[NR[mi,]], divedata$ind.type[NR[mi,]], quantile, 0.5)
    NDR[mi,,3] <- tapply(divedata$dive.rate[NR[mi,]], divedata$ind.type[NR[mi,]], quantile, 0.975)
    
    # How many p-values above 2% sig. level? Which parameter has the worst p-value?
    pars <- names(nparams) #names(MLM$all.full@coef)
    pvals <- summary(MLM$all.full)@coef[,c("Pr(z)")]
    mtab$drop_num[mi] <- sum(pvals>0.02)
    mtab$drop_par[mi] <- pars[which.max(pvals)]
    mtab$drop_pvl[mi] <- max(pvals)
    
    # Store parameter estimates (natural parameter space)
    pars_est <- names(MLM$all.full@coef)
    pars_tab <- names(mtab)
    pars_list <- names(nparams)
    FI <- which(pars_tab=="k_gamma")
    for(p in FI:length(pars_tab)) {
      if(sum(pars_tab[p]==pars_est)>0) {
        mtab[mi,p] <- nparams[[which(pars_tab[p]==pars_list)]] # MLM$all.full@coef[p]
      }
    }
    
    # Store predictions for diving and fluking cost for each individual
    
    tObj <- getResults(MLM$all.full, preddata1, returnAll=T)
    
    for(ci in 1:3) { # Loop over mean, lower and upper
      
      DC[mi,,ci] <- tObj$BMR[,ci]
      FC[mi,,ci] <- tObj$FC[,ci]
      
      if(!is.na(DC[mi,1,ci]) & !is.na(FC[mi,1,ci])) {
        
        # Time to recovery
        
        for(w in 1:6) {
          dive.dur <-10
          eb <- pmax(dive.dur*tObj$BMR[w,ci]-1,0)
          if(eb>0) {
            beta1 <- pmin(tObj$nparams$r*eb,tObj$nparams$imax)
            tr <- eb/beta1
            tt <- 0
            nn <- 0
            for(k in 2:100) {
              y <-beta1*exp(-tt[k-1]/tr)
              tt[k] <- tt[k-1] + 1/y
              nn[k] <- nn[k-1]+1
            }
            trec <- approx(x=nn, y=tt, xout=eb)$y

            TR[mi,w,ci] <- trec
          } else {
            TR[mi,w,ci] <- 0
          }
        }
      }
    }
    
  }
  
  mtab$dAIC <- mtab$AIC-min(mtab$AIC)
  View(mtab[order(mtab$AIC),])
  
  
## Interpret best model
####################################################
  
  load("Step4_Model_3_fit_(6).Rd")
  temp <- fit_mlik
  temp[is.na(fit_convergence) | fit_convergence==1 | fit_convergence==52] <- max(fit_mlik)
  MLM <- MLM_list[[which.min(temp)]] # select best converged model
  
  ### Predicted IBR for the whole data set (including data that were excluded for model fitting)
  
    tObj <- getResults(MLM$all.full, Xdata, returnAll=T)
    nparams<-tObj$nparams 
    nparams.l <- tObj$nparams.l
    nparams.u <- tObj$nparams.u
    D <- PredD(nparams, Ydata, Xdata, transform=F) 
    nb <- tObj$nb[,2] # diving cost of each IBI (breaths per IBI)
    
    myRsq(Ydata$ibr[Xdata$dive.filter], D[Xdata$dive.filter]) # 0.3058543
    tBool <- divedata$dive.depth < 31 & !fBool & bBool
    myRsq(Ydata$ibr[tBool], D[tBool]) # 0.3012832 -> similar R2 for non-focal whales during baseline
  
    save(D, nb, file="Step4_Model_3_expected_IBR.Rd")
    
  ### Get model estimates
  
  summary(MLM$all.full) # Note: all parameter estimates in log-scale (working parameter values)
  #           Estimate Std. Error  z value     Pr(z)    
  #  k_gamma  2.655938   0.043960   60.4165 < 2.2e-16 *** # Observation error
  #  DR0      1.061950   0.010532  100.8292 < 2.2e-16 *** # Intercept for non-recovery rate (NRR)
  #  DR1      0.058990   0.013336    4.4233 9.721e-06 *** # Calf effect for NRR
  #  DR2     -0.085063   0.022212   -3.8296 0.0001284 *** # Small body size effect for NRR
  #  DR3     -0.065313   0.010834   -6.0283 1.657e-09 *** # Large body size effect for NRR
  #  BMR0    -1.199130   0.077113  -15.5504 < 2.2e-16 *** # Intercept for diving cost (DC)
  #  BMR1    -4.613646  19.840195   -0.2325 0.8161183     # Calf effect for NRR
  #  BMR2     1.185174   0.119565    9.9124 < 2.2e-16 *** # Small body size effect for DC
  #  FC0     -2.457146   0.019891 -123.5324 < 2.2e-16 *** # Fluke stroke cost
  #  r       -1.369986   0.022093  -62.0112 < 2.2e-16 *** # Increase in initial breathing rate
  #  imax     1.693292   0.025930   65.3028 < 2.2e-16 *** # Maximum initial breathing rate


  ### Transform and interpret
  n <- length(summary(MLM$all.full)@coef[,1])
  mu <- unlist(nparams)[1:n]
  l <- unlist(nparams.l)[1:n]
  u <- unlist(nparams.u)[1:n]
  
  # Fluke stroke cost
  exp(mu[names(mu)=="FC0"]) # 0.08567912 
  exp(l[names(mu)=="FC0"]) # 0.08240319
  exp(u[names(mu)=="FC0"]) # 0.08908529
  
  # Diving cost for small individuals
  exp(sum(mu[names(mu)=="BMR0" | names(mu)=="BMR2"])) #  0.9861411
  exp(sum(l[names(mu)=="BMR0" | names(mu)=="BMR2"])) # 0.670698
  exp(sum(u[names(mu)=="BMR0" | names(mu)=="BMR2"])) # 1.449944
  
  # NRR for small individuals
  exp(sum(mu[names(mu)=="DR0" | names(mu)=="DR2"])) # 2.656173
  exp(sum(l[names(mu)=="DR0" | names(mu)=="DR2"])) # 2.49106
  exp(sum(u[names(mu)=="DR0" | names(mu)=="DR2"])) #  2.832229
  
  # Diving cost for medium individuals
  exp(mu[names(mu)=="BMR0"]) # 0.3014565
  exp(l[names(mu)=="BMR0"]) # 0.2591709
  exp(u[names(mu)=="BMR0"]) #  0.3506413
  
  # NRR for medium individuals
  exp(mu[names(mu)=="DR0"]) # 2.892004
  exp(l[names(mu)=="DR0"]) # ) # 2.832917
  exp(u[names(mu)=="DR0"]) #  2.952323
  
  # Effect of calf on NRR (factor increase)
  exp(mu[names(mu)=="DR1"]) # 1.060765 
  exp(l[names(mu)=="DR1"]) # 1.033397 
  exp(u[names(mu)=="DR1"]) # 1.088858 
  
  ### What proportion of diving cost was explained by locomotion?
  
    # Predict diving cost excluding the cost of fluke strokes
    Xdata2 <- Xdata
    Xdata2$dive.fluken <- 0
    tObj2 <- getResults(MLM$all.full, Xdata2, returnAll=T)
    nb2 <- tObj2$nb[,2]
    tBool <- divedata$dive.depth > 31 & bBool & fBool
    tprop <- (1-nb2[tBool]/nb[tBool]) 
    summary(tprop) # proportion of net diving cost explained by locomotion 
    # (given other characteristics of the the dataset, i.e., dive duration and individuals)
    # Min.  1st Qu.   Median  Mean   3rd Qu.   Max. 
    # 0.2058  0.8298  0.8713  0.8362  0.9975  0.9990 
    

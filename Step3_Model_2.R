## R-tools
####################################################

  library(mgcv)
  library(geepack)
  library(MuMIn)
  library(mvtnorm)

## Data
####################################################
  
  load("Step1_dive_recovery_data.Rd")
  divedata <- divedata2
  remove(divedata2)
  
  divedata$ind2 <- as.character(divedata$ind)
  divedata$ind2[divedata$ind=="gm09_137c"] <- "gm09_137b"
  divedata$ind2[divedata$ind=="gm09_138b"] <- "gm09_138a"
  divedata$ind2[divedata$ind=="gm13_169b"] <- "gm13_169a"
  divedata$ind2[divedata$ind=="gm14_180b"] <- "gm14_180a"
  divedata$ind2 <- as.factor(divedata$ind2)
  length(unique(divedata$ind)) # 17
  length(unique(divedata$ind2)) # 13
  
  divedata$ind_num <- as.numeric(as.factor(divedata$ind))
  divedata$ind.med <- as.numeric(divedata$ind.small==0 & divedata$ind.large==0)
  
## Sensitivity analysis on time window selection
####################################################
  
  mtab <- data.frame(twin=seq(5,20,0.5))
  mtab$n <- NA
  mtab$pd <- NA
  mtab$Rsq <- NA
  mtab$Rsq2 <- NA
  
  mtab$mu_deep.logging_tw <- NA
  mtab$mu_deep.postdur_tw <- NA
  mtab$mu_deep.dur <- NA
  mtab$mu_deep.dur2 <- NA
  mtab$mu_deep.fluken <- NA
  
  mtab$se_deep.logging_tw <- NA
  mtab$se_deep.postdur_tw <- NA
  mtab$se_deep.dur <- NA
  mtab$se_deep.dur2 <- NA
  mtab$se_deep.fluken <- NA
  
  mtab$deep.logging_tw <- NA
  mtab$deep.postdur_tw <- NA
  mtab$deep.dur <- NA
  mtab$deep.dur2 <- NA
  mtab$deep.fluken <- NA
  
  
  fitNow <- T
  for(m in 1:length(mtab$twin)) {
  
  divedata$deep.surfn_tw <- NA
  divedata$deep.postdur_tw <- NA
  divedata$deep.logging_tw <- NA
  divedata$deep.lognum_tw <- NA
  divedata$deep.Session_tw <- NA
  
  tempI <- unique(divedata$deep.index[!is.na(divedata$deep.index)])
  
  twin <- mtab$twin[m]
  
  for(j in 1:length(tempI)) {
    
    k <- which(divedata$deep.filter & divedata$deep.index==tempI[j])
    
    # Filter for dive cycle
    tBool <- divedata$deep.index==tempI[j] & !is.na(divedata$deep.index)
    # Filter for all post
    pBool <- tBool & divedata$deep.post
    # Filter for post time window
    pBool_tw <- pBool & (divedata$deep.tsince < (twin))
    
    divedata$deep.surfn_tw[tBool] <- sum(pBool_tw)
    divedata$deep.postdur_tw[tBool] <- pmin(twin, divedata$deep.postdur[k]) # in minutes
    divedata$deep.logging_tw[tBool] <- sum(divedata$surf.dur*pBool_tw*divedata$logging)# in seconds

    # How much logging period exceeds end of 5 min interval?
    if(sum(divedata$surf.dur*pBool_tw*divedata$logging)>0) {
      
      loggingT <- divedata$surf.dur[pBool_tw & divedata$logging]
      ST <- divedata$dive.end[k]
      # end time of last breath/surface interval
      ET_temp <- max(divedata$dive.end[pBool_tw & divedata$logging]+divedata$surf.dur[pBool_tw & divedata$logging])
      if(ET_temp > (ST+5*60)) {
        print(paste("correcting logging for # ", k))
        loggingT[length(loggingT)] <- loggingT[length(loggingT)]-(ET_temp-(ST+twin*60))
      }
      divedata$deep.logging_tw[tBool] <- sum(loggingT)
      
      # estimated # breaths during logging
      smoothIBI <- runmed(divedata$dive.dur[pBool_tw],k=3,endrule="median")
      #preIBI <- divedata$dive.dur[pBool_tw & divedata$logging]
      preIBI <- smoothIBI[divedata$logging[pBool_tw]]
      divedata$deep.lognum_tw[tBool] <- sum(loggingT/preIBI)
    } else {
      divedata$deep.lognum_tw[tBool] <- 0
    }
    divedata$deep.logging_tw[tBool] <- divedata$deep.logging_tw[tBool]/60
    
    sessions <- unique(c(divedata$Session[k],divedata$Session[pBool_tw]))
    exposures <- sessions[sessions!="Baseline" & sessions!="Post"]
    if(length(exposures)>0) {divedata$deep.Session_tw[tBool] <- exposures[1]
    } else {divedata$deep.Session_tw[tBool] <- sessions[1]}
    
    
  }
  
  filter3 <- divedata$deep.first & divedata$deep.postdur > (40/60)
  filter3 <- filter3 & (divedata$deep.Session_tw=="Baseline" | divedata$deep.Session_tw=="Post")
  filter3 <- filter3 & !is.na(divedata$deep.surfn_tw) & !is.na(divedata$deep.logging_tw) & !is.na(divedata$deep.postdur_tw) 
  filter3 <- filter3 & !is.na(divedata$deep.dur) & !is.na(divedata$deep.fluken)
  
  mtab$pd[m] <- 1-sum(divedata$deep.postdur_tw[filter3]==twin)/sum(filter3)

    if(fitNow) {
  
      fit1 <- try(gamm(deep.surfn_tw ~ deep.lognum_tw + deep.postdur_tw + deep.dur,
                   family=quasipoisson(link = "identity"),
                   data=divedata[filter3,],
                   random=list(ind2=~1, deep.index=~1),
                   correlation=corAR1()),silent=T)
    
      if(class(fit1)[1]!="try-error") {
        
        temp <- summary(fit1$gam)
      
        mtab$n[m] <- sum(filter3)
        mtab$Rsq[m] <- temp$r.sq
        
        # store estimates
        mtab$mu_deep.logging_tw[m] <- temp$p.table["deep.lognum_tw",1] #temp$p.table["deep.logging_tw",1]
        mtab$mu_deep.postdur_tw[m] <- temp$p.table["deep.postdur_tw",1]
        mtab$mu_deep.dur[m] <- temp$p.table["deep.dur",1]
        
        # store standard errors
        mtab$se_deep.logging_tw[m] <- temp$p.table["deep.lognum_tw",2] #temp$p.table["deep.logging_tw",2]
        mtab$se_deep.postdur_tw[m] <- temp$p.table["deep.postdur_tw",2]
        mtab$se_deep.dur[m] <- temp$p.table["deep.dur",2]
        
        # store test statistics
        mtab$deep.logging_tw[m] <- temp$p.table["deep.lognum_tw",3] #temp$p.table["deep.logging_tw",3]
        mtab$deep.postdur_tw[m] <- temp$p.table["deep.postdur_tw",3]
        mtab$deep.dur[m] <- temp$p.table["deep.dur",3]
      }
      
      fit2 <- try(gamm(deep.surfn_tw ~ deep.lognum_tw + deep.postdur_tw + deep.dur + deep.fluken,
                       family=quasipoisson(link ="identity"),
                       data=divedata[filter3,],
                       random=list(ind2=~1, deep.index=~1),
                       correlation=corAR1()),silent=T)
      
      if(class(fit2)[1]!="try-error") {
        
        temp <- summary(fit2$gam)
        
        mtab$Rsq2[m] <- temp$r.sq
        
        mtab$mu_deep.fluken[m] <- temp$p.table["deep.fluken",1]
        mtab$se_deep.fluken[m] <- temp$p.table["deep.fluken",2]
        mtab$deep.fluken[m] <- temp$p.table["deep.fluken",3]
        
        mtab$mu_deep.dur2[m] <- temp$p.table["deep.dur",1]
        mtab$se_deep.dur2[m] <- temp$p.table["deep.dur",2]
        mtab$deep.dur2[m] <- temp$p.table["deep.dur",3]
      }
    }
  }
  
  ### Plot coefficients as a function of time window
  
  twin_final <- 10 # selected window
  par(mfrow=c(2,2), mar=c(4,4,3,1))
  
  plot(mtab$twin, mtab$mu_deep.postdur_tw, type="l", ylim=c(1.2,2.6),
       xlab="Post-dive period (min)", ylab="Breaths per PDSI min",
       main="a. Post-dive duration (min)")
  lines(mtab$twin, mtab$mu_deep.postdur_tw-2*mtab$se_deep.postdur_tw, lty=2)
  lines(mtab$twin, mtab$mu_deep.postdur_tw+2*mtab$se_deep.postdur_tw, lty=2)
  abline(v=twin_final, col="blue")
  
  plot(mtab$twin, mtab$mu_deep.logging_tw, type="l", ylim=c(-0.62,0),
       xlab="Post-dive period (min)", ylab="Breaths per missed breath",
       main="b. Missed breaths")
  lines(mtab$twin, mtab$mu_deep.logging_tw-2*mtab$se_deep.logging_tw, lty=2)
  lines(mtab$twin, mtab$mu_deep.logging_tw+2*mtab$se_deep.logging_tw, lty=2)
  abline(v=twin_final, col="blue")
  
  plot(mtab$twin, mtab$mu_deep.dur, type="l", ylim=c(1.4,2.3),
       xlab="Post-dive period (min)", ylab="Breaths per diving min",
       main="c. Dive duration (min)")
  lines(mtab$twin, mtab$mu_deep.dur-2*mtab$se_deep.dur, lty=2)
  lines(mtab$twin, mtab$mu_deep.dur+2*mtab$se_deep.dur, lty=2)
  abline(v=twin_final, col="blue")
  
  plot(mtab$twin, mtab$mu_deep.fluken, type="l", ylim=c(0.03,0.11),
       xlab="Post-dive period (min)", ylab="Breaths per stroke",
       main="d. Number of strokes")
  lines(mtab$twin, mtab$mu_deep.fluken-2*mtab$se_deep.fluken, lty=2)
  lines(mtab$twin, mtab$mu_deep.fluken+2*mtab$se_deep.fluken, lty=2)
  abline(v=twin_final, col="blue")

  
  
## Select time window and re-do data for final analysis
####################################################
  
  divedata$deep.surfn_tw <- NA
  divedata$deep.postdur_tw <- NA
  divedata$deep.logging_tw <- NA
  divedata$deep.lognum_tw <- NA
  divedata$deep.Session_tw <- NA
  
  twin <- 10
  tempI <- unique(divedata$deep.index[!is.na(divedata$deep.index)])
  
  for(j in 1:length(tempI)) {
    
    k <- which(divedata$deep.filter & divedata$deep.index==tempI[j])
    
    # Filter for dive cycle
    tBool <- divedata$deep.index==tempI[j] & !is.na(divedata$deep.index)
    # Filter for all post
    pBool <- tBool & divedata$deep.post
    # Filter for post time window
    pBool_tw <- pBool & (divedata$deep.tsince < (twin))
    
    divedata$deep.surfn_tw[tBool] <- sum(pBool_tw)
    divedata$deep.postdur_tw[tBool] <- pmin(twin, divedata$deep.postdur[k]) # in minutes
    divedata$deep.logging_tw[tBool] <- sum(divedata$surf.dur*pBool_tw*divedata$logging)# in seconds
    
    # How much logging period exceeds end of 5 min interval?
    if(sum(divedata$surf.dur*pBool_tw*divedata$logging)>0) {
      
      loggingT <- divedata$surf.dur[pBool_tw & divedata$logging]
      ST <- divedata$dive.end[k]
      # end time of last breath/surface interval
      ET_temp <- max(divedata$dive.end[pBool_tw & divedata$logging]+divedata$surf.dur[pBool_tw & divedata$logging])
      if(ET_temp > (ST+5*60)) {
        print(paste("correcting logging for # ", k))
        loggingT[length(loggingT)] <- loggingT[length(loggingT)]-(ET_temp-(ST+twin*60))
      }
      divedata$deep.logging_tw[tBool] <- sum(loggingT)
      
      # estimated # breaths during logging
      smoothIBI <- runmed(divedata$dive.dur[pBool_tw],k=3,endrule="median")
      #preIBI <- divedata$dive.dur[pBool_tw & divedata$logging]
      preIBI <- smoothIBI[divedata$logging[pBool_tw]]
      divedata$deep.lognum_tw[tBool] <- sum(loggingT/preIBI)
    } else {
      divedata$deep.lognum_tw[tBool] <- 0
    }
    divedata$deep.logging_tw[tBool] <- divedata$deep.logging_tw[tBool]/60
    
    sessions <- unique(c(divedata$Session[k],divedata$Session[pBool_tw]))
    exposures <- sessions[sessions!="Baseline" & sessions!="Post"]
    if(length(exposures)>0) {divedata$deep.Session_tw[tBool] <- exposures[1]
    } else {divedata$deep.Session_tw[tBool] <- sessions[1]}
    
    
  }
  
  filter3 <- divedata$deep.first & divedata$deep.postdur > (40/60)
  filter3 <- filter3 & (divedata$deep.Session_tw=="Baseline" | divedata$deep.Session_tw=="Post")
  sum(filter3) # 133
  filter3 <- filter3 & !is.na(divedata$deep.surfn_tw) & !is.na(divedata$deep.logging_tw) & !is.na(divedata$deep.postdur_tw) 
  filter3 <- filter3 & !is.na(divedata$deep.dur) & !is.na(divedata$deep.fluken)
  sum(filter3) # 133

  
  
## Model selection
####################################################


 basemodel <- "deep.surfn_tw ~ deep.lognum_tw + deep.postdur_tw" 
 
 mycovariates <- c("deep.dur", 
                   "deep.fluken",
                   
                   "deep.dur:deep.fluken",
                   "deep.flukerate",
                   
                   "ind.calf:deep.dur",
                   "ind.small:deep.dur",
                   "ind.large:deep.dur",
                   
                   "ind.calf:deep.postdur_tw",
                   "ind.small:deep.postdur_tw",
                   "ind.large:deep.postdur_tw")#,
                   
                   #"ind.calf:deep.fluken",
                   #"ind.small:deep.fluken",
                   #"ind.large:deep.fluken")
 
 myformulas <- c("")
 kk <- 1
 for(j in 1:6) {
   temp <- combn(mycovariates,j)
   for(k in 1:(dim(temp)[2])) {
     myformulas[kk] <- paste(temp[,k],collapse=" + ")
     kk <- kk+1
   }
 }
 
 myformulas <- c(basemodel, paste(basemodel, myformulas, sep=" + "))
 length(myformulas) # 466
 
 #### 
 gamlist <- list()
 lmelist <- list()
 for(m in 1:length(myformulas)) {
   fit0 <- try(gamm(as.formula(myformulas[m]),
                    family=quasipoisson(link = "identity"),
                    data=divedata[filter3,],
                    random=list(ind2=~1),
                    correlation=corAR1()),silent=T)
   
   if(class(fit0)[1]!="try-error") {
     gamlist[[m]] <- summary(fit0$gam)
     lmelist[[m]] <- summary(fit0$lme)
   } else {
     print("try-error")
     gamlist[[m]] <- fit0
     lmelist[[m]] <- fit0}
 }
 
 mtab <- data.frame(ID=1:length(myformulas))
 mtab$model <- myformulas
 mtab$Rsq <- NA
 mtab$AIC <- NA
 mtab$pa <- FALSE
 
 mtab$intercept <- NA
 
 mtab$deep.dur <- NA
 mtab$p.deep.dur <- NA
 
 mtab$deep.fluken <- NA
 mtab$p.deep.fluken <- NA
 
 mtab$deep.flukerate <- NA
 mtab$p.deep.flukerate <- NA
 
 mtab$deep.dur.fluken <- NA
 mtab$p.deep.dur.fluken <- NA
 
 # deep.dur interactions
 mtab$ind.calf.dur <- NA
 mtab$p.ind.calf.dur <- NA
 
 mtab$ind.small.dur <- NA
 mtab$p.ind.small.dur <- NA
 
 mtab$ind.large.dur <- NA
 mtab$p.ind.large.dur <- NA
 
 mtab$p.ind.calf.postdur <- NA
 mtab$p.ind.small.postdur <- NA
 mtab$p.ind.large.postdur <- NA
 
 
 for(m in 1:length(mtab$ID)) {
   
   if(class(gamlist[[m]])[1]!="try-error") {
     
     mtab$Rsq[m] <- gamlist[[m]]$r.sq
     mtab$AIC[m] <- AIC(lmelist[[m]])
     mtab$pa[m] <- all(gamlist[[m]]$pTerms.pv < 0.05)
     
     mtab$intercept[m] <- gamlist[[m]]$p.coeff["(Intercept)"]
     
     mtab$deep.dur[m] <- gamlist[[m]]$p.coeff["deep.dur"]
     mtab$deep.fluken[m] <- gamlist[[m]]$p.coeff["deep.fluken"]
     mtab$deep.flukerate[m] <- gamlist[[m]]$p.coeff["deep.flukerate"]
     
     mtab$p.deep.dur[m] <- gamlist[[m]]$pTerms.pv["deep.dur"]
     mtab$p.deep.fluken[m] <- gamlist[[m]]$pTerms.pv["deep.fluken"]
     mtab$p.deep.flukerate[m] <- gamlist[[m]]$pTerms.pv["deep.flukerate"]
     
     # postdur interactions
     
     # estimates
     if(!is.na(match("deep.dur:deep.fluken", strsplit(myformulas[m]," ")[[1]]))) {mtab$deep.dur.fluken[m] <- as.numeric(na.omit(c(gamlist[[m]]$p.coeff["deep.dur:deep.fluken"], gamlist[[m]]$p.coeff["deep.fluken:deep.dur"])))}
     # p-values
     if(!is.na(match("deep.dur:deep.fluken", strsplit(myformulas[m]," ")[[1]]))) {mtab$p.deep.dur.fluken[m] <- as.numeric(na.omit(c(gamlist[[m]]$pTerms.pv["deep.dur:deep.fluken"], gamlist[[m]]$pTerms.pv["deep.fluken:deep.dur"])))}

     # deep.dur interactions
     
     # estimates
     if(!is.na(match("ind.calf:deep.dur", strsplit(myformulas[m]," ")[[1]]))) {mtab$ind.calf.dur[m] <- as.numeric(na.omit(c(gamlist[[m]]$p.coeff["deep.dur:ind.calf"], gamlist[[m]]$p.coeff["ind.calf:deep.dur"])))}
     if(!is.na(match("ind.small:deep.dur", strsplit(myformulas[m]," ")[[1]]))) {mtab$ind.small.dur[m] <- as.numeric(na.omit(c(gamlist[[m]]$p.coeff["deep.dur:ind.small"], gamlist[[m]]$p.coeff["ind.small:deep.dur"])))}
     if(!is.na(match("ind.large:deep.dur", strsplit(myformulas[m]," ")[[1]]))) {mtab$ind.large.dur[m] <- as.numeric(na.omit(c(gamlist[[m]]$p.coeff["deep.dur:ind.large"], gamlist[[m]]$p.coeff["ind.large:deep.dur"])))}
     # p-values
     if(!is.na(match("ind.calf:deep.dur", strsplit(myformulas[m]," ")[[1]]))) {mtab$p.ind.calf.dur[m] <- as.numeric(na.omit(c(gamlist[[m]]$pTerms.pv["deep.dur:ind.calf"], gamlist[[m]]$pTerms.pv["ind.calf:deep.dur"])))}
     if(!is.na(match("ind.small:deep.dur", strsplit(myformulas[m]," ")[[1]]))) {mtab$p.ind.small.dur[m] <- as.numeric(na.omit(c(gamlist[[m]]$pTerms.pv["deep.dur:ind.small"], gamlist[[m]]$pTerms.pv["ind.small:deep.dur"])))}
     if(!is.na(match("ind.large:deep.dur", strsplit(myformulas[m]," ")[[1]]))) {mtab$p.ind.large.dur[m] <- as.numeric(na.omit(c(gamlist[[m]]$pTerms.pv["deep.dur:ind.large"], gamlist[[m]]$pTerms.pv["ind.large:deep.dur"])))}

     if(!is.na(match("ind.calf:deep.postdur_tw", strsplit(myformulas[m]," ")[[1]]))) {mtab$p.ind.calf.postdur[m] <- as.numeric(na.omit(c(gamlist[[m]]$pTerms.pv["deep.postdur_tw:ind.calf"], gamlist[[m]]$pTerms.pv["ind.calf:deep.postdur_tw"])))}
     if(!is.na(match("ind.small:deep.postdur_tw", strsplit(myformulas[m]," ")[[1]]))) {mtab$p.ind.small.postdur[m] <- as.numeric(na.omit(c(gamlist[[m]]$pTerms.pv["deep.postdur_tw:ind.small"], gamlist[[m]]$pTerms.pv["ind.small:deep.postdur_tw"])))}
     if(!is.na(match("ind.large:deep.postdur_tw", strsplit(myformulas[m]," ")[[1]]))) {mtab$p.ind.large.postdur[m] <- as.numeric(na.omit(c(gamlist[[m]]$pTerms.pv["deep.postdur_tw:ind.large"], gamlist[[m]]$pTerms.pv["ind.large:deep.postdur_tw"])))}
     
   }
 }
 

 ### Model selection with all covariates
 
   mtab$dAIC <- mtab$AIC-min(mtab$AIC)
   mtab <- mtab[order(mtab$AIC),]
   View(mtab[mtab$dAIC<2,])
 
 ### Model selection without fluke rate
 
   tBool <- is.na(mtab$deep.fluken) & is.na(mtab$deep.dur.fluken) & is.na(mtab$deep.flukerate)
   mtab2 <- mtab[tBool,]
   mtab2$dAIC <- mtab2$AIC-min(mtab2$AIC)
   View(mtab2[mtab2$dAIC<2,])
 
 ### Check correlation between dive duration and fluke number

   tBool <-   divedata$deep.first & (divedata$deep.Session_tw=="Baseline" | divedata$deep.Session_tw=="Post")
   plot(divedata$deep.dur[tBool], divedata$deep.fluken[tBool])
   cor(divedata$deep.dur[tBool], divedata$deep.fluken[tBool], method = "spearman")
   #  0.8604748

   
## Interpret best models
####################################################

 fit1 <- gamm(deep.surfn_tw ~ deep.lognum_tw + deep.postdur_tw + deep.fluken,
              family=quasipoisson(link = "identity"),
              data=divedata[filter3,],
              random=list(ind2=~1),
              correlation=corAR1())
 summary(fit1$gam)
 summary(fit1$gam)$pTerms.table
 
 fit2 <- gamm(deep.surfn_tw ~ deep.lognum_tw + deep.postdur_tw + deep.dur,
              family=quasipoisson(link = "identity"),
              data=divedata[filter3,],
              random=list(ind2=~1),
              correlation=corAR1())
 summary(fit2$gam)
 summary(fit2$gam)$pTerms.table
 
  ### Set up data frame for model prediction
   
   preddata0 <- data.frame(deep.fluken=seq(0,320,1))
   preddata0$deep.lognum_tw <- 0 #mean(divedata$deep.lognum_tw[filter3])
   preddata0$deep.postdur_tw <- 0 #mean(divedata$deep.postdur_tw[filter3])# 
   
   preddata1 <- data.frame(deep.dur=seq(1,13,0.1))
   preddata1$deep.lognum_tw <-  0 #mean(divedata$deep.logging_tw[filter3])
   preddata1$deep.postdur_tw <- 0 #mean(divedata$deep.postdur_tw[filter3])#
   preddata1$ind.calf <- 0

  ### Plot predictions
 
     par(mfrow=c(1,3), mar=c(4, 4, 2, 1), lwd=1)
    
      # Dive duration
      
      nbreaths0 <- divedata$deep.surfn_tw
      nbreaths <- nbreaths0+divedata$deep.lognum_tw
      tBool <- filter3
      tBool2 <- tBool & divedata$ind.calf==1
      plot(divedata$deep.dur[tBool],
           nbreaths0[tBool]/divedata$deep.postdur_tw[tBool],
           pch=divedata$pch[tBool], col=NA, ylim=c(1.5, 7.5),
           xlab="Dive duration (min)", ylab="Breathing rate (min-1)")
      segments(x0=divedata$deep.dur[tBool], 
               x1=divedata$deep.dur[tBool],
               y0=nbreaths0[tBool]/divedata$deep.postdur_tw[tBool],
               y1=nbreaths[tBool]/divedata$deep.postdur_tw[tBool],
               col=divedata$ind.col[tBool])
      points(divedata$deep.dur[tBool],
             nbreaths0[tBool]/divedata$deep.postdur_tw[tBool], 
             pch=divedata$pch[tBool], col=divedata$ind.col[tBool])
      points(divedata$deep.dur[tBool2],
             nbreaths0[tBool2]/divedata$deep.postdur_tw[tBool2], 
             pch=3, col="white")
      title("a.")

      preddata <- divedata[filter3,c("deep.lognum_tw","deep.postdur_tw","ind.calf")]
      preddata$deep.dur <- 0
      preddata$ind.calf <- 0
      preds <- predict(fit2$gam, newdata=preddata, se.fit=T, type="response")
      breath_detrend <- preds$fit
      
      tBool <- filter3
      tBool2 <- tBool & divedata$ind.calf==1
      plot(divedata$deep.dur[tBool],
           divedata$deep.surfn_tw[tBool]-breath_detrend,
           pch=divedata$pch[tBool], col=divedata$ind.col[tBool],
           xlab="Dive duration (min)", ylab="Recovery breaths")
      points(divedata$deep.dur[tBool2],
             divedata$deep.surfn_tw[tBool2]-breath_detrend[divedata$ind.calf[tBool]==1], pch=3, col="white")
      
      preds <- predict(fit2$gam, newdata=preddata1, se.fit=T, type="response")
      lines(preddata1$deep.dur, (preds$fit))
      lines(preddata1$deep.dur, (preds$fit+2*preds$se.fit), lty=2)
      lines(preddata1$deep.dur, (preds$fit-2*preds$se.fit), lty=2)
    
      title("b.")

      # Plot as a function # fluke strokes
      
      preddata <- divedata[filter3,c("deep.lognum_tw","deep.postdur_tw")]
      preddata$deep.dur <- 0
      preddata$deep.fluken <- 0
      preddata$ind.calf <- 0
      preds <- predict(fit1$gam, newdata=preddata, se.fit=T, type="response")
      breath_detrend <- preds$fit
      
      tBool <- filter3
      tBool2 <- tBool & divedata$ind.calf==1
      plot(divedata$deep.fluken[tBool],
           divedata$deep.surfn_tw[tBool]-breath_detrend,
           pch=divedata$pch[tBool], col=divedata$ind.col[tBool],
           xlab="# Fluke strokes", ylab="Recovery breaths")
      points(divedata$deep.fluken[tBool2],
             divedata$deep.surfn_tw[tBool2]-breath_detrend[divedata$ind.calf[tBool]==1], pch=3, col="white")
      
      table(cut(divedata$deep.dur,seq(0,14,1)), cut(divedata$deep.fluken,seq(0,350,50)))
      
      temp_cols <- "black" #c("lightgrey", "darkgrey", "black")
      D <- 0#c(2,6,10)
      for(j in 1:length(D)) {
        D1 <- D-2
        D2 <- D+2
        preddata <- preddata0
        preddata$deep.dur <- D[j] #preddata$deep.fluken/FR[j]
        summary(preddata$deep.dur)
        preds <- predict(fit1$gam, newdata=preddata, se.fit=T, type="response")
        y <- (preds$fit)-coef(fit1$gam)[1]
        n <- length(y)
        lines(preddata$deep.fluken, y, col=temp_cols[j])
        lines(preddata$deep.fluken, y-2*preds$se.fit, lty=2, col=temp_cols[j])
        lines(preddata$deep.fluken, y+2*preds$se.fit, lty=2, col=temp_cols[j])
      }
      title("c.")
      
  

  

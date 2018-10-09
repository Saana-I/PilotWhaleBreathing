## R-tools
####################################################

 library(mgcv)
 getCI <- function(b, preddata, n.rep=10000) {
   
   ## extract parameter estiamtes and cov matrix...
   beta <- coef(b)
   Vb <- vcov(b)
   
   ## simulate replicate beta vectors from posterior...
   Cv <- chol(Vb)
   nb <- length(beta)
   br <- t(Cv) %*% matrix(rnorm(n.rep*nb),nb,n.rep) + beta
   
   ## turn these into replicate linear predictors...
   Xp <- predict(b,newdata=preddata,type="lpmatrix")
   lp <- Xp%*%br
   fv <- exp(lp) ## ... finally, replicate expected value vectors
   
   ## now simulate from Gamma deviates with mean as in fv
   ## and estimated scale...
   scale_par <- summary(b)$scale
   yr <- matrix(rgamma(fv*0,shape=1/scale_par,scale=fv*scale_par),nrow(fv),ncol(fv))
   
   ## compute 95% prediction interval...
   PI <- apply(yr,1,quantile,prob=c(.025,0.975))
   
   return(PI)
 }
 
 
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
 
 
## Fit models
####################################################
 
   ### Prediction model
   fit1 <- gamm(dive.rate ~ te(deep.tsince, deep.dur, bs="tp") + s(deep.tnext, bs="tp"), 
                family=Gamma(link="log"),
                data=divedata[filter1,],
                random=list(ind2=~1, deep.index=~1),
                correlation=corAR1())
   
   fit_tsince <- fit1
   save(fit_tsince, file="Step2_Model_1_fit.Rd")
 
   ### Estimate baseline recovered IBR
   preddata <- data.frame(deep.dur=divedata$deep.dur[divedata$deep.filter])
   preddata$deep.tnext <- 30
   preddata$deep.tsince <- 30
   preds0 <- predict(fit1$gam, preddata, type="response", se.fit=T)
   mean(preds0$fit) # 3.398772
   mean(preds0$se.fit) # 0.4670192
   
   ### Compare to other models
   
       # deep.tsince as main effect
       temp <- gamm(dive.rate ~ s(deep.tsince, bs="tp") + s(deep.dur, bs="tp") + s(deep.tnext, bs="tp"), 
                    family=Gamma(link="log"),
                    data=divedata[filter1,],
                    random=list(ind2=~1, deep.index=~1),
                    correlation=corAR1())
       summary(temp$gam)
     
       # include next deep.dur as main effect
       temp <- gamm(dive.rate ~ te(deep.tsince, deep.dur, bs="tp") + s(deep.tnext, bs="tp")
                    + s(deep.ndur, bs="tp"), 
                    family=Gamma(link="log"),
                    data=divedata[filter1,],
                    random=list(ind2=~1, deep.index=~1),
                    correlation=corAR1())
       summary(temp$gam)
    
       # include next deep.dur as tensor effect
       temp <- gamm(dive.rate ~ te(deep.tsince, deep.dur, bs="tp")
                    + te(deep.tnext, deep.ndur, bs="tp"), 
                    family=Gamma(link="log"),
                    data=divedata[filter1,],
                    random=list(ind2=~1, deep.index=~1),
                    correlation=corAR1())
       summary(temp$gam)

   
  ### Check serial correlation 
   
   temp <-  gamm(dive.rate ~ te(deep.tsince, deep.dur, bs="tp") + s(deep.tnext, bs="tp"), 
                  family=Gamma(link="log"),
                  data=divedata[filter1,],
                  random=list(ind2=~1, deep.index=~1))
   
    par(mfrow=c(2,2), mar=c(4, 4, 4, 1), lwd=0.5)
      
    acf(residuals(temp$lme, type="normalized"), main="") # 
    title("No AR1")
    pacf(residuals(temp$lme, type="normalized"), main="") # 
      
    acf(residuals(fit1$lme, type="normalized"), main="") # 
    title("AR1")
    pacf(residuals(fit1$lme, type="normalized"), main="") # 

  

## Model prediction
####################################################    
    
  ylims <- c(0,15)
  par(mfrow=c(2,2), mar=c(4, 4, 2, 1), lwd=0.5)
  
  ### Pre-dive rates
    
      plot(divedata$deep.tnext[filter1], divedata$dive.rate[filter1],
           col=NA, xlab="Time to next dive (min)", ylab="IBR (# breaths min-1)",
           xlim=c(30,0), ylim=ylims) #, yaxt="n"
      title("a.")
      
      tempI <- unique(divedata$deep.index[filter1])
      for(j in tempI) {
        tBool <- divedata$deep.index==j & filter1
        lines(divedata$deep.tnext[tBool], divedata$dive.rate[tBool], col="grey")
      }
      grid(col="darkgrey", lty=2)
      preddata <- data.frame(deep.tnext=seq(0,30,0.1))
      preddata$deep.tsince <-  30 
      preddata$deep.dur <- 2
      preddata$deep.ndur <- 2 
      preds <- predict(fit1$gam, preddata, type="response", se.fit=T)
      lines(preddata$deep.tnext, preds$fit, type="l", lty=1, lwd=1)
      lines(preddata$deep.tnext, preds$fit-2*preds$se.fit, type="l", lty=2, lwd=1)
      lines(preddata$deep.tnext, preds$fit+2*preds$se.fit, type="l", lty=2, lwd=1)

      preds$fit[preddata$deep.tnext==0] # 4.635166 - IBR just before a dive
      preds$se.fit[preddata$deep.tnext==0] # 0.6485845 

  ### Post-dive rates
      
      plot(divedata$deep.tsince[filter1], divedata$dive.rate[filter1],
           col=NA, xlab="Time since dive (min)", ylab="IBR (# breaths min-1)",
           xlim=c(0,30), ylim=ylims) # , yaxt="n", ylim=ylims)
      grid(col="darkgrey", lty=2)
      title("b.")
      
      tempI <- unique(divedata$deep.index[filter1])
      for(j in tempI) {
        tBool <- divedata$deep.index==j & filter1
        lines(divedata$deep.tsince[tBool], divedata$dive.rate[tBool], col="grey")
      }
      preddata <- data.frame(deep.tsince=seq(0,30,0.1))
      preddata$deep.tnext <- 30 
      preddata$deep.fluken <- mean(divedata$deep.fluken[filter1], na.rm=T)
      preddata$deep.dur <- 2 
      preds <- predict(fit1$gam, preddata, type="response", se.fit=T)
      lines(preddata$deep.tsince, preds$fit, type="l", lty=1, lwd=1)
      lines(preddata$deep.tsince, preds$fit-2*preds$se.fit, type="l", lty=2, lwd=1)
      lines(preddata$deep.tsince, preds$fit+2*preds$se.fit, type="l", lty=2, lwd=1)
      
      preds$fit[preddata$deep.tsince==0] # 3.42951 - IBR following a 2-min dive
      preds$se.fit[preddata$deep.tsince==0] # 0.3169234 
 
      preddata$deep.dur <- 10 # mean(divedata$deep.dur[filter1], na.rm=T)
      preds <- predict(fit1$gam, preddata, type="response", se.fit=T)
      lines(preddata$deep.tsince, preds$fit, type="l", lty=1, lwd=1, col="red")
      lines(preddata$deep.tsince, preds$fit-2*preds$se.fit, type="l", lty=2, lwd=1, col="red")
      lines(preddata$deep.tsince, preds$fit+2*preds$se.fit, type="l", lty=2, lwd=1, col="red")

      preds$fit[preddata$deep.tsince==0] # 7.678576  - IBR following a 10-min dive
      preds$se.fit[preddata$deep.tsince==0] # 0.606905 
      
  ### Initial post-dive rates
    
      tBool <- filter1 & divedata$deep.first
      tBool2 <- tBool & divedata$ind.calf==1
      plot(divedata$deep.dur[tBool], 
           divedata$dive.rate[tBool], xlab="Dive duration (min)", 
           ylab="Initial IBR (# breaths min-1)", col=divedata$ind.col[tBool], pch=divedata$pch[tBool])
      grid(col="darkgrey")
      points(divedata$deep.dur[tBool2], divedata$dive.rate[tBool2], pch=3, col="white")
      title("c.")
      
      beta1 <- 0
      beta1.l <- 0
      beta1.h <- 0
      tr <- 0
      tr2 <- 0
      tr3 <- 0
      TH <- 5 #0.08
      TH2 <- 4 #0.07
      TH3 <- 3.4 #0.06
      durvec <- seq(1,14,0.5)
      
      for(j in 1:length(durvec)) {
        
        preddata <- data.frame(deep.tsince=seq(0,30,0.1))
        preddata$deep.tnext <- 30 # mean(divedata$deep.fluken[filter1], na.rm=T)
        preddata$deep.fluken <- mean(divedata$deep.fluken[filter1], na.rm=T)
        preddata$deep.dur <- durvec[j]
        preds <- predict(fit1$gam, preddata, type="response", se.fit=T)
        beta1[j] <- preds$fit[1]
        beta1.l[j] <- preds$fit[1]-2*preds$se.fit[1]
        beta1.h[j] <- preds$fit[1]+2*preds$se.fit[1]
        tr[j] <- preddata$deep.tsince[min(which(preds$fit<TH))]
        tr2[j] <- preddata$deep.tsince[min(which(preds$fit<TH2))]
        tr3[j] <- preddata$deep.tsince[min(which(preds$fit<TH3))]
      }
      lines(durvec, beta1,lty=1)
      lines(durvec, beta1.l,lty=2)
      lines(durvec, beta1.h,lty=2)
 
      mean(diff(beta1)/diff(durvec)) # 0.523014
      
      # recovery time after max duration dive
      preddata <- data.frame(deep.tsince=seq(0,30,0.1))
      preddata$deep.tnext <- 30 # mean(divedata$deep.fluken[filter1], na.rm=T)
      preddata$deep.dur <- max(divedata$deep.dur, na.rm=T)
      preds <- predict(fit1$gam, preddata, type="response", se.fit=T)
      preddata$deep.tsince[min(which(preds$fit<TH3))]
      # 9.5
      
  ### Predicted recovery time
    
    plot(divedata$deep.dur[tBool], divedata$deep.postdur[tBool], 
         pch=divedata$pch[tBool],col=divedata$ind.col[tBool],
         xlab="Dive duration (min)", 
         ylab="Recovery time (min)", 
         ylim=c(0,15))
    grid(col="darkgrey")
    points(divedata$deep.dur[tBool2], divedata$deep.postdur[tBool2], 
           pch=3, col="white")
    lines(durvec, tr3, lty=1)
    text(durvec[durvec==7], tr3[durvec==7]-0.8, paste("IBR <",TH3), cex=0.8)
    title("d.")
    
    mean(diff(tr3[durvec>=2])/diff(durvec[durvec>=2])) # 0.775
    
    # time to recovery following a 2, 10 and 14-min dives
    tr3[durvec==2] # 0.2
    tr3[durvec==10] # 7.3
    tr3[durvec==14] # 9.5
    
  

  ### Predict post-dive rate following different duration of dives
    
    par(mfrow=c(3,2), mar=c(4, 4.5, 2, 0.5), lwd=0.5)

    for(j in c(3,6,8,10,12)) {
      
      tBool <- filter1 & divedata$deep.dur>((j-1)) & divedata$deep.dur<((j+1))
      plot(divedata$deep.tsince[tBool], divedata$dive.rate[tBool],
           col="grey", xlab="Time since dive (min)", ylab=expression(paste("IBR (min"^"-1", ")")),
           xlim=c(0,15), ylim=c(0,15))
      grid(col="darkblue")
      
      preddata <- data.frame(deep.tsince=seq(0,30,0.1))
      preddata$deep.tnext <- 30 # mean(divedata$deep.fluken[filter1], na.rm=T)
      preddata$deep.fluken <- mean(divedata$deep.fluken[filter1], na.rm=T)
      preddata$deep.dur <- j
      preds <- predict(fit1$gam, preddata, type="response", se.fit=T)
      lines(preddata$deep.tsince, preds$fit, type="l", lty=1, lwd=1)
      title(paste(j-1,"-",j+1, " min dive", sep=""))
    }
    
    plot(divedata$deep.tsince[filter1], divedata$dive.rate[filter1],
         col="grey", xlab="Time since dive (min)", ylab=expression(paste("IBR (min"^"-1", ")")),
         xlim=c(0,15), ylim=c(0,15))
    
    for(j in 2:14) {
      
      preddata <- data.frame(deep.tsince=seq(0,30,0.1))
      preddata$deep.tnext <- 30 # mean(divedata$deep.fluken[filter1], na.rm=T)
      preddata$deep.fluken <- mean(divedata$deep.fluken[filter1], na.rm=T)
      preddata$deep.dur <- j
      preds <- predict(fit1$gam, preddata, type="response", se.fit=T)
      lines(preddata$deep.tsince, preds$fit, type="l", lty=1, lwd=1)
    }
    title("2-14min")


## R-tools 
####################################################

  library(mgcv)
  
  ### Functions from Fahlman et al 2016 - based on measurements on bottlenose dolphins

  # function to estimate TLC from TBM
  TLC_fun <- function(TBM) {return(0.135*TBM^0.92)}

  # function to estimate tidal volume as a proportion of TLC, given breath number
  # assuming TLC_fun for TLC relationship with body mass
  b_VTp_fun <- function(n) {return((33.7+21.1*exp(-0.111*n))/100)}
  
  # function to estimate O2 fraction (expired-inspired), given a breath number
  # should be independent of mass
  b_DO2_fun <- function(n) {return(0.2094-(14.3+2.81*(1-exp(-0.0495*n)))/100)}  
  
  # Calculate VTp (% TLC), DO2 (O2 fraction), VO2 (per breath)
  # with 95% measurement error
  br_VO2_fun <- function(n, its, TBM) {
    
    TLC <- TLC_fun(TBM)
    
    VT_a <- rnorm(its, mean=33.7, sd=1.18)
    VT_b <- rnorm(its, mean=21.1, sd=1.7)
    VT_c <- rnorm(its, mean=-0.111, sd=0.022)
    
    D_a <- rnorm(its, mean=14.3, sd=0.2)
    D_b <- rnorm(its, mean=2.81, sd=0.43)
    D_c <- rnorm(its, mean=-0.0495, sd=0.0176)
    
    VTp.l <- NA
    VTp.u <- NA
    DO2.l <- NA
    DO2.u <- NA
    VO2.l <- NA
    VO2.u <- NA
    
    for(j in 1:length(n)) {
      
      
      VTp <- (VT_a+VT_b*exp(VT_c*n[j]))/100
      DO2 <- 0.2094-(D_a+D_b*(1-exp(D_c*n[j])))/100  
      VTp.l[j] <- quantile(VTp, 0.025)
      VTp.u[j] <- quantile(VTp, 0.975)
      DO2.l[j] <- quantile(DO2, 0.025)
      DO2.u[j] <- quantile(DO2, 0.975)
      VO2.l[j] <- quantile(TLC*VTp*DO2,0.025)
      VO2.u[j] <- quantile(TLC*VTp*DO2,0.975)
      
    }
    
    return(list(VTp.l=VTp.l,
                VTp.u=VTp.u,
                DO2.l=DO2.l,
                DO2.u=DO2.u,
                VO2.l=VO2.l,
                VO2.u=VO2.u))
  }
  
## Load data and Model 1 fit (objct 'fit_tsince')
####################################################

  load("Step1_dive_recovery_data.Rd")
  divedata <- divedata2
  remove(divedata2)
  
  load("Step2_Model_1_fit.Rd")
  
## Set constants for simulation
####################################################
  
  TBM <- 1000 # body mass, kg
  its <- 10000  # number of iterations 
  nmax <- 100 # maximum number of breaths calculated post-dive
  maxT <- 15 # max recovery time allowed
  deep.dur <- 10 # duration of dive
  
  
## Simulations
####################################################
  
    ### Make the same calculations for different values for total lung capacity (TLC)
  
    TLC1 <- 40 
    VTp.l <- b_VTp_fun(1:nmax) # tidal volume (% TLC)
    DO2.l <- b_DO2_fun(1:nmax) # o2 extracted (%)
    VO2.l <- TLC1*VTp.l*DO2.l  # O2 update (l per breath)

    TLC3 <- 100
    VTp.u <- b_VTp_fun(1:nmax) 
    DO2.u <- b_DO2_fun(1:nmax) 
    VO2.u <- TLC3*VTp.u*DO2.u  
  
    TLC_fun <- function(TBM) {return(0.135*TBM^0.92)} # TLC based on total body mass (TBM)
    TLC2 <- TLC_fun(TBM)
    VTp.mu <- b_VTp_fun(1:nmax)
    DO2.mu <- b_DO2_fun(1:nmax) 
    VO2.mu <- TLC2*VTp.mu*DO2.mu
  
    temp <- br_VO2_fun(n=1:nmax, its=its, TBM=TBM) # function to retrieve values with 95% CI, uses TLC_fun
    VTp.l <- temp$VTp.l
    VTp.u <- temp$VTp.u
    DO2.l <- temp$DO2.l
    DO2.u <- temp$DO2.u
  
    ### Predict from Model 1 (fit_tsince)
    
    # Predict as a function of time
    preddata <- data.frame(deep.tsince=seq(0,maxT,0.1))
    preddata$deep.tnext <-  30 # mean(divedata$deep.fluken[filter1], na.rm=T)
    preddata$deep.dur <- deep.dur
    preds <- predict(fit_tsince$gam, preddata, type="response", se.fit=T)
    minval <- min(preds$fit-2*preds$se.fit)
    maxval <- max(preds$fit+2*preds$se.fit)

    # Predict values breath-by-breath
    deep.tsince <- 0
    n <- 1  
    ibr <- preds$fit[1]
    ibr.se <- preds$se.fit[1]
    while(deep.tsince[n] <= maxT) {
      n <- n+1 # add one breath
      deep.tsince[n] <- deep.tsince[n-1]+(1/ibr[n-1])
      preddata <- data.frame(deep.tsince=deep.tsince[n])
      preddata$deep.tnext <-  30 
      preddata$deep.dur <- deep.dur
      preds <- predict(fit_tsince$gam, preddata, type="response", se.fit=T)
      ibr[n] <- preds$fit
      ibr.se[n] <- preds$se.fit
    }
  
  
## Plot predictions
####################################################
  
  par(mfrow=c(2,2), mar=c(4, 4, 1, 3), lwd=1)

  mylabsize <- 1
  mylabpos <- 2.3
  
  # o2 difference
  plot(1:n, ibr, type="l", lty=1, col=1,
       ylim=c(minval, maxval),  yaxt="n", 
       ylab="", xlab="")
  axis(side=4)
  mtext(side=2,expression(paste("O"[2]," extracted (%)")), #expression(paste("O"^"2","extracted (%)"))
        line=mylabpos,cex=mylabsize)
  mtext(side=1,"Breath number",line=mylabpos,cex=mylabsize)
  segments(1:n, 1:n, y0=ibr-2*ibr.se, y1=ibr+2*ibr.se, col="grey")

  par(new=T)
  plot(1:n, DO2.mu[1:n]*100, xaxt="n",xlab="",ylab="", 
       xlim=c(1,n), type="l", col="blue", ylim=100*c(min(DO2.l),max(DO2.u)))
  lines(1:n, DO2.l[1:n]*100, lty=2, col="blue")
  lines(1:n, DO2.u[1:n]*100, lty=2, col="blue")
  axis(side=2, col="blue")
  title("a.")
  
  # VT
  plot(1:n, ibr, type="l", lty=1,
       ylim=c(minval, maxval), yaxt="n", 
       ylab="", xlab="")
  segments(1:n, 1:n, y0=ibr-qnorm(0.975)*ibr.se, y1=ibr+qnorm(0.975)*ibr.se, col="grey")
  axis(side=4)
  mtext(side=2,"Tidal volume (% TLC)",line=mylabpos,cex=mylabsize)
  mtext(side=1,"Breath number",line=mylabpos,cex=mylabsize)
  
  par(new=T)
  plot(1:n, VTp.mu[1:n]*100, xaxt="n", xlab="",ylab="", 
       xlim=c(1,n), type="l", col="blue", ylim=c(min(VTp.l),max(VTp.u))*100)
  lines(1:n, VTp.l[1:n]*100, lty=2, col="blue")
  lines(1:n, VTp.u[1:n]*100, lty=2, col="blue")
  axis(side=2, col="blue")
  title("b.")
  
  # o2 uptake
  plot(deep.tsince, ibr, type="l", lty=1,
       ylim=c(minval, maxval), yaxt="n", 
       ylab="", xlab="")
  segments(deep.tsince, deep.tsince, y0=ibr-qnorm(0.975)*ibr.se, y1=ibr+qnorm(0.975)*ibr.se, col="grey")
  axis(side=4)
  mtext(side=2, expression(paste("O"[2]," uptake (l per breath)")),line=mylabpos,cex=mylabsize)
  mtext(side=1,"Time since dive (min)",line=mylabpos,cex=mylabsize)
  
  par(new=T)
  plot(deep.tsince, VO2.mu[1:n], xaxt="n", xlab="",ylab="", 
       type="l", col="blue", ylim=c(min(VO2.l),max(VO2.u)))
  lines(deep.tsince, VO2.l[1:n], lty=1, col="blue")
  lines(deep.tsince, VO2.u[1:n], lty=1, col="blue")
  axis(side=2, col="blue")
  text(deep.tsince[n]-4, VO2.l[n], expression(paste("TLC ", 40, " ml kg"^"-1","")))
  text(deep.tsince[n]-4, VO2.mu[n], expression(paste("TLC ", 77.7, " ml kg"^"-1","")))
  text(deep.tsince[n]-4, VO2.u[n], expression(paste("TLC ", 100, " ml kg"^"-1","")))
  title("c.")
  
  # o2 uptake
  plot(deep.tsince, ibr, type="l", lty=1,
       ylim=c(minval, maxval),  yaxt="n", 
       ylab="", xlab="")
  segments(deep.tsince, deep.tsince, y0=ibr-qnorm(0.975)*ibr.se, y1=ibr+qnorm(0.975)*ibr.se, col="grey")
  axis(side=4)
  mtext(side=2,expression(paste("Total O2 (ml O2 kg"^"-1",")")),line=mylabpos,cex=mylabsize) # 
  mtext(side=1,"Time since dive (min)",line=mylabpos,cex=mylabsize)
  
  par(new=T)
  plot(deep.tsince, cumsum(VO2.mu)[1:n], xaxt="n", yaxt="n", xlab="",ylab="", 
       type="l", col="blue", ylim=c(min(VO2.l),max(cumsum(VO2.u)[1:n])))
  lines(deep.tsince, cumsum(VO2.l)[1:n], lty=1, col="blue")
  lines(deep.tsince, cumsum(VO2.u)[1:n], lty=1, col="blue")
  axis(side=2, col="blue")
  title("d.")
  text(deep.tsince[n]-4, cumsum(VO2.l[1:n])[n], expression(paste("TLC ", 40, " ml kg"^"-1","")))
  text(deep.tsince[n]-4, cumsum(VO2.mu[1:n])[n], expression(paste("TLC ", 77.7, " ml kg"^"-1","")))
  text(deep.tsince[n]-4, cumsum(VO2.u[1:n])[n], expression(paste("TLC ", 100, " ml kg"^"-1","")))


  
  
## Interpretation 
####################################################
  
  # Total O2 uptake after 10 min
  tempI <- min(which(deep.tsince>10))
  c(cumsum(VO2.l)[tempI],cumsum(VO2.mu)[tempI],cumsum(VO2.u)[tempI]) 
  # 35.10113 68.17024 87.75282

  # Oxygen consumption rate over the dive cycle
  cumsum(VO2.l)[tempI]/(deep.tsince[tempI]+deep.dur) # 1.729076
  cumsum(VO2.mu)[tempI]/(deep.tsince[tempI]+deep.dur) # 3.358055
  cumsum(VO2.u)[tempI]/(deep.tsince[tempI]+deep.dur) # 4.322689
  
  # Mean uptake per breath
  c(mean(VO2.l[1:tempI]),mean(VO2.mu[1:tempI]),mean(VO2.u[1:tempI]))
  # 0.7468325 1.4504307 1.8670813
  
  # Uptake after 10 min
  c(VO2.l[tempI],VO2.mu[tempI],VO2.u[tempI])
  # 0.5551464 1.0781554 1.3878661
  
  # Fixed uptake at the end of the interval & post-dive breathing rate (3.4)
  (VO2.l[tempI]*3.4*10)/20 # 0.943749
  (VO2.mu[tempI]*3.4*10)/20 # 1.832864 
  (VO2.u[tempI]*3.4*10)/20 # 2.359372
  
  ((VO2.l[tempI]*3.4*10)/20)/(cumsum(VO2.l)[tempI]/(deep.tsince[tempI]+deep.dur)) # 0.5458112
  ((VO2.mu[tempI]*3.4*10)/20)/(cumsum(VO2.mu)[tempI]/(deep.tsince[tempI]+deep.dur)) # 0.5458112
  ((VO2.u[tempI]*3.4*10)/20)/(cumsum(VO2.u)[tempI]/(deep.tsince[tempI]+deep.dur)) # 0.5458112
  
  (1-0.5458112)*100 # 45.41888 lower estimates
  
  # Daily energy expenditure assuming 
  # 3 ml O2 kg-1 min-1 oxygen consumption 
  # 20.1 J per ml O2, and  4.1868 joules per calorie. 
  3*60*24 # 4320 ml o2 kg-1  per day
  3*60*24*20.1/1000 # 86.8 kJ per day
  3*60*24*20.1/4.2/1000 # 20.7 kCal per day

  # Average per-breath uptake over the 10-min interval
  c(mean(VO2.l[1:tempI]),mean(VO2.mu[1:tempI]),mean(VO2.u[1:tempI]))
  # 0.7468325 1.4504307 1.8670813
  
  # Diving metabolic rate, assuming 1.7 breaths per min net diving cost
  c(mean(VO2.l[1:tempI]),mean(VO2.mu[1:tempI]),mean(VO2.u[1:tempI]))*1.7
  # 1.269615 2.465732 3.174038
  
  # Stroking cost (ml O2 kg-1)
  c(mean(VO2.l[1:tempI]),mean(VO2.mu[1:tempI]),mean(VO2.u[1:tempI]))*0.086
  # 0.0642276 0.1247370 0.1605690
  
  # Stroking cost (Joules)
  c(mean(VO2.l[1:tempI]),mean(VO2.mu[1:tempI]),mean(VO2.u[1:tempI]))*0.086*20.1
  # 1.290975 2.507215 3.227437
  
  # Stroking cost (ml o2 kg-1 per dive with fluke rate of 18.3)
  c(mean(VO2.l[1:tempI]),mean(VO2.mu[1:tempI]),mean(VO2.u[1:tempI]))*0.086*18.3
  # 1.175365 2.282688 2.938412
  
  # Stroking cost (as proportion of net diving cost)
  (c(mean(VO2.l[1:tempI]),mean(VO2.mu[1:tempI]),mean(VO2.u[1:tempI]))*0.086*18.3)/(
    c(mean(VO2.l[1:tempI]),mean(VO2.mu[1:tempI]),mean(VO2.u[1:tempI]))*1.7)
  #  0.9257647 0.9257647 0.9257647

  


  

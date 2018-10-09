## R-tools
####################################################

 library(mgcv)

## Data
####################################################

  # load dive-by-dive data
  divedata <- read.csv("divedata.csv")
  indsummary <- read.csv("indsummary.csv")
  etab <- read.csv("etab.csv")
  
  divedata$ind <- as.character(divedata$ind)
  divedata$ind.type <- as.character(divedata$ind.type)
  divedata$ind.col <- as.character(divedata$ind.col)
  divedata$Session <- as.character(divedata$Session)
  divedata$dive.GMTtime <- as.POSIXct(as.character(divedata$dive.GMTtime), tz="GMT")
  
  # load dive profile data
  load("diveprofiles.Rd")
  whales <- unique(divedata$ind)
  
  divedata$ind.calf <- as.numeric(divedata$ind.calf)
  divedata$ind.small <- as.numeric(divedata$ind.small)
  divedata$ind.large <- as.numeric(divedata$ind.large)

## Set up dive thresholds
####################################################

  durTH <- seq(from=1/6, to=13, by=1/6)*60 # 10s res.
  depthTH <- seq(from=1, to=600, by=1)     # 1m res
  flukeTH <- seq(from=1, to=320, by=1)     # 1stroke res
  
  myq <- seq(0, 1, by=0.05)
  
  ### Duration
  
    delta <- matrix(NA, length(durTH), length(myq))
    for(j in 1:length(durTH)) {
      
      fBool <- divedata$dive.dur > durTH[j]
      ftab <- divedata[fBool,c("ind", "dive.sfromtot", "dive.end")]
      n <- dim(ftab)[1]
      ftab$next.st <- NA
      ftab$next.st[-n] <- ftab$dive.sfromtot[-1]
      ftab$next.ind <- NA
      ftab$next.ind[-n] <- ftab$ind[-1]
      ftab$diffs <- ftab$next.st-ftab$dive.end
      ftab$diffs[ftab$ind!=ftab$next.ind] <- NA
      
      # calculate median, lower and upper
      delta[j,] <- quantile(ftab$diffs, myq, na.rm=T)
    }
    delta_dur <- delta
  
    
  ### Depth
    
    delta <- matrix(NA, length(depthTH), length(myq))
    delta_depth2 <- delta
    names(quantile(divedata$dive.depth))
    
    for(j in 1:length(depthTH)) {
      
      fBool <- divedata$dive.depth > depthTH[j]
      ftab <- divedata[fBool,c("ind", "ind.type", "dive.sfromtot", "dive.end")]
      n <- dim(ftab)[1]
      ftab$next.st <- NA
      ftab$next.st[-n] <- ftab$dive.sfromtot[-1]
      ftab$next.ind <- NA
      ftab$next.ind[-n] <- ftab$ind[-1]
      ftab$diffs <- ftab$next.st-ftab$dive.end
      ftab$diffs[ftab$ind!=ftab$next.ind] <- NA
      
      # calculate median, lower and upper
      delta[j,] <- quantile(ftab$diffs, myq, na.rm=T)
      for(k in 1:length(myq)) { # calculate individual-average quantile
        delta_depth2[j,k] <- median(tapply(ftab$diffs, ftab$ind, quantile, myq[k], na.rm=T), na.rm=T)
      }
    }
    delta_depth <- delta

    
  ### Select thresholds

    divedata$logging <- divedata$surf.dur > 5
    medianIBI <- 1/median(1/divedata$dive.dur)    
    postdurTH <- 2*medianIBI
    depth0 <- depthTH[min(which(delta_depth[,1]>=postdurTH))] 
    dur0 <- durTH[min(which(delta_dur[,1]>=postdurTH))] # 210s = 3.5min
    
    print(depth0) # 31  m
    print(dur0) # 210 s
    print(dur0/60) # 3.5 min
    
  ### Plot - Figure 1 
    
    par(mfrow=c(1,3), mar=c(4, 4, 2, 2))
    
    tBool <- (divedata$Session=="Baseline" | divedata$Session=="Post")
    tBool2 <- tBool & divedata$ind.calf==1
    
    tempcol <- rep("black", length(tBool))
    tempcol[divedata$ind.calf==1] <- "darkgrey"
    
    
    plot(divedata$dive.depth[tBool], divedata$dive.dur[tBool]/60, col=NA,
         xlab="Dive depth (m)", ylab="Dive duration (min)")
    grid(col="grey")
    points(divedata$dive.depth[tBool], divedata$dive.dur[tBool]/60, cex=1, 
           pch=divedata$pch[tBool], col=tempcol[tBool])
    #points(divedata$dive.depth[tBool2], divedata$dive.dur[tBool2]/60, cex=0.7, col="white",pch=3)
    abline(v=depth0)
    title("a.")
    
    plot(depthTH, delta_depth[,1]/60, col=NA,
         ylim=c(0,30), xlim=c(0,550),
         xlab="Dive depth threshold (m)", 
         ylab="Post-dive surface interval (min)")
    grid(col="darkgrey")
    title("b.")
    
    lines(depthTH, delta_depth[,myq==0]/60, lty=2) # minimum
    text(150,1+delta_depth[depthTH==150,myq==0]/60,"Min")
    
    lines(depthTH, delta_depth[,myq==0.05]/60, lty=2) 
    text(150,1+delta_depth[depthTH==150,myq==0.05]/60,"5%")
    
    lines(depthTH, delta_depth[,myq==0.25]/60, lty=2)
    text(150,1+delta_depth[depthTH==150,myq==0.25]/60,"25%")
    
    lines(depthTH, delta_depth[,myq==0.5]/60, lty=1)
    text(150,1+delta_depth[depthTH==150,myq==0.5]/60,"50%")
    
    rug(divedata$dive.depth)
    abline(v=depth0)
    
    
    plot(durTH/60, delta_dur[,1]/60, col=NA,
         ylim=c(0,30), xlim=c(0,11), xlab="Dive duration threshold (min)", 
         ylab="Post-dive surface interval (min)")
    grid(col="darkgrey")
    lines(durTH/60, delta_dur[,myq==0]/60, lty=2)
    text(4.5,-0.6+delta_dur[durTH/60==5,myq==0]/60,"Min")
    title("c.")
    
    lines(durTH/60, delta_dur[,myq==0.05]/60, lty=2)
    text(4.5,0.5+delta_dur[durTH/60==5,myq==0.05]/60,"5%")
    
    lines(durTH/60, delta_dur[,myq==0.25]/60, lty=2)
    text(4.5,1+delta_dur[durTH/60==5,myq==0.25]/60,"25%")
    
    lines(durTH/60, delta_dur[,myq==0.5]/60, lty=1)
    text(4.5,1.1+delta_dur[durTH/60==5,myq==0.5]/60,"50%")
    
    rug(divedata$dive.dur/60)


## Make dive cycle data
####################################################

  # deep.filter: filter for dive cycles that include deep dive 
  divedata$deep.filter <- divedata$dive.depth>=depth0
  divedata$deep.post <- FALSE   # whether datapoint is within post-interval
  divedata$deep.first <- FALSE  # whether first breath after dive
  
  # deep.index: unique index for each dive cycle (dive + post dives)
  divedata$deep.index <- NA
  divedata$deep.index[divedata$deep.filter] <- 1:sum(divedata$deep.filter)
  
  divedata$deep.tsince <- NA  # time since end of dive (NA when unknown)
  divedata$deep.tnext <- NA   # time to next dive (NA when unknown)
  
  divedata$deep.postdur <- NA # post dur: duration of interval (also for intervals when tnext is unknown)
  divedata$deep.surfmax <- NA # maximum surface interval duration
  divedata$deep.logging <- NA # time spent logging (surf.dur > medianIBI/2)
  divedata$deep.surfn <- NA   # total number of breaths during PDSI
  divedata$deep.surfIBI <- NA   # median IBI during surface
  
  # 5min into post
  divedata$deep.surfn5 <- NA
  divedata$deep.logging5 <- NA # time spent logging (surf.dur > medianIBI/2)
  
  divedata$deep.depth <- NA   # dive depth
  divedata$deep.dur <- NA     # dive duration
  divedata$deep.fluken <- NA  # dive fluke rate
  
  divedata$deep.ndepth <- NA # next dive depth
  divedata$deep.ndur <- NA # next dur
  
  divedata$deep.Session <- NA
  divedata$deep.Session5 <- NA
  divedata$deep.Session30 <- NA
  divedata$deep.last_in_record <- FALSE
  
  tempI <- which(divedata$deep.filter)
  
  for(j in 1:length(tempI)) {
    
    cwhale <- divedata$ind[tempI[j]]
    cwhaleET <- max(divedata$dive.end[divedata$ind==cwhale])
    temp <- tempList[[which(whales==cwhale)]]
    
    ST0 <- divedata$dive.sfromtot[tempI[j]]
    ST <- divedata$dive.end[tempI[j]]
    
    # time to next dive
    remove(ET)
    timetonext <- FALSE
    isLast <- FALSE
    if(tempI[j] == tempI[length(tempI)]){
      ET <- cwhaleET # end time for the very last dive cycle
    } else {
      if(divedata$ind[tempI[j]]==divedata$ind[tempI[j+1]]) { # when there are more dives in the tag record,
        ET <- divedata$dive.sfromtot[tempI[j+1]] # dive cycle ends to the start of the next deep dive
        timetonext <- TRUE # indicator that there is a dive following this one
        ETdur <- divedata$dive.depth[tempI[j+1]]
        ETdepth <- divedata$dive.dur[tempI[j+1]]
        } else {
        ET <- cwhaleET # end time for the last dive cycle for each tag record
        ETdur <- NA
        ETdepth <- NA
        isLast <- TRUE
        print(cwhale)
      }
    }
    
    # Filter for dive cycle
    tBool <- divedata$ind==divedata$ind[tempI[j]] & divedata$dive.sfromtot>=ST0 & divedata$dive.sfromtot<ET 
    # Filter for post-dive interval
    pBool <- divedata$ind==divedata$ind[tempI[j]] & divedata$dive.sfromtot>=ST & divedata$dive.sfromtot<ET 
  
    divedata$deep.index[tBool] <- divedata$deep.index[tempI[j]] # unique index for each dive cycle
    divedata$deep.post[pBool] <- TRUE # data rows to be included in GAMM
    divedata$deep.first[pBool][1] <- TRUE
    divedata$deep.last_in_record[tBool] <- isLast
    
    divedata$deep.tsince[tBool] <- divedata$dive.sfromtot[tBool]-ST # time since the ST of post-dive
    divedata$deep.postdur[tBool] <- ET-ST # duration of the post-dive interval

    # Filter for post-dive interval <5min 
    pBool5 <- pBool & divedata$dive.sfromtot < (ST+5*60)
    # Filter for post-dive interval <30min 
    pBool30 <- pBool & divedata$dive.sfromtot < (ST+30*60)
    
    if(timetonext) { # if the post-dive interval ends in another dive, record time to the next dive
      divedata$deep.tnext[tBool] <- ET-divedata$dive.sfromtot[tBool]
      divedata$deep.ndepth[tBool] <- ETdur # next dive depth
      divedata$deep.ndur[tBool] <- ETdepth # next dur
      
    }
    
    divedata$deep.surfn[tBool] <- sum(pBool)
    divedata$deep.surfn5[tBool] <- sum(pBool5)
    
    divedata$deep.logging[tBool] <- sum(divedata$surf.dur*pBool*divedata$logging)
    divedata$deep.logging5[tBool] <- sum(divedata$surf.dur*pBool5*divedata$logging) 
    
    # How much logging period exceeds end of 5 min interval?
    if(sum(divedata$surf.dur*pBool5*divedata$logging)>0) {
      # end time of last breath/surface interval
      ET_temp <- max(divedata$dive.end[pBool5 & divedata$logging]+divedata$surf.dur[pBool5 & divedata$logging])
      if(ET_temp > (ST+5*60)) {
        print(paste("correcting logging for # ", tempI[j]))
        divedata$deep.logging5[tBool] <- sum(divedata$surf.dur*pBool5*divedata$logging)-(ET_temp-(ST+5*60))
      }
    }
    
    if(sum(pBool)>0) {
    divedata$deep.surfmax[tBool] <- max(divedata$surf.dur[pBool])   # maximum surface duration
    divedata$deep.surfIBI[tBool] <- median(1/(divedata$dive.dur[pBool]))
    }
    
    divedata$deep.depth[tBool] <- divedata$dive.depth[tempI[j]]
    divedata$deep.dur[tBool] <- divedata$dive.dur[tempI[j]]
    divedata$deep.fluken[tBool] <- divedata$dive.fluken[tempI[j]]
    
    sessions <- unique(divedata$Session[tBool])
    exposures <- sessions[sessions!="Baseline" & sessions!="Post"]
    if(length(exposures)>0) {divedata$deep.Session[tBool] <- exposures[1]
    } else {divedata$deep.Session[tBool] <- sessions[1]}
    
    sessions <- unique(c(divedata$Session[tempI[j]],divedata$Session[pBool5]))
    exposures <- sessions[sessions!="Baseline" & sessions!="Post"]
    if(length(exposures)>0) {divedata$deep.Session5[tBool] <- exposures[1]
    } else {divedata$deep.Session5[tBool] <- sessions[1]}
    
    sessions <- unique(c(divedata$Session[tempI[j]],divedata$Session[pBool30]))
    exposures <- sessions[sessions!="Baseline" & sessions!="Post"]
    if(length(exposures)>0) {divedata$deep.Session30[tBool] <- exposures[1]
    } else {divedata$deep.Session30[tBool] <- sessions[1]}
    
  }
  

  divedata$deep.postBR <- divedata$deep.surfn/(divedata$deep.postdur/60) 
  divedata$deep.postBR5 <- divedata$deep.surfn5/pmin(5, divedata$deep.postdur/60)
  
  # inverse-transform dive & calculate fluke rates
  divedata$dive.rate <- 1/(divedata$dive.dur/60)
  divedata$dive.flukerate <- divedata$dive.fluken/(divedata$dive.dur/60)
  divedata$deep.flukerate <- divedata$deep.fluken/(divedata$deep.dur/60)
  

  
## Data filtering & save
####################################################

    # 1. Loading & recovery ibr time series models
    postwindow1 <- 30
    filter1 <- divedata$deep.post & divedata$deep.postdur>40
    filter1 <- filter1 & divedata$deep.tsince<(postwindow1*60) & divedata$deep.tnext<(postwindow1*60) 
    filter1 <- filter1 & (divedata$deep.Session30=="Baseline" | divedata$deep.Session30=="Post")
 
    # 2. Dive-by-dive analysis for # of breaths (assume 30 min window)
    filter2 <- divedata$deep.first & divedata$deep.postdur>40
    filter2 <- filter2 & (divedata$deep.Session30=="Baseline" | divedata$deep.Session30=="Post")

   
   keeps <-   c("ind","dive.sfromtot","Session",
                "dive.dur","dive.rate", "dive.depth",
                "dive.end", "surf.dur", "logging", 
                "dive.fluken","dive.flukerate",
                "ind.type","ind.calf","ind.small","ind.med", "ind.large","ind.TBM",
                
                "pch","ind.col",
                
                "deep.filter","deep.post","deep.first","deep.index",
                "deep.tsince","deep.tnext",
                "deep.postdur","deep.surfn","deep.postBR","deep.logging","deep.surfIBI",
                "deep.surfn5","deep.postBR5","deep.logging5", 
                
                "deep.depth","deep.dur", "deep.ndepth", "deep.ndur",
                "deep.fluken","deep.flukerate", 
                "deep.Session")
   
   divedata2 <- divedata[,keeps]
  
   divedata2$deep.tsince <- divedata$deep.tsince/60
   divedata2$deep.tnext <- divedata$deep.tnext/60
   divedata2$deep.postdur <- divedata$deep.postdur/60
   divedata2$deep.dur <- divedata$deep.dur/60
   divedata2$deep.ndur<- divedata$deep.ndur/60
     
   save(divedata2, 
        filter1, filter2,
        postwindow1,
        file="Step1_dive_recovery_data.Rd")
   
   
## Summarize data per individual
####################################################
   
  # load dive-by-dive data
  load("alldata.Rd")
  divedata0 <- divedata
  load("Step1_dive_recovery_data.Rd")
  divedata <- divedata2
  remove(divedata2)

  itab <- indsummary[,c("ind","tagstart_str")] 
  itab$h_total <- NA
  itab$h_baseline_post <- NA
  itab$exposures <- ""
  
  itab$type <- indsummary$type
  itab$TBM <- indsummary$TBM.m
  itab$flukefs_m <- indsummary$flukefs_m
  itab$col <-indsummary$col
    
  itab$IBR <- NA # median IBR (min-1)
  
  itab$d_n <- NA  # number of dives with depth depth > TH
  itab$d_np <- NA  # number of dives per hour
  itab$d_tp <- NA  # proportion of time spent in dives > TH
  
  itab$d_depth_mean <- NA # mean dive depth
  itab$d_depth_max <- NA # maximum dive depth
  
  itab$d_dur_mean <- NA # mean dive duration
  itab$d_dur_sd <- NA   # sd of dive duration
  itab$d_dur_max <- NA # max dive duration

  cycleh <- (divedata$dive.dur+divedata$surf.dur)/60/60
  
  whales <- indsummary$ind
  
  for(w in 1:length(whales)) {
    
    wBool <- divedata$ind==whales[w]
    eBool <- wBool & (divedata$Session=="Baseline" | divedata$Session=="Post")
    sBool <- wBool & eBool & !divedata$deep.filter
    dBool <- wBool & eBool & divedata$deep.filter
      
    itab$h_total[w] <- sum(cycleh[wBool])
    itab$h_baseline_post[w] <- sum(cycleh[eBool])
    itab$exposures[w] <- paste(unique(as.character(divedata$Session[wBool])), collapse=" ")
    
    itab$IBR[w] <- median(1/(divedata$dive.dur[sBool]/60)) # median IBR (min-1)

    itab$d_n[w] <- sum(dBool)  # number of dives with depth depth > TH
    itab$d_np[w] <- sum(dBool)/sum(eBool*cycleh) # number of dives per hour
    itab$d_tp[w] <- 100*sum(dBool*cycleh)/sum(eBool*cycleh)  # proportion of time spent in dives > TH
    
    if(sum(dBool)>0) {
      
    itab$d_depth_mean[w] <- mean(divedata$dive.depth[dBool]) # mean dive depth
    itab$d_depth_max[w] <- max(divedata$dive.depth[dBool]) # maximum dive depth
    
    itab$d_dur_mean[w] <- mean(divedata$dive.dur[dBool]/60) # mean dive duration (min)
    itab$d_dur_sd[w] <- sd(divedata$dive.dur[dBool]/60)   # sd of dive duration (min)
    itab$d_dur_max[w] <- max(divedata$dive.dur[dBool]/60) # max dive duration
    
    max_dur_index <- which(dBool)[which.max(divedata$dive.dur[dBool])]
    max_depth_index <- which(dBool)[which.max(divedata$dive.depth[dBool])]
    
    }
  }
  
  itab$col <- indsummary$col
  itab$pch <- NA
  itab$pch[itab$type=="S" | itab$type=="SC"] <- unique(divedata$pch[divedata$ind.type=="S"])
  itab$pch[itab$type=="M" | itab$type=="MC"] <- unique(divedata$pch[divedata$ind.type=="M"])
  itab$pch[itab$type=="L" | itab$type=="LC"] <- unique(divedata$pch[divedata$ind.type=="L"])
  itab$type.col <- "black"
  itab$type.col[itab$type=="SC" | itab$type=="MC" | itab$type=="LC"] <- "darkgrey"
  

  ### Supplementary figure about individual size classes

  l <- indsummary$flukefs_m-(indsummary$flukefs_m-indsummary$flukefs_l)/2
  u <- indsummary$flukefs_m+(indsummary$flukefs_u-indsummary$flukefs_m)/2
  indsummary$TBM.l <- exp(log(l/3.56)/-0.29) 
  indsummary$TBM.m <- exp(log(indsummary$flukefs_m/3.56)/-0.29) 
  indsummary$TBM.u <- exp(log(u/3.56)/-0.29) 
  
  par(mfrow=c(1,1), mar=c(4,4,1,4))
  
  mypch <- 14+as.numeric(indsummary$Size2)
  
  tempI <- order(indsummary$Size2, indsummary$TBM.m)
  tBool <- indsummary$Calf[tempI]==1
  
  plot(1:length(tempI), indsummary$TBM.m[tempI], 
       pch=mypch[tempI], ylim=c(100,2000),
       xlab="Individual", ylab="Estimated TBM (kg)")
  segments(x0=1:length(tempI),x1=1:length(tempI),
           y0=indsummary$TBM.l[tempI], y1=indsummary$TBM.u[tempI])
  points((1:length(tempI))[tBool], indsummary$TBM.m[tempI][tBool],col="white",pch=3)
  
  legend(1,2000, c("Small","Medium","Large","W. Calf","Dorsal fin size"), 
         pch=c(14+(1:3),3,5),col=c(rep("black",4),"darkgrey"))
  
  par(new=T)
  plot(1:length(tempI), indsummary$Fin.width[tempI],
       xlab="",ylab="",yaxt="n",xaxt="n", pch=5, col="darkgrey")
  axis(4,col="darkgrey")
  mtext(side=4,line=3, "Dorsal fin size (cm)")
  
  

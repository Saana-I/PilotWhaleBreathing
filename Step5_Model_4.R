## R-tools
####################################################

  library(mgcv)
  
  myRsq <- function(obs, pred, type="CoD") {
    n <- length(obs)
    SStot <- sum((obs-mean(obs))^2)
    SSreg <- sum((pred-mean(obs))^2)
    SSres <- sum((obs-pred)^2)
    if(type=="CoD") {return(1-SSres/SStot)}
    if(type=="EV") {return((SSreg/n)/(SStot/n))}
  }

## Data
####################################################
  
  divedata <- read.csv("divedata.csv")
  load("Step4_Model_3_expected_IBR.Rd")
  
  ### Transform data
  divedata$ind <- as.character(divedata$ind)
  divedata$ind.type <- as.character(divedata$ind.type)
  divedata$ind.col <- as.character(divedata$ind.col)
  divedata$Session <- as.character(divedata$Session)
  divedata$dive.GMTtime <- as.POSIXct(as.character(divedata$dive.GMTtime), tz="GMT")

  divedata$dive.rate <- 1/(divedata$dive.dur/60) # response variable
  divedata$pred.ibr <- D # expected IBR
  divedata$dive.flukerate <- divedata$dive.fluken/(divedata$dive.dur/60)
  divedata$ind.calf <- as.numeric(divedata$ind.calf)
  divedata$ind.small <- as.numeric(divedata$ind.small)
  divedata$ind.large <- as.numeric(divedata$ind.large)
  
  ### Set up data filters
  bBool <- divedata$dive.SIL+divedata$dive.PB_BBN+divedata$dive.MLFAS+divedata$dive.PB_SON+divedata$dive.PB_KWF+divedata$dive.PB_KWM+divedata$dive.PB_HW
  bBool <- as.numeric(bBool)==0  
  fBool <- divedata$ind!="gm09_137c" & divedata$ind!="gm09_138b" & divedata$ind!="gm13_169b" & divedata$ind!="gm14_180b"
  filter4 <- fBool & divedata$dive.depth < 31 & divedata$dive.PB_HW==0
  
  ### Set up covariate data
  dur_TH <- 5 # minutes
  RL_TH <- 70 # dB SPL
  divedata$RL_max2 <- divedata$RL_max-RL_TH
  divedata$RL_max[is.na(divedata$RL_max)] <- 0
  divedata$RL_max2[is.na(divedata$RL_max2)] <- 0
  divedata$dive.MLFAS <- divedata$dive.MFAS + divedata$dive.LFAS
  divedata$tsince.MLFAS <- pmin(divedata$tsince.LFAS, divedata$tsince.MFAS)
  divedata$SIL <- as.numeric(divedata$tsince.SIL<=dur_TH)
  divedata$LFAS <- as.numeric(divedata$tsince.LFAS<=dur_TH)
  divedata$MFAS <- as.numeric(divedata$tsince.MFAS<=dur_TH)
  divedata$MLFAS <- as.numeric(divedata$tsince.MLFAS<=dur_TH)
  divedata$PB_BBN <- as.numeric(divedata$tsince.PB_BBN<=dur_TH)
  divedata$PB_SON <- as.numeric(divedata$tsince.PB_SON<=dur_TH)
  divedata$PB_KWM <- as.numeric(divedata$tsince.PB_KWM<=dur_TH)
  divedata$PB_KWF <- as.numeric(divedata$tsince.PB_KWF<=dur_TH)
  divedata$RL_max_LFAS <- divedata$RL_max2
  divedata$RL_max_LFAS[!divedata$LFAS] <- 0
  divedata$RL_max_MFAS <- divedata$RL_max2
  divedata$RL_max_MFAS[!divedata$MFAS] <- 0
  

  ### Order effects
  whales <- unique(divedata$ind)
  firstRL <- rep(FALSE, length(divedata$ind))
  for(w in 1:length(whales)) {
    tempI <- which(divedata$ind==whales[w])
    tFirst <- TRUE
    firstRL[tempI[1]] <- tFirst
    for(j in 2:length(tempI)) {
      if(divedata$RL_max2[tempI[j-1]]>0 & divedata$RL_max2[tempI[j]]==0) {
        tFirst <- FALSE
      }
      firstRL[tempI[j]] <- tFirst
    }
  }
  
  divedata$RL_max_first <- divedata$RL_max2
  divedata$RL_max_first[!firstRL] <- 0
  
  divedata$MLFAS_first <- divedata$MLFAS
  divedata$MLFAS_first[!firstRL] <- 0
  
  divedata$MLFAS_subs <- rep(0, length(divedata$MLFAS))
  divedata$MLFAS_subs[divedata$MLFAS & !firstRL] <- 1
  
  
  
## Compare alternative full models
####################################################
  
  # Full model 1 (combined effect of LFAS and MLFAS)
  temp <- gamm(dive.rate ~ pred.ibr + SIL + PB_BBN + PB_SON + PB_KWF + PB_KWM + MLFAS_subs
               + s(RL_max2, bs="ts",k=3),
               family=Gamma(link = "identity"),
               data=divedata[filter4,],
               random=list(ind=~1),
               correlation=corAR1())
  summary(temp$gam)
  # R-sq.(adj) =  0.328   
  anova(temp$gam)
  # pred.ibr    1 1034.445 < 2e-16
  # SIL         1    1.780 0.18217
  # PB_BBN      1    0.325 0.56840
  # PB_SON      1    0.396 0.52922
  # PB_KWF      1    0.370 0.54319
  # PB_KWM      1    7.276 0.00700
  # MLFAS_subs  1    8.081 0.00448
  
  # Approximate significance of smooth terms:
  #             edf     Ref.df    F  p-value
  # s(RL_max2) 0.8849 2.0000 3.844 0.00316
  
  # Full model 2 (separate effect of LFAS and MLFAS)
  temp <- gamm(dive.rate ~ pred.ibr + SIL + PB_BBN + PB_SON + PB_KWF + PB_KWM + MLFAS_subs
               + s(RL_max_LFAS, bs="ts",k=3) + s(RL_max_MFAS, bs="ts",k=3),
               family=Gamma(link = "identity"),
               data=divedata[filter4,],
               random=list(ind=~1),
               correlation=corAR1())
  summary(temp$gam)
  # R-sq.(adj) =  0.329   
  anova(temp$gam)
  
  #pred.ibr    1 1027.945 < 2e-16
  #SIL         1    1.792 0.18076
  #PB_BBN      1    0.320 0.57146
  #PB_SON      1    0.394 0.53001
  #PB_KWF      1    0.372 0.54216
  #PB_KWM      1    7.268 0.00703
  #MLFAS_subs  1    7.239 0.00714
  
  # Approximate significance of smooth terms:
  # edf            Ref.df    F  p-value
  # s(RL_max_LFAS) 9.167e-01 2.000e+00 5.408 0.000592
  # s(RL_max_MFAS) 8.146e-06 2.000e+00 0.000 0.808154
  
  # Full model 3 (separate effect of combined & LFAS)
  temp <- gamm(dive.rate ~ pred.ibr + SIL + PB_BBN + PB_SON + PB_KWF + PB_KWM + MLFAS_subs
               + s(RL_max2, bs="ts",k=3) + s(RL_max_MFAS, bs="ts",k=3),
               family=Gamma(link = "identity"),
               data=divedata[filter4,],
               random=list(ind=~1),
               correlation=corAR1())
  summary(temp$gam)
  # R-sq.(adj) =  0.329   
  anova(temp$gam)
  # pred.ibr    1 1029.086 < 2e-16
  # SIL         1    1.797 0.18007
  # PB_BBN      1    0.322 0.57046
  # PB_SON      1    0.395 0.52968
  # PB_KWF      1    0.365 0.54557
  # PB_KWM      1    7.269 0.00703
  # MLFAS_subs  1    5.999 0.01433
  
  # Approximate significance of smooth terms:
  #                edf     Ref.df    F  p-value
  # s(RL_max2)     0.8980 2.0000 4.350 0.00176
  # s(RL_max_MFAS) 0.6167 2.0000 0.829 0.10047
  
## Model selection (step-wise, backwards)
####################################################

  # Start from full model 2 (separate effect of LFAS and MLFAS)
  temp <- gamm(dive.rate ~ pred.ibr + SIL + PB_BBN + PB_SON + PB_KWF + PB_KWM + MLFAS_subs
               + s(RL_max_LFAS, bs="ts",k=3) + s(RL_max_MFAS, bs="ts",k=3),
               family=Gamma(link = "identity"),
               data=divedata[filter4,],
               random=list(ind=~1),
               correlation=corAR1())
  anova(temp$gam)
  #pred.ibr    1 1027.945 < 2e-16
  #SIL         1    1.792 0.18076
  #PB_BBN      1    0.320 0.57146
  #PB_SON      1    0.394 0.53001
  #PB_KWF      1    0.372 0.54216
  #PB_KWM      1    7.268 0.00703
  #MLFAS_subs  1    7.239 0.00714
  
  # Approximate significance of smooth terms:
  #                edf     Ref.df    F  p-value
  # s(RL_max_LFAS) 9.167e-01 2.000e+00 5.408 0.000592
  # s(RL_max_MFAS) 8.146e-06 2.000e+00 0.000 0.808154 << drop this next
  
  # Step 1
  temp <- gamm(dive.rate ~ pred.ibr + SIL + PB_BBN + PB_SON + PB_KWF + PB_KWM + MLFAS_subs
               + s(RL_max_LFAS, bs="ts",k=3),
               family=Gamma(link = "identity"),
               data=divedata[filter4,],
               random=list(ind=~1),
               correlation=corAR1())
  anova(temp$gam)
  
  # pred.ibr    1 1027.947 < 2e-16
  # SIL         1    1.792 0.18076
  # PB_BBN      1    0.320 0.57146 << drop this next
  # PB_SON      1    0.394 0.53001
  # PB_KWF      1    0.372 0.54216
  # PB_KWM      1    7.268 0.00703
  # MLFAS_subs  1    7.239 0.00714
  
  # Approximate significance of smooth terms:
  #             edf     Ref.df    F  p-value
  # s(RL_max_LFAS) 0.9167 2.0000 5.408 0.000592

  # Step 2
  temp <- gamm(dive.rate ~ pred.ibr + SIL + PB_SON + PB_KWF + PB_KWM + MLFAS_subs
               + s(RL_max_LFAS, bs="ts",k=3),
               family=Gamma(link = "identity"),
               data=divedata[filter4,],
               random=list(ind=~1),
               correlation=corAR1())
  anova(temp$gam)
  # pred.ibr    1 1028.275 < 2e-16
  # SIL         1    1.798 0.17993
  # PB_SON      1    0.389 0.53284 
  # PB_KWF      1    0.365 0.54549 << drop this next
  # PB_KWM      1    7.238 0.00715
  # MLFAS_subs  1    7.250 0.00710
  
  # Approximate significance of smooth terms:
  #             edf     Ref.df    F  p-value
  # s(RL_max_LFAS) 0.9168 2.0000 5.41 0.000591
  
  # Step 3
  temp <- gamm(dive.rate ~ pred.ibr + SIL + PB_SON + PB_KWM + MLFAS_subs
               + s(RL_max_LFAS, bs="ts",k=3),
               family=Gamma(link = "identity"),
               data=divedata[filter4,],
               random=list(ind=~1),
               correlation=corAR1())
  anova(temp$gam)
  # pred.ibr    1 1027.007 < 2e-16
  # SIL         1    1.868 0.17176
  # PB_SON      1    0.389 0.53290 << drop this next
  # PB_KWM      1    7.235 0.00716
  # MLFAS_subs  1    7.326 0.00681
  
  # Approximate significance of smooth terms:
  #             edf     Ref.df    F  p-value
  # s(RL_max_LFAS) 0.916  2.000 5.355 0.000626
  
  # Step 4
  temp <- gamm(dive.rate ~ pred.ibr + SIL + PB_KWM + MLFAS_subs
               + s(RL_max_LFAS, bs="ts",k=3),
               family=Gamma(link = "identity"),
               data=divedata[filter4,],
               random=list(ind=~1),
               correlation=corAR1())
  anova(temp$gam)
  # pred.ibr    1 1027.353 < 2e-16
  # SIL         1    1.867 0.17181 << drop this next
  # PB_KWM      1    7.064 0.00787
  # MLFAS_subs  1    7.325 0.00681
  
  # Approximate significance of smooth terms:
  #             edf     Ref.df    F  p-value
  # s(RL_max_LFAS) 0.9159 2.0000 5.353 0.000626
  
## Best model interpretation
####################################################
  
  fit1 <- gamm(dive.rate ~ pred.ibr + PB_KWM + MLFAS_subs
               + s(RL_max_LFAS, bs="ts",k=3),
               family=Gamma(link = "identity"),
               data=divedata[filter4,],
               random=list(ind=~1),
               correlation=corAR1())
  summary(fit1$gam)
  # R-sq.(adj) =  0.329   
  anova(fit1$gam)
  # pred.ibr    1 1068.236 < 2e-16
  # PB_KWM      1    7.089 0.00777
  # MLFAS_subs  1    7.117 0.00764
  
  # Approximate significance of smooth terms:
  #             edf     Ref.df    F  p-value
  # s(RL_max_LFAS) 0.9189 2.0000 5.557 0.000503
  
  ### Check model fit
  bBool <- divedata$Session=="Baseline"
  eBool <- (divedata$dive.SIL+divedata$dive.PB_BBN+divedata$dive.MLFAS+divedata$dive.PB_SON+divedata$dive.PB_KWF+divedata$dive.PB_KWM+divedata$dive.PB_HW)>0.5
  pBool <- !bBool & !eBool
  myRsq(divedata$dive.rate[filter4 & bBool], predict(fit1$gam, type="response")[bBool[filter4]]) # 0.3238379
  myRsq(divedata$dive.rate[filter4 & !bBool], predict(fit1$gam, type="response")[!bBool[filter4]]) # 0.3311242
  myRsq(divedata$dive.rate[filter4 & eBool], predict(fit1$gam, type="response")[eBool[filter4]]) # 0.4304414
  myRsq(divedata$dive.rate[filter4 & pBool], predict(fit1$gam, type="response")[pBool[filter4]]) # 0.2724515
  
  ### Plot fitted values
  plot(divedata$RL_max[filter4], divedata$dive.rate[filter4], 
       xlab="Received SPL max", ylab="IBR (min-1)",
       col="grey", ylim=c(0,15))
  preddata <- divedata[filter4,]
  preddata$MLFAS_subs <- 0
  fitted_vals <- predict(fit1$gam, newdata=preddata,type="response")
  points(divedata$RL_max[filter4], fitted_vals, col="red")
  abline(h=3, col="blue")
  
  ### Model predictions - baseline, PB_KWM and max RL_max_LFAS (180)
  preddata <- data.frame(RL_max_LFAS = c(0,0,180-RL_TH,0))
  preddata$MLFAS_subs <- c(0,0,0,1)
  preddata$PB_KWM <- c(0,1,0,0)
  preddata$pred.ibr <- median(divedata$pred.ibr[filter4])
  preddata$preds <- predict(fit1$gam, newdata=preddata, type="response", se.fit=F)
  preddata$preds_se <- predict(fit1$gam, newdata=preddata, type="response", se.fit=T)$se
  print(preddata)
  #  RL_max_LFAS MLFAS_subs PB_KWM pred.ibr    preds  preds_se
  # 1           0          0      0 3.099511 3.085339 0.1259566
  # 2           0          0      1 3.099511 2.687673 0.1926381
  # 3         110          0      0 3.099511 2.612215 0.1883660
  # 4           0          1      0 3.099511 3.454662 0.1857803
  
  
## Plot predicted breathing rate
####################################################

  ### Figure panels margins
  par(mfrow=c(2,2), mar=c(5, 4.5, 3, 1), lwd=0.8, cex=0.7, cex.lab=1.1)

  ### 1. Predict for breathing rate during non-recovery periods
  
  ### filter_new: recovery periods defined as rec_TH x above
  ### minimum predicted for each tag
  rec_TH <- 1.1 # multiplier for recovery
  tBool <- rep(FALSE, length(divedata$ind))
  whales <- unique(divedata$ind)
  tTH <- 0
  for(w in 1:length(whales)) {
    wBool <- divedata$ind==whales[w]
    tTH[w] <- min(divedata$pred.ibr[wBool])
    tBool[wBool & divedata$pred.ibr > rec_TH*tTH[w]] <- TRUE
  }
  
  filter_new <- filter4 & !tBool 
  ibr_median <- median(tapply(divedata$pred.ibr[filter_new], divedata$ind[filter_new], median))
  
  preddata <- data.frame(RL_max_LFAS = 0:round(max(divedata$RL_max_LFAS,na.rm=T)))
  preddata$MLFAS_subs <- 0
  preddata$PB_KWM <- 0
  preddata$pred.ibr <- ibr_median
  preds <- predict(fit1$gam, newdata=preddata, type="response", se.fit=T)
  
  tval <- c(40,50,60) # baseline, PB_KWM, MLFAS_subs
  temp <- divedata$RL_max
  temp[divedata$RL_max_LFAS==0] <- NA
  temp[divedata$Session=="Baseline"] <- tval[1]
  temp[divedata$PB_KWM==1] <- tval[2]
  temp[divedata$MLFAS_subs==1] <- tval[3]
  
  plot(temp[filter_new], divedata$dive.rate[filter_new], col="grey", 
       xlab=expression(paste(plain("      SPL"),  
                             scriptstyle("max"),plain(" (dB re 1"), mu,"Pa)", sep="")), 
       ylab=expression(paste("Breathing rate (min"^"-1",")")), 
       xaxt="n", ylim=c(0,15), xlim=c(tval[1],180))
  temp2 <- divedata$RL_max
  temp2[divedata$RL_max_LFAS==0] <- NA
  lines(temp2[filter_new], divedata$dive.rate[filter_new], col="grey")
  
  lines(preddata$RL_max+RL_TH, preds$fit, type="l")
  lines(preddata$RL_max+RL_TH, preds$fit-qnorm(0.975)*preds$se.fit, col=1, lty=2)
  lines(preddata$RL_max+RL_TH, preds$fit+qnorm(0.975)*preds$se.fit, col=1, lty=2)
  axis(1,at=seq(80,180,10))
  axis(1, at=tval, labels=c("Baseline","PB_KWM", "SON_2"), las=2)
  title("a. Non-recovery periods")
  
  # baseline response
  points(tval[1], preds$fit[1], pch=15)
  bval <- preds$fit[1]
  segments(x0=tval[1], x1=tval[1],
           y0=preds$fit[1]-qnorm(0.975)*preds$se.fit[1],
           y1=preds$fit[1]+qnorm(0.975)*preds$se.fit[1], col=1)
  
  preds$fit[1] # 2.845524 
  preds$se.fit[1] # 0.126282 
  preds$fit[1]-qnorm(0.975)*preds$se.fit[1] # 2.598016 
  preds$fit[1]+qnorm(0.975)*preds$se.fit[1] # 3.093032 
  
  preds$fit[which.max(preddata$RL_max)] # 2.3724 
  preds$se.fit[which.max(preddata$RL_max)] # 0.1883354 
  max(divedata$RL_max[divedata$LFAS>0], na.rm=T) # 179.8408
  max(preddata$RL_max) # 180
  
  ## PB_KWM
  tvalI <- 2
  preddata2 <- preddata[1,]
  preddata2$PB_KWM <- 1
  preds2 <- predict(fit1$gam, newdata=preddata2, type="response", se.fit=T)
  points(tval[tvalI], preds2$fit, pch=15)
  segments(x0=tval[tvalI], x1=tval[tvalI],
           y0=preds2$fit-qnorm(0.975)*preds2$se.fit,
           y1=preds2$fit+qnorm(0.975)*preds2$se.fit, col=1)
  
  ## SON-2
  tvalI <- 3
  preddata2 <- preddata[1,]
  preddata2$MLFAS_subs <- 1
  preds2 <- predict(fit1$gam, newdata=preddata2, type="response", se.fit=T)
  points(tval[tvalI], preds2$fit, pch=15)
  segments(x0=tval[tvalI], x1=tval[tvalI],
           y0=preds2$fit-qnorm(0.975)*preds2$se.fit,
           y1=preds2$fit+qnorm(0.975)*preds2$se.fit, col=1)
  
  #abline(h=3, col="blue", lty=3)
  abline(h=bval, col="blue", lty=3)
  
  ### 2. Predict for breathing rate during recovery
  
  filter_new <- filter4 & tBool 
  ibr_median <- median(tapply(divedata$pred.ibr[filter_new], divedata$ind[filter_new], median))
  preddata <- data.frame(RL_max_LFAS = 0:round(max(divedata$RL_max_LFAS,na.rm=T)))
  preddata$MLFAS_subs <- 0
  preddata$PB_KWM <- 0
  preddata$pred.ibr <- ibr_median
  preds <- predict(fit1$gam, newdata=preddata, type="response", se.fit=T)
  
  #####################
  temp <- divedata$RL_max
  tval <- c(40,50,60) # baseline, PB_KWM, MLFAS_subs
  temp <- divedata$RL_max
  temp[divedata$RL_max_LFAS==0] <- NA
  temp[divedata$Session=="Baseline"] <- tval[1]
  temp[divedata$PB_KWM==1] <- tval[2]
  temp[divedata$MLFAS_subs==1] <- tval[3]
  
  plot(temp[filter_new], divedata$dive.rate[filter_new], col="grey", 
       xlab=expression(paste(plain("      SPL"),  
                             scriptstyle("max"),plain(" (dB re 1"), mu,"Pa)", sep="")), 
       ylab=expression(paste("Breathing rate (min"^"-1",")")), 
       xaxt="n", ylim=c(0,15), xlim=c(tval[1],180))
  temp2 <- divedata$RL_max
  temp2[divedata$RL_max_LFAS==0] <- NA
  lines(temp2[filter_new], divedata$dive.rate[filter_new], col="grey")
  
  lines(preddata$RL_max+RL_TH, preds$fit, type="l")
  lines(preddata$RL_max+RL_TH, preds$fit-qnorm(0.975)*preds$se.fit, col=1, lty=2)
  lines(preddata$RL_max+RL_TH, preds$fit+qnorm(0.975)*preds$se.fit, col=1, lty=2)
  axis(1,at=seq(80,180,10))
  axis(1, at=tval, labels=c("Baseline","PB_KWM", "SON_2"), las=2)
  title("b. Recovery periods")
  
  # baseline response
  points(tval[1], preds$fit[1], pch=15)
  bval <- preds$fit[1]
  segments(x0=tval[1], x1=tval[1],
           y0=preds$fit[1]-qnorm(0.975)*preds$se.fit[1],
           y1=preds$fit[1]+qnorm(0.975)*preds$se.fit[1], col=1)
  
  preds$fit[1] # 3.838985 
  preds$se.fit[1] # 0.1277029 
  preds$fit[1]-qnorm(0.975)*preds$se.fit[1] # 3.588692 
  preds$fit[1]+qnorm(0.975)*preds$se.fit[1] # 4.089279 
  
  preds$fit[which.max(preddata$RL_max)] # 3.365861 
  preds$se.fit[which.max(preddata$RL_max)] # 0.1903124 
  max(preddata$RL_max_LFAS) # 110
  
  ## PB_KWM
  tvalI <- 2
  preddata2 <- preddata[1,]
  preddata2$PB_KWM <- 1
  preds2 <- predict(fit1$gam, newdata=preddata2, type="response", se.fit=T)
  points(tval[tvalI], preds2$fit, pch=15)
  segments(x0=tval[tvalI], x1=tval[tvalI],
           y0=preds2$fit-qnorm(0.975)*preds2$se.fit,
           y1=preds2$fit+qnorm(0.975)*preds2$se.fit, col=1)
  
  ## SON-2
  tvalI <- 3
  preddata2 <- preddata[1,]
  preddata2$MLFAS_subs <- 1
  preds2 <- predict(fit1$gam, newdata=preddata2, type="response", se.fit=T)
  points(tval[tvalI], preds2$fit, pch=15)
  segments(x0=tval[tvalI], x1=tval[tvalI],
           y0=preds2$fit-qnorm(0.975)*preds2$se.fit,
           y1=preds2$fit+qnorm(0.975)*preds2$se.fit, col=1)
  
  #abline(h=3, col="blue", lty=3)
  abline(h=bval, col="blue", lty=3)
  
  
  
  
  ### 3. Predict for breathing rate during non-recovery periods - correct for expected rate
  
  ### filter_new: recovery periods defined as rec_TH x above
  ### minimum predicted for each tag
  rec_TH <- 1.1 # multiplier for recovery
  tBool <- rep(FALSE, length(divedata$ind))
  whales <- unique(divedata$ind)
  tTH <- 0
  for(w in 1:length(whales)) {
    wBool <- divedata$ind==whales[w]
    tTH[w] <- min(divedata$pred.ibr[wBool])
    tBool[wBool & divedata$pred.ibr > rec_TH*tTH[w]] <- TRUE
  }
  
  filter_new <- filter4 & !tBool 
  ibr_median <- median(tapply(divedata$pred.ibr[filter_new], divedata$ind[filter_new], median))
  
  preddata <- data.frame(RL_max_LFAS = 0:round(max(divedata$RL_max_LFAS,na.rm=T)))
  preddata$MLFAS_subs <- 0
  preddata$PB_KWM <- 0
  preddata$pred.ibr <- ibr_median
  preds <- predict(fit1$gam, newdata=preddata, type="response", se.fit=T)
  preds$fit <- preds$fit-preddata$pred.ibr
  
  tval <- c(40,50,60) # baseline, PB_KWM, MLFAS_subs
  temp <- divedata$RL_max
  temp[divedata$RL_max_LFAS==0] <- NA
  temp[divedata$Session=="Baseline"] <- tval[1]
  temp[divedata$PB_KWM==1] <- tval[2]
  temp[divedata$MLFAS_subs==1] <- tval[3]
  
  res <- divedata$dive.rate[filter_new]-divedata$pred.ibr[filter_new]
  plot(temp[filter_new], res, col="grey", 
       xlab=expression(paste(plain("      SPL"),  
                             scriptstyle("max"),plain(" (dB re 1"), mu,"Pa)", sep="")), 
       ylab=expression(paste("Obs - expected IBR (min"^"-1",")")), 
       xaxt="n", ylim=c(-3,10), xlim=c(tval[1],180))
  temp2 <- divedata$RL_max
  temp2[divedata$RL_max_LFAS==0] <- NA
  lines(temp2[filter_new], res, col="grey")
  
  lines(preddata$RL_max+RL_TH, preds$fit, type="l")
  lines(preddata$RL_max+RL_TH, preds$fit-qnorm(0.975)*preds$se.fit, col=1, lty=2)
  lines(preddata$RL_max+RL_TH, preds$fit+qnorm(0.975)*preds$se.fit, col=1, lty=2)
  axis(1,at=seq(80,180,10))
  axis(1, at=tval, labels=c("Baseline","PB_KWM", "SON_2"), las=2)
  title("a. Non-recovery periods")
  
  # baseline response
  points(tval[1], preds$fit[1], pch=15)
  bval <- preds$fit[1]
  segments(x0=tval[1], x1=tval[1],
           y0=preds$fit[1]-qnorm(0.975)*preds$se.fit[1],
           y1=preds$fit[1]+qnorm(0.975)*preds$se.fit[1], col=1)
  
  ## PB_KWM
  tvalI <- 2
  preddata2 <- preddata[1,]
  preddata2$PB_KWM <- 1
  preds2 <- predict(fit1$gam, newdata=preddata2, type="response", se.fit=T)
  preds2$fit <- preds2$fit-preddata2$pred.ibr
  points(tval[tvalI], preds2$fit, pch=15)
  segments(x0=tval[tvalI], x1=tval[tvalI],
           y0=preds2$fit-qnorm(0.975)*preds2$se.fit,
           y1=preds2$fit+qnorm(0.975)*preds2$se.fit, col=1)
  
  ## SON-2
  tvalI <- 3
  preddata2 <- preddata[1,]
  preddata2$MLFAS_subs <- 1
  preds2 <- predict(fit1$gam, newdata=preddata2, type="response", se.fit=T)
  preds2$fit <- preds2$fit-preddata2$pred.ibr
  points(tval[tvalI], preds2$fit, pch=15)
  segments(x0=tval[tvalI], x1=tval[tvalI],
           y0=preds2$fit-qnorm(0.975)*preds2$se.fit,
           y1=preds2$fit+qnorm(0.975)*preds2$se.fit, col=1)
  
  #abline(h=3, col="blue", lty=3)
  abline(h=bval, col="blue", lty=3)
  
  ### 4. Predict for breathing rate during recovery - correct for expected rate
  
  filter_new <- filter4 & tBool 
  ibr_median <- median(tapply(divedata$pred.ibr[filter_new], divedata$ind[filter_new], median))
  preddata <- data.frame(RL_max_LFAS = 0:round(max(divedata$RL_max_LFAS,na.rm=T)))
  preddata$MLFAS_subs <- 0
  preddata$PB_KWM <- 0
  preddata$pred.ibr <- ibr_median
  preds <- predict(fit1$gam, newdata=preddata, type="response", se.fit=T)
  preds$fit <- preds$fit-preddata$pred.ibr
  
  #####################
  temp <- divedata$RL_max
  tval <- c(40,50,60) # baseline, PB_KWM, MLFAS_subs
  temp <- divedata$RL_max
  temp[divedata$RL_max_LFAS==0] <- NA
  temp[divedata$Session=="Baseline"] <- tval[1]
  temp[divedata$PB_KWM==1] <- tval[2]
  temp[divedata$MLFAS_subs==1] <- tval[3]
  
  
  res <- divedata$dive.rate[filter_new]-divedata$pred.ibr[filter_new]
  plot(temp[filter_new], res, col="grey", 
       xlab=expression(paste(plain("      SPL"),  
                             scriptstyle("max"),plain(" (dB re 1"), mu,"Pa)", sep="")), 
       ylab=expression(paste("Obs - expected IBR (min"^"-1",")")), 
       xaxt="n", ylim=c(-3,10), xlim=c(tval[1],180))
  temp2 <- divedata$RL_max
  temp2[divedata$RL_max_LFAS==0] <- NA
  lines(temp2[filter_new], res, col="grey")
  
  lines(preddata$RL_max+RL_TH, preds$fit, type="l")
  lines(preddata$RL_max+RL_TH, preds$fit-qnorm(0.975)*preds$se.fit, col=1, lty=2)
  lines(preddata$RL_max+RL_TH, preds$fit+qnorm(0.975)*preds$se.fit, col=1, lty=2)
  axis(1,at=seq(80,180,10))
  axis(1, at=tval, labels=c("Baseline","PB_KWM", "SON_2"), las=2)
  title("b. Recovery periods")
  
  # baseline response
  points(tval[1], preds$fit[1], pch=15)
  bval <- preds$fit[1]
  segments(x0=tval[1], x1=tval[1],
           y0=preds$fit[1]-qnorm(0.975)*preds$se.fit[1],
           y1=preds$fit[1]+qnorm(0.975)*preds$se.fit[1], col=1)
  
  ## PB_KWM
  tvalI <- 2
  preddata2 <- preddata[1,]
  preddata2$PB_KWM <- 1
  preds2 <- predict(fit1$gam, newdata=preddata2, type="response", se.fit=T)
  preds2$fit <- preds2$fit-preddata2$pred.ibr
  points(tval[tvalI], preds2$fit, pch=15)
  segments(x0=tval[tvalI], x1=tval[tvalI],
           y0=preds2$fit-qnorm(0.975)*preds2$se.fit,
           y1=preds2$fit+qnorm(0.975)*preds2$se.fit, col=1)
  
  ## SON-2
  tvalI <- 3
  preddata2 <- preddata[1,]
  preddata2$MLFAS_subs <- 1
  preds2 <- predict(fit1$gam, newdata=preddata2, type="response", se.fit=T)
  preds2$fit <- preds2$fit-preddata2$pred.ibr
  points(tval[tvalI], preds2$fit, pch=15)
  segments(x0=tval[tvalI], x1=tval[tvalI],
           y0=preds2$fit-qnorm(0.975)*preds2$se.fit,
           y1=preds2$fit+qnorm(0.975)*preds2$se.fit, col=1)
  
  #abline(h=3, col="blue", lty=3)
  abline(h=bval, col="blue", lty=3)
  
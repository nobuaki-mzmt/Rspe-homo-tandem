## R. speratus homosexual tandem analysis
## N. Mizumoto
## 7/12/2021 ~
## script with R projects

## This file is for analyzing datasets during tandems
## 1. compare duration of tandems (or sep search) across pair combinations
## 2. examine leader-follower interactions: interactions are estimated using the correlation between acceleration and inter-individual distance. (Franks and Richardson 2006; Mizumoto and Bourguignon in prep).

# Environment -------------------------------------------------------------
{
  RPROJ <- list(PROJHOME = normalizePath(getwd()))
  attach(RPROJ)
  rm(RPROJ)
}

# Setup -------------------------------------------------------------------
{
  rm(list = ls())
  source(file.path(PROJHOME, "/scripts/Functions.R"))
  
  ## Data
  rda.place <- file.path(PROJHOME, "/data/", "rda")
  f.name <- paste0(rda.place, "/AllData-scaled.rda")
  load(f.name)
  
  ## Parameters
  interaction.threshold <- 7
  min.tandem.sec <- 2
  threshold.moving.dis.tandem <- 30
}

# Computation --------------------------------------------------------------
{
  Plot = F
  
  ## get only during tandem
  df.tan <- df[df$scheme != "sep",]
  
  ## calcurate step length
  df.tan$sl <- c(NA, sqrt( diff(df.tan$x)^2 + diff(df.tan$y)^2))
  df.tan$sl[df.tan$time == 0] <- NA
  
  ## calculate accerelation
  df.tan$acc <- c(diff(df.tan$sl*5), NA)
  df.tan$acc[df.tan$time == 0 | df.tan$time == 300] <- NA
  
  ## main calculation
  ind.names <- unique(df.tan$name) 
  res.tan.time = res.sep.time = res.tandem = res.all <- NULL
  for(i in 1:length(ind.names)){
    print(paste(i, "/", length(ind.names), "->", ind.names[i]))
    
    df.temp <- df.tan[df.tan$name == ind.names[i],]
    df.temp <- na.omit(df.temp)
    df.follower <- df.temp[df.temp$role == "Follower",]
    df.leader <- df.temp[df.temp$role == "Leader",]
    if(dim(df.temp)[1] < 2998){print(paste("not enough data", dim(df.temp)[1], "next"));next;}
    x.follower = df.follower$x
    y.follower = df.follower$y
    x.leader = df.leader$x
    y.leader = df.leader$y
    sl.follower = df.follower$sl
    sl.leader = df.leader$sl
    ind.dis <- sqrt( (x.follower - x.leader)^2 + (y.follower - y.leader)^2 )
    
    ## definition of tandem running
    {
      # 1.1. not separation less than min.sec --> tandem
      # 1.2. tandem shorter than min.sec --> separation
      {
        interaction <- ind.dis <= interaction.threshold
        non.interaction <- !interaction
        separation <- ind.dis > interaction.threshold * 2
        tandem = !tandem.smoothing2(!interaction, separation, min.tandem.sec*5)
        tandem <- tandem.smoothing(tandem, !tandem, min.tandem.sec*5)
      }
      
      # 2. tandem: both need to move more than thresh_dis or thresh_speed during interaction
      if(sum(tandem)>0){
        tan.timing <- which(tandem)[c(T, diff(which(tandem))>1)]
        tan.end <- which(tandem)[c(diff(which(tandem))>1,T)]
        for(j in 1:length(tan.timing)){
          if(sum(sl.follower[tan.timing[j]:tan.end[j]], na.rm=T)< threshold.moving.dis.tandem){
            tandem[tan.timing[j]:tan.end[j]] <- F
          }
          if(sum(sl.leader[tan.timing[j]:tan.end[j]], na.rm=T)< threshold.moving.dis.tandem){
            tandem[tan.timing[j]:tan.end[j]] <- F
          }
        }
      }
    }
    
    ## Predict role change
    {
      lv <- length(x.leader)
      angle.leader = atan2(y.leader[2:lv] - y.leader[2:lv-1] , x.leader[2:lv] - x.leader[2:lv-1])
      angle.follower = atan2(y.follower[2:lv] - y.follower[2:lv-1] , x.follower[2:lv] - x.follower[2:lv-1])
      angle.l.to.f = atan2(y.follower - y.leader , x.follower - x.leader)[2:lv-1]
      angle.f.to.l = atan2(y.leader - y.follower , x.leader - x.follower)[2:lv-1]
      r.dir.l.to.f <- apply(cbind(abs(angle.leader - angle.l.to.f), pi*2-abs(angle.leader - angle.l.to.f)), 1, min)
      r.dir.f.to.l <- apply(cbind(abs(angle.follower - angle.f.to.l), pi*2-abs(angle.follower - angle.f.to.l)), 1, min)
      
      correct.role <- r.dir.f.to.l < pi/2 & abs(pi-r.dir.l.to.f) < pi/2
      reverse.role <- r.dir.l.to.f < pi/2 & abs(pi-r.dir.f.to.l) < pi/2
      correct.role <- (c(correct.role,F) & tandem)
      reverse.role <- (c(reverse.role,F) & tandem)
      
      x.leader.vec <- x.leader[2:lv] - x.leader[2:lv-1]
      y.leader.vec <- y.leader[2:lv] - y.leader[2:lv-1]
      leader.vec.length <- sqrt(x.leader.vec^2+y.leader.vec^2)
      x.leader.vec <- x.leader.vec/leader.vec.length
      y.leader.vec <- y.leader.vec/leader.vec.length
      
      x.follower.vec <- x.follower[2:lv] - x.follower[2:lv-1]
      y.follower.vec <- y.follower[2:lv] - y.follower[2:lv-1]
      follower.vec.length <- sqrt(x.follower.vec^2+y.follower.vec^2)
      x.follower.vec <- x.follower.vec/follower.vec.length
      y.follower.vec <- y.follower.vec/follower.vec.length
      
      x.center <- (x.follower + x.leader)/2
      y.center <- (y.follower + y.leader)/2
      x.leader.rotation.vec <- x.leader[2:lv-1] - x.center[2:lv-1]
      y.leader.rotation.vec <- y.leader[2:lv-1] - y.center[2:lv-1]
      leader.rotation.vec.length <- sqrt(x.leader.rotation.vec^2+y.leader.rotation.vec^2)
      x.leader.rotation.vec <- x.leader.rotation.vec/leader.rotation.vec.length
      y.leader.rotation.vec <- y.leader.rotation.vec/leader.rotation.vec.length
      
      x.follower.rotation.vec <- x.follower[2:lv-1] - x.center[2:lv-1]
      y.follower.rotation.vec <- y.follower[2:lv-1] - y.center[2:lv-1]
      follower.rotation.vec.length <- sqrt(x.follower.rotation.vec^2+y.follower.rotation.vec^2)
      x.follower.rotation.vec <- x.follower.rotation.vec/follower.rotation.vec.length
      y.follower.rotation.vec <- y.follower.rotation.vec/follower.rotation.vec.length
      
      x.leader.ang.moment <- x.leader.vec*x.leader.rotation.vec - y.leader.vec*y.leader.rotation.vec
      y.leader.ang.moment <- x.leader.vec*y.leader.rotation.vec + y.leader.vec*x.leader.rotation.vec
      x.follower.ang.moment <- x.follower.vec*x.follower.rotation.vec - y.follower.vec*y.follower.rotation.vec
      y.follower.ang.moment <- x.follower.vec*y.follower.rotation.vec + y.follower.vec*x.follower.rotation.vec
      leader.ang.moment.length <- sqrt(x.leader.ang.moment^2+y.leader.ang.moment^2)
      follower.ang.moment.length <- sqrt(x.follower.ang.moment^2+y.follower.ang.moment^2)
      
      polarization <- sqrt(
        (x.leader.vec + x.follower.vec)^2+
          (y.leader.vec + y.follower.vec)^2
      )/2
      
      rotation <- sqrt(
        (x.leader.ang.moment + x.follower.ang.moment)^2+
          (y.leader.ang.moment + y.follower.ang.moment)^2
      )/2
      if(Plot){
        png(file.path(PROJHOME, "plot/tandem-reverse/", paste0(df.temp$name[1], ".png")))
        par(mfrow=c(2,1), pin=c(6,2))
        plot((tandem)*0, col=(tandem)*2, type="p", ylim=c(-1,1), main=ind.names[i])
        points((correct.role - reverse.role), type="l")
        
        plot(polarization, col=2, type="l", ylim=c(0,1))
        points(rotation, col=4, type="l")
        
        dev.off()
      }
    }
    
    rotation[is.na(rotation)] <- 0
    competition <- tandem & ind.dis < 4 & c(rotation, F) > 0.5
    competition <- tandem.smoothing(competition, !competition, 5)
    par(mfrow=c(2,1), pin=c(6,2))
    plot((tandem)*0, col=(tandem)*2, type="p", ylim=c(-1,1), main=ind.names[i])
    points((correct.role - reverse.role), type="l")
    points(competition*0-0.2, col=competition*4)
    
    plot(polarization, col=2, type="l", ylim=c(0,1))
    points(rotation, col=4, type="l")
    
    
    ## get scheme
    {
      scheme <- tandem
      scheme[interaction] <- "i"
      scheme[non.interaction] <- "r"
      scheme[tandem] <- "t"
      
      if(sum(tandem)>0){
        tan.end <- which(tandem)[c(diff(which(tandem))>1,T)]
        if(tan.end[1]==length(tandem)){  }else {
          sep.begin <- tan.end[tan.end<length(tandem)]+1
          
          ## sep until next interaction
          tan.timing <- which(interaction)[c(T, diff(which(interaction))>1)]
          if(length(tan.timing)==1){ sep.end <- NULL} else {
            sep.end <- tan.timing[2:(length(tan.timing))]-1}
          if(length(sep.begin) - length(sep.end) == 1){ sep.end <- c(sep.end,dim(tandem)[1])}
          for(j in 1:length(sep.begin)){
            for(k in 1:length(sep.end)){
              if(sep.end[k]>sep.begin[j]){
                scheme[sep.begin[j]:sep.end[k]] <- "s"
                break;
              }
            }
          }
        }
      }
    }
    
    
    ## Plots
    if(Plot){
    png(file.path(PROJHOME, "plot/tandem/", paste0(df.temp$name[1], ".png")))
    par(pin=c(6,3))
    plot(ind.dis, ylim=c(0,15), type="l", main=ind.names[i])
    points(sl.follower, type="l", col=alpha("blue", 0.5))
    points(sl.leader, type="l", col=alpha("red", 0.5))
    points(tandem*2, type="l", col="green")
    dev.off()
    }
    
    ## Output
    ## tandem survival time
    {
      Tan.time = Cens = fem.speed = f.m.speed <- NULL
      if(sum(scheme=="t")>0){
        tan.timing <- which(tandem)[c(T, diff(which(tandem))>1)]
        tan.end <- which(tandem)[c(diff(which(tandem))>1,T)]
        for(j in 1:length(tan.timing)){
          Tan.time = c(Tan.time, tan.end[j] - tan.timing[j] + 1)
          f.speed = sl.follower[(tan.timing[j]+1):tan.end[j]] * 5
          fem.speed =c(fem.speed, mean(f.speed))

          if(tan.timing[j]==1 || tan.end[j]==length(tandem)){
            Cens = c(Cens, 0)
          } else {
            Cens = c(Cens, 1)
          }
        } 
        res.tan.time <- rbind(res.tan.time, data.frame(Video = df.temp$name[1],
                                                       Treat = df.temp$treat[1],
                                                       Colony = df.temp$colony[1],
                                                       Tan.time= Tan.time/5, Cens, 
                                                       F.speed = fem.speed))
      }
    }
    
    ## sep search survival time
    {
      if(sum(scheme=="s")>0){
        Sep.time = Cens =Sep.dis <- NULL
        sep.timing <- which(scheme=="s")[c(T, diff(which(scheme=="s"))>1)]
        sep.end <- which(scheme=="s")[c(diff(which(scheme=="s"))>1,T)]
        for(j in 1:length(sep.timing)){
          Sep.time = c(Sep.time, sep.end[j] - sep.timing[j] + 1)
          if(sep.timing[j]==1 || sep.end[j]==length(tandem)){
            Cens = c(Cens, 0)
          } else {
            Cens = c(Cens, 1)
          }
          Sep.dis = c(Sep.dis, mean(ind.dis[sep.timing[j]:sep.end[j]]))
        }
        res.sep.time <- rbind(res.sep.time, data.frame(Video = df.temp$name[1],
                                                       Treat = df.temp$treat[1],
                                                       Colony = df.temp$colony[1],
                                                       Sep.time=Sep.time/5, Cens,
                                                       Sep.dis = Sep.dis))
      }
    }
    
    ## summarized individual tandem data
    res.tandem = rbind( res.tandem, data.frame(
      Video = df.temp$name[1],
      Treat = df.temp$treat[1],
      Colony = df.temp$colony[1],
      Tandem.Duration = sum(tandem)/5,
      Tandem.Proportion = sum(tandem)/length(tandem)
    ))
    
    ## all data
    res.all = rbind(res.all,
                    data.frame(
                      df.follower[, c(1:2, 4:5, 7:8)],
                      x.follower = df.follower$x,
                      y.follower = df.follower$y,
                      x.leader = df.leader$x,
                      y.leader = df.leader$y,
                      sl.follower = df.follower$sl,
                      sl.leader = df.leader$sl,
                      acc.follower = df.follower$acc,
                      acc.leader = df.leader$acc,
                      ind.dis = sqrt( (x.follower - x.leader)^2 + (y.follower - y.leader)^2 ),
                      scheme,
                      polarization= c(polarization,NA),
                      rotation = c(rotation,NA)
                    )
                    )
    
  }
}


## 1. Survival analysis for tandem duration and sep duration
{
  m <- coxme(Surv(Tan.time, Cens) ~ Treat + (1|Colony/Video), data = res.tan.time)
  summary(m)
  Anova(m)
  
  ## ggplots
  ggsurvplot(
    fit = survfit(Surv(Tan.time, Cens) ~ Treat, 
                  type = "kaplan-meier", 
                  data = res.tan.time),
    conf.int = T,
    xlab = "Duration (sec)", 
    ylab = "Tandem Prob",
    xlim = c(0,300),
    palette = viridis(3)[c(1:3)],
    legend = c(0.8,0.8)
    
  )
  
  res.sep.time$Treat <-  as.factor(res.sep.time$Treat)
  m <- coxme(Surv(Sep.time, Cens) ~ Treat + (1|Colony/Video), data = res.sep.time)
  summary(m)
  Anova(m)
  
  ## ggplots
  ggsurvplot(
    fit = survfit(Surv(Sep.time, Cens) ~ Treat, 
                  type = "kaplan-meier", 
                  data = res.sep.time),
    conf.int = T,
    xlab = "Duration (sec)", 
    ylab = "Tandem Prob",
    xlim = c(0,60),
    palette = viridis(3)[c(1:3)],
    legend = c(0.8,0.8)
    
  )
  
  ggplot( res.tandem, aes(x=Treat, y= Tandem.Proportion)) +
    #geom_boxplot()
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, binwidth=.01)+
    geom_point(alpha = 0.4, position = position_jitter(.05), size = 2) +
    geom_boxplot(outlier.shape = NA, alpha=0.3, width =.2) +
    guides(fill = F)+
    scale_fill_viridis(discrete=T)+
    ylim(c(0,1.01))+
    theme_bw()
  
  
  ggplot( res.tandem, aes(x=Treat, y= Tandem.Proportion, fill=Treat)) +
    geom_violin(trim=F, alpha = 0.4) +
    geom_boxplot(outlier.shape = NA, alpha=0.3, width =.2) +
    guides(fill = F)+
    scale_fill_viridis(discrete=T)+
    theme_bw()
  
  
  r <- lmer(Tandem.Proportion ~ Treat + (1|Colony), data=res.tandem)
  Anova(r)
  
  
  res.tandem$Treat <-  as.factor(res.tandem$Treat)
  m <- coxme(Surv(Tandem.Duration) ~ Treat + (1|Colony), data = res.tandem)
  summary(m)
  Anova(m)
  
  ## ggplots
  ggsurvplot(
    fit = survfit(Surv(Tandem.Duration) ~ Treat, 
                  type = "kaplan-meier", 
                  data = res.tandem),
    conf.int = T,
    xlab = "Duration (sec)", 
    ylab = "Tandem Prob",
    xlim = c(0,300),
    palette = viridis(3)[c(1:3)],
    legend = c(0.8,0.8)
    
  )
}

## 2. Correalation between acceleration and inter-individual distance
df.temp <- data.frame(
  res.all[, c(2,3,4,6)],
  role = rep(c("leader", 'follower'), each=dim(res.all)[1]),
  scheme = res.all$scheme,
  acc = c(res.all$acc.leader, res.all$acc.follower),
  ind.dis = res.all$ind.dis)

df.temp <- df.temp[sample(1:dim(df.temp)[1], 60000, replace = F),]
ggplot(df.temp[df.temp$scheme == 't',], aes(x = ind.dis, y=acc, col=role))+
  geom_point(alpha=0.1)+
  facet_wrap(vars(treat)) + 
  ylim(c(-20,20)) +
  stat_smooth(se =T, method = "lm")

dim(df.temp)

  

##

library(ggridges)
ggplot(res.all[res.all$scheme=="t",], aes(x=polarization, y=treat)) +
  geom_density_ridges(fill=2, alpha=0.2, draw_baseline = T)
ggplot(res.all[res.all$scheme=="t",], aes(x=rotation, y=treat)) +
  geom_density_ridges(fill=2, alpha=0.2, draw_baseline = T)
ggplot(res.all[res.all$scheme=="t",], aes(x=ind.dis, y=treat)) +
  geom_density_ridges(fill=2, alpha=0.2, draw_baseline = T) +
  geom_vline(xintercept  = 4)

  
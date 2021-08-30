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
                      scheme
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

  
  
  
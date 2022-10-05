## R. speratus homosexual tandem analysis
## N. Mizumoto

#------------------------------------------------------------------------------#
# This file is for sensitivity analysis
#------------------------------------------------------------------------------#
rm(list = ls())
{
  source("scripts/Source.R")
  Sensitivity.analysis.Tandem()
  Sensitivity.analysis.Competition()
  Sensitivity.analysis.MovementParameter(time.window = 30)
  Sensitivity.analysis.MovementParameter(time.window = 90)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Sensitiviety analysis for the thresholds of tandem run
#------------------------------------------------------------------------------#
Sensitivity.analysis.Tandem <- function(){
  Tandem.SA <- function(interaction.threshold, min.tandem.sec, threshold.moving.dis.tandem){
    load(paste0(rda.place, "/AllData-scaled.rda"))
    
    ## get data only during tandem (not separation search)
    df.tan <- df[df$scheme != "sep",]
    
    ## calcurate step length
    df.tan$sl <- c(NA, sqrt( diff(df.tan$x)^2 + diff(df.tan$y)^2))
    df.tan$sl[df.tan$time == 0] <- NA
    
    ## main calculation
    ind.names <- unique(df.tan$name) 
    df.tandem.sum <- NULL
    for(i in 1:length(ind.names)){
      print(paste(i, "/", length(ind.names), "->", ind.names[i]))
      
      {
        df.temp <- df.tan[df.tan$name == ind.names[i],]
        df.temp <- na.omit(df.temp)
        df.follower <- df.temp[df.temp$role == "Follower",]
        df.leader <- df.temp[df.temp$role == "Leader",]
        #if(dim(df.temp)[1] < 2998){print(paste("not enough data", dim(df.temp)[1], "next"));next;}
        x.follower = df.follower$x
        y.follower = df.follower$y
        x.leader = df.leader$x
        y.leader = df.leader$y
        sl.follower = df.follower$sl
        sl.leader = df.leader$sl
        ind.dis <- sqrt( (x.follower - x.leader)^2 + (y.follower - y.leader)^2 )
      }
      
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
        # determine leader or follower based on the relative position:
        # assuming that follower is at the behind of leader's heading direction
        # leader is at the front of follower's heading direction
        lv <- length(x.leader)
        angle.leader = atan2(y.leader[2:lv] - y.leader[2:lv-1] , x.leader[2:lv] - x.leader[2:lv-1])
        angle.follower = atan2(y.follower[2:lv] - y.follower[2:lv-1] , x.follower[2:lv] - x.follower[2:lv-1])
        angle.l.to.f = atan2(y.follower - y.leader , x.follower - x.leader)[2:lv-1]
        angle.f.to.l = atan2(y.leader - y.follower , x.leader - x.follower)[2:lv-1]
        r.dir.l.to.f <- apply(cbind(abs(angle.leader - angle.l.to.f), pi*2-abs(angle.leader - angle.l.to.f)), 1, min)
        r.dir.f.to.l <- apply(cbind(abs(angle.follower - angle.f.to.l), pi*2-abs(angle.follower - angle.f.to.l)), 1, min)
        
        # check if the estimated role is consistent with the pre-allocated leader/follower role
        correct.role <- r.dir.f.to.l < pi/2 & abs(pi-r.dir.l.to.f) < pi/2
        reverse.role <- r.dir.l.to.f < pi/2 & abs(pi-r.dir.f.to.l) < pi/2
        correct.role <- (c(correct.role,F) & tandem)
        reverse.role <- (c(reverse.role,F) & tandem)
        
      }
      
      ## definition of competition over follower position
      {
        # 1. the distance between individuals needed to be smaller than 4mm
        competition.dis <- ind.dis < 4
        
        # 2. the rotation index of a pair needed to be larger than 0.5 
        {
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
          
          rotation <- sqrt(
            (x.leader.ang.moment + x.follower.ang.moment)^2+
              (y.leader.ang.moment + y.follower.ang.moment)^2
          )/2
          
          polarization <- sqrt(
            (x.leader.vec + x.follower.vec)^2+
              (y.leader.vec + y.follower.vec)^2
          )/2
          
          rotation[is.na(rotation)] <- 0
        }
        
        # 3. the above two conditions needed to persist for more than 2 seconds
        {
          competition <- tandem & competition.dis & c(rotation, F) > 0.5
          competition <- tandem.smoothing(competition, !competition, min.tandem.sec*5)
        }
      }
      
      ## Output
      ## summarized individual tandem data
      df.tandem.sum = rbind( df.tandem.sum, data.frame(
        Video = df.temp$name[1],
        Treat = df.temp$treat[1],
        Colony = df.temp$colony[1],
        Tandem.Duration = sum(tandem&!competition)/5,
        Tandem.Proportion = sum(tandem&!competition)/length(tandem),
        Competition = sum(competition)/5,
        Competition.Proportion = sum(competition)/length(tandem),
        TandemComp.Proportion = sum(tandem)/length(tandem)
      ))
    }
    return(df.tandem.sum)
  }
  
  Tandem.SA.plot <- function(param = c(2,2,2)){
    print(paste("interaction dis:", p_interaction[param[1]],
                "minimum tandem time:", p_min.tandem[param[2]],
                "moved distance required:", p_moved.dis[param[3]]))
    
    df_SAtemp <- Tandem.SA(p_interaction[param[1]], p_min.tandem[param[2]], p_moved.dis[param[3]])
    df_SAtemp$Treat <- factor(df_SAtemp$Treat, levels = c("FM", "FF", "MM"))
    
    e = 0.01
    y <- df_SAtemp$TandemComp.Proportion
    df_SAtemp$logit.TandemComp.Proportion <- (log((y+e)/(1-y+e)))
    r <- lmer(logit.TandemComp.Proportion ~ Treat + (1|Colony), data=df_SAtemp)
    res <- Anova(r)
    
    fig_label <- paste0(p_interaction[param[1]], 
                        "mm, ",  p_min.tandem[param[2]], 
                        "sec, ", p_moved.dis[param[3]], 
                        "mm; LMM, P = ", round(res$`Pr(>Chisq)`,3))
    
    ggplot(df_SAtemp, aes(x=Treat, y=TandemComp.Proportion))+
      geom_dotplot(binaxis = "y", stackdir='center', dotsize=1, binwidth = 0.01)+
      ylim(c(0,1))+
      stat_summary(fun.data=mean_se, 
                   geom="pointrange", color="red") +
      theme_bw() + ylab("Prop of time engaging in tandems") +
      xlab("")+
      theme(legend.position = "none", aspect.ratio = 0.75) +
      ggtitle(fig_label)
    
    file_label = paste0(p_interaction[param[1]], 
           "mm",  p_min.tandem[param[2]], 
           "sec", p_moved.dis[param[3]], 
           "mm")
    
    ggsave(filename = file.path("plot/", paste0("SA_",file_label, ".pdf")),
           width=4, height = 3, family="PT Sans")
  }
  
  p_interaction <- c(6.5, 7, 7.5)
  p_min.tandem  <- c(1, 2, 3)
  p_moved.dis   <- c(20, 30, 40)
  
  Tandem.SA.plot(c(2,2,2))
  Tandem.SA.plot(c(1,2,2))
  Tandem.SA.plot(c(3,2,2))
  Tandem.SA.plot(c(2,1,2))
  Tandem.SA.plot(c(2,3,2))
  Tandem.SA.plot(c(2,2,1))
  Tandem.SA.plot(c(2,2,3))
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Sensitiviety analysis for the thresholds of competition
#------------------------------------------------------------------------------#
Sensitivity.analysis.Competition <- function(){
  compeittion_dis_thresh = 4
  rotation_index_thresh = 0.5
  competition_time_thresh = 2
  
  p_comp.dis = c(3.5, 4, 4.5)
  p_rot.ind = c(0.4, 0.5, 0.6)
  p_comp.time = c(1.5, 2, 2.5)
  
  Competition.SA <- function(compeittion_dis_thresh, rotation_index_thresh, competition_time_thresh){
    load(paste0(rda.place, "/AllData-scaled.rda"))
    ## get data only during tandem (not separation search)
    df.tan <- df[df$scheme != "sep",]
    
    ## calcurate step length
    df.tan$sl <- c(NA, sqrt( diff(df.tan$x)^2 + diff(df.tan$y)^2))
    df.tan$sl[df.tan$time == 0] <- NA
    
    ## main calculation
    ind.names <- unique(df.tan$name) 
    df.tandem.sum = df.tandem = df.separate.duration = df.tandem.duration <- NULL
    for(i in 1:length(ind.names)){
      {
        df.temp <- df.tan[df.tan$name == ind.names[i],]
        df.temp <- na.omit(df.temp)
        df.follower <- df.temp[df.temp$role == "Follower",]
        df.leader <- df.temp[df.temp$role == "Leader",]
        #if(dim(df.temp)[1] < 2998){print(paste("not enough data", dim(df.temp)[1], "next"));next;}
        x.follower = df.follower$x
        y.follower = df.follower$y
        x.leader = df.leader$x
        y.leader = df.leader$y
        sl.follower = df.follower$sl
        sl.leader = df.leader$sl
        ind.dis <- sqrt( (x.follower - x.leader)^2 + (y.follower - y.leader)^2 )
      }
      
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
        # determine leader or follower based on the relative position:
        # assuming that follower is at the behind of leader's heading direction
        # leader is at the front of follower's heading direction
        lv <- length(x.leader)
        angle.leader = atan2(y.leader[2:lv] - y.leader[2:lv-1] , x.leader[2:lv] - x.leader[2:lv-1])
        angle.follower = atan2(y.follower[2:lv] - y.follower[2:lv-1] , x.follower[2:lv] - x.follower[2:lv-1])
        angle.l.to.f = atan2(y.follower - y.leader , x.follower - x.leader)[2:lv-1]
        angle.f.to.l = atan2(y.leader - y.follower , x.leader - x.follower)[2:lv-1]
        r.dir.l.to.f <- apply(cbind(abs(angle.leader - angle.l.to.f), pi*2-abs(angle.leader - angle.l.to.f)), 1, min)
        r.dir.f.to.l <- apply(cbind(abs(angle.follower - angle.f.to.l), pi*2-abs(angle.follower - angle.f.to.l)), 1, min)
        
        # check if the estimated role is consistent with the pre-allocated leader/follower role
        correct.role <- r.dir.f.to.l < pi/2 & abs(pi-r.dir.l.to.f) < pi/2
        reverse.role <- r.dir.l.to.f < pi/2 & abs(pi-r.dir.f.to.l) < pi/2
        correct.role <- (c(correct.role,F) & tandem)
        reverse.role <- (c(reverse.role,F) & tandem)
        
      }
      
      ## definition of competition over follower position
      {
        # 1. the distance between individuals needed to be smaller than 4mm
        competition.dis <- ind.dis < compeittion_dis_thresh
        
        # 2. the rotation index of a pair needed to be larger than 0.5 
        {
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
          
          rotation <- sqrt(
            (x.leader.ang.moment + x.follower.ang.moment)^2+
              (y.leader.ang.moment + y.follower.ang.moment)^2
          )/2
          
          polarization <- sqrt(
            (x.leader.vec + x.follower.vec)^2+
              (y.leader.vec + y.follower.vec)^2
          )/2
          
          rotation[is.na(rotation)] <- 0
        }
        
        # 3. the above two conditions needed to persist for more than 2 seconds
        {
          competition <- tandem & competition.dis & c(rotation, F) > rotation_index_thresh
          competition <- tandem.smoothing(competition, !competition, competition_time_thresh*5)
        }
      }
      
      ## Output
      ## summarized individual tandem data
      df.tandem.sum = rbind( df.tandem.sum, data.frame(
        Video = df.temp$name[1],
        Treat = df.temp$treat[1],
        Colony = df.temp$colony[1],
        Tandem.Duration = sum(tandem&!competition)/5,
        Tandem.Proportion = sum(tandem&!competition)/length(tandem),
        Competition = sum(competition)/5,
        Competition.Proportion = sum(competition)/length(tandem),
        TandemComp.Proportion = sum(tandem)/length(tandem)
      ))
    }
    return(df.tandem.sum)
  }
  
  
  Competition.SA.plot <- function(param = c(2,2,2)){
    print(paste("compeition dis:", p_comp.dis[param[1]],
                "rotation index:", p_rot.ind[param[2]],
                "competition time:", p_comp.time[param[3]]))
    df_SAtemp <- Competition.SA(p_comp.dis[param[1]], p_rot.ind[param[2]], p_comp.time[param[3]])
    df_SAtemp$Treat <- factor(df_SAtemp$Treat, levels = c("FM", "FF", "MM"))
    return(table(subset(df_SAtemp, Competition>5)$Treat))
  }
  print(Competition.SA.plot(c(2,2,2)))
  print(Competition.SA.plot(c(1,2,2)))
  print(Competition.SA.plot(c(3,2,2)))
  print(Competition.SA.plot(c(2,1,2)))
  print(Competition.SA.plot(c(2,3,2)))
  print(Competition.SA.plot(c(2,2,1)))
  print(Competition.SA.plot(c(2,2,3)))
  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Sensitiviety analysis for time time window to obtain movemen parameters
#------------------------------------------------------------------------------#
Sensitivity.analysis.MovementParameter <- function(time.window = 60){
  print(paste("focused time window = ", time.window))
  load(paste0(rda.place, "/AllData-scaled.rda"))
  
  ## Obtain moved distance for every second as a proxy of movement speed (df.sec)
  {
    ## calcurate step length
    df <- na.omit(df)
    df$sl <- c(NA, sqrt( diff(df$x)^2 + diff(df$y)^2))
    df$sl[df$time == 0] <- NA
    
    df$angle <- c(NA, angle_cal(df$x, df$y, dim(df)[1]), NA)
    df$angle[df$time == 0] <- NA
    df$angle[df$time == 300] <- NA
    
    ind.names <- unique(df$name) 
    scheme.names <- unique(df$scheme) 
    
    ## step length was calculated in 5FPS
    ## then sum up every second to obtain mm/sec
    df.sec <- data.frame()
    for(i in 1:length(ind.names)){
      for(j in 1:3){
        df.temp <- df[df$name == ind.names[i] & df$scheme == scheme.names[j],]
        if(scheme.names[j] == "tan-L" && df.temp$role[1] == "Follower"){ next; }
        if(scheme.names[j] == "tan-F" && df.temp$role[1] == "Leader"){ next; }
        
        if(dim(df.temp)[1] < 1801){
          df.temp <- df.temp[df.temp$time <= floor(max(df.temp$time)),]
        }
        
        df.temp <- df.temp[2:dim(df.temp)[1],]
        df.temp.out <- df.temp[df.temp$time %% 1==0,1:8]
        df.temp.out$speed <- df.temp$sl[(df.temp$time*10) %% 10 == 2] + 
          df.temp$sl[(df.temp$time*10) %% 10 == 4] +
          df.temp$sl[(df.temp$time*10) %% 10 == 4] +
          df.temp$sl[(df.temp$time*10) %% 10 == 8] +
          df.temp$sl[(df.temp$time*10) %% 10 == 0]
        df.sec <- rbind(df.sec, df.temp.out)
      }
    }
    df.sec$scheme[df.sec$scheme != "sep"] <- "tan"
  }
  
  ## Change the focused time window and obtain movement parameters
  {
    {
      df.sec.sep <- df.sec[df.sec$time < (time.window+1) & df.sec$scheme =="sep",]
      df.sep.sec.sum <- 
        data.frame(
          df.sec.sep[df.sec.sep$time == 1, 1:7],
          speed.mean = tapply(df.sec.sep$speed, df.sec.sep$name, mean)[df.sec.sep[df.sec.sep$time == 1,]$name]
        )
      
      df.sep.sec.sum$treat <- factor(df.sep.sec.sum$treat, levels = c("FM", "FF", "MM"))
    }
    {
      df$angle <- c(NA, angle_cal(df$x, df$y, dim(df)[1]), NA)
      df$angle[df$time == 0] <- NA
      df$angle[df$time == 300] <- NA
      
      df.sep.focus <- df[df$scheme=="sep" & df$time <= time.window,]
      df.sep.turn <- df.sep.focus[df.sep.focus$time==0, 1:7]
      
      ind.names <- unique(df$name) 
      mu = rho = c()
      for(i in 1:length(ind.names)){
        df.temp <- df.sep.focus[df.sep.focus$name == ind.names[i],]
        mle_wrpcauchy <- wrpcauchy.ml(na.omit(df.temp$angle), 0, 0, acc=1e-015)
        mu <- c(mu, as.numeric(mle_wrpcauchy[1]))
        rho <- c(rho, as.numeric(mle_wrpcauchy[2]))
      }
      df.sep.turn$mu = mu
      df.sep.turn$rho = rho
    }
  }
  
  ## Obtain parameters for simulations
  {
    df.temp <- data.frame(
      treat = rep(c("FF", "FM", "MM"), 2),
      role = rep( c("F", "L"), each=3),
      speed = round(as.vector(tapply(df.sep.sec.sum$speed.mean, df.sep.sec.sum[,c("treat", "role")], mean)), 3),
      sinuosity =round( as.vector(tapply(df.sep.turn$rho, df.sep.turn[,c("treat", "role")], mean)), 3)
    )
    print(df.temp)
  }
}
#------------------------------------------------------------------------------#

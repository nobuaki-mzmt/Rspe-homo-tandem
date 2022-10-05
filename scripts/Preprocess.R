## R. speratus homosexual tandem analysis
## N. Mizumoto

#------------------------------------------------------------------------------#
# This file is for preprocess all data for movement analysis (tandem and separation search)
#------------------------------------------------------------------------------#
rm(list = ls())
{
  source("scripts/Source.R")
  convert.rda()
  scale.data()
  Tandem.analysis()
  Movement.analysis()
}

#------------------------------------------------------------------------------#
# This function reads all csv files, then
# 1) convert them into R data.frame and saved as .rda files
# 2) output raw trajectories in .png file
#------------------------------------------------------------------------------#
convert.rda <- function(Plot = F,      # See overall trajectories
                        Dataframe = T # Integrate all in one dataframe 
                        ){
  ## data
  rawdata <- list.files("data/raw_movement", full.names = TRUE, pattern = ".csv")
  dataname <- list.files("data/raw_movement", full.names = F, pattern = ".csv")
  
  ## options
  options(warn = 0)
  
  df <- data.frame()
  ## plots
  for(v in 1:length(rawdata)){
    
    # file info
    d <- data.frame(fread(rawdata[v], header=T))
    d <- d[1:(5*60*30+1),]
    
    species = substr(dataname[v], 1, 2)
    treat = substr(dataname[v], 4, 5)
    role = substr(dataname[v], 7, 7)
    if(role=="F"){role ="Follower"} else {role="Leader"}
    colony = substr(dataname[v], 9, 9)
    id = substr(dataname[v], 9, 10)
    scheme = substr(dataname[v], 12, 14) 
    name <- paste(species, treat, id, role, sep="-")
    print(paste(v, "/", length(rawdata), "->", name))
    
    # plot
    if(Plot){
      png(file.path("plot/raw-trajectories/", str_replace(dataname[v], "csv", "png")))
      par(mfrow=c(1,1), pin=c(3,2))
      plot( d[,2], d[,3], type="l", col=1, xlim=c(0,1800), ylim=c(0,1200),
            main=paste(name, scheme),
            xlab="x(mm)", ylab="y(mm)")
      dev.off()
    }
    
    # datafrmae
    d[,1] <- d[,1]/30
    d <- d[seq(1,9001,6),]
    if(Dataframe){
      if(scheme == "tan"){
        dftemp1 <- data.frame(species, treat, role="Leader", colony, id,
                              scheme = paste0(scheme, "-L"), name, time = d[,1], x=d[,2], y=d[,3])
        df <- rbind(df,dftemp1)
        dftemp2 <- data.frame(species, treat, role="Follower", colony, id,
                              scheme = paste0(scheme, "-F"), name, time = d[,1], x=d[,4], y=d[,5])
        df <- rbind(df,dftemp2)
      } else {
        dftemp1 <- data.frame(species, treat, role, colony, id,
                              scheme = scheme, name, time = d[,1], x=d[,2], y=d[,3])
        df <- rbind(df,dftemp1)
      }
    }
  }
  f.name <- paste0(rda.place, "/AllData.rda")
  save(df, file = f.name)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function reads AllData.rda, then
# 1) perform scaling:
#  assuming that termites moved throughout the arena, scale the movement range as arena size
#  scaling is performed for each observed arena (pool tandem and separation for each individual)
#  also check speed time development to detect tracking error
# 2) save scaled data as AllData-scaled.rda
#------------------------------------------------------------------------------#
scale.data <- function(Plot = F){
  dir.create("plot/scale-trajectories", showWarnings = FALSE)
  load(paste0(rda.place, "/AllData.rda"))
  ind.names <- unique(df$name) 
  
  for(v in 1:length(ind.names)){
    df.temp <- df[df$name == ind.names[v],]
    print(paste(v, "/", length(ind.names), "->", ind.names[v]))
    
    x <- df.temp$x
    y <- df.temp$y
    xL <- max(x, na.rm=T) - min(x, na.rm=T)
    yL <- max(y, na.rm=T) - min(y, na.rm=T)
    df.temp$x <- (x-min(x, na.rm=T))/xL * arena.size
    df.temp$y <- (y-min(y, na.rm=T))/yL * arena.size
    
    # plot
    if(Plot){
      p1 <- ggplot(data=df.temp, aes(x=x, y=y, col=as.factor(scheme), group=as.factor(scheme))) +
        geom_path() + 
        theme_bw() +
        theme(legend.position = "none", aspect.ratio = 1)+ 
        xlim(c(0,140)) + ylim(c(0,140)) +
        ggtitle(paste(v, ind.names[v]))
      df.temp$speed = c(NA, sqrt( diff(df.temp$x)^2 + diff(df.temp$y)^2 ) )
      p2 <- ggplot(data=df.temp, aes(x=time, y=speed, col=as.factor(scheme), group=as.factor(scheme))) +
        geom_path() + 
        theme_bw() +
        theme(legend.position = "none", aspect.ratio = 1)+ 
        coord_cartesian(xlim=c(0,300), ylim=c(0, 30)) +
        ggtitle(paste(v, ind.names[v]))
      ggsave(file.path("plot/scale-trajectories/", paste0(df.temp$name[1], ".png")), 
             plot = grid.arrange(p1, p2, nrow = 1), width = 6, height = 3)
    }
    
    # data.frame
    df[df$name == ind.names[v],] <- df.temp
  }  
  f.name <- paste0(rda.place, "/AllData-scaled.rda")
  save(df, file = f.name)
}
#------------------------------------------------------------------------------#

## Caution! threshold for competition was 2 sec in the paper, but calucuated by 1 sec in the code (fixed 220310 but need to confirm the results were the same)
#------------------------------------------------------------------------------#
# This function reads AllData-scaled.rda, then
# 1) Determine if a pair is doing tandem or not for each frame.
# 2) Determine if a pair is competiting over followe position for each frame.
# 3) save dataframe for each step (df_tandem.rda),
#    for summary for each pair (df_sum_tandem.rda),
#    and for duration of separation search (df_separate_duration.rda)
#------------------------------------------------------------------------------#
Tandem.analysis <- function(Plot = F){
  load(paste0(rda.place, "/AllData-scaled.rda"))
  
  ## get data only during tandem (not separation search)
  df.tan <- df[df$scheme != "sep",]
  
  ## calcurate step length
  df.tan$sl <- c(NA, sqrt( diff(df.tan$x)^2 + diff(df.tan$y)^2))
  df.tan$sl[df.tan$time == 0] <- NA
  
  ## calculate accerelation
  df.tan$acc <- c(diff(df.tan$sl*5), NA)
  df.tan$acc[df.tan$time == 0 | df.tan$time == 300] <- NA
  
  ## main calculation
  ind.names <- unique(df.tan$name) 
  df.tandem.sum = df.tandem = df.separate.duration = df.tandem.duration <- NULL
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
        
        if(Plot){
          png(file.path("plot/tandem-reverse/", paste0(df.temp$name[1], ".png")))
          par(mfrow=c(2,1), pin=c(6,2))
          plot((tandem)*0, col=(tandem)*2, type="p", ylim=c(-1,1), main=ind.names[i])
          points((correct.role - reverse.role), type="l")
          
          plot(polarization, col=2, type="l", ylim=c(0,1))
          points(rotation, col=4, type="l")
          
          dev.off()
        }
        
        rotation[is.na(rotation)] <- 0
      }
      
      # 3. the above two conditions needed to persist for more than 2 seconds
      {
        competition <- tandem & competition.dis & c(rotation, F) > 0.5
        competition <- tandem.smoothing(competition, !competition, min.tandem.sec*5)
      }
      
      if(Plot){
        par(mfrow=c(2,1), pin=c(6,2))
        plot((tandem)*0, col=(tandem)*2, type="p", ylim=c(-1,1), main=ind.names[i])
        points((correct.role - reverse.role), type="l")
        points(competition*0-0.2, col=competition*4)
        
        plot(polarization, col=2, type="l", ylim=c(0,1))
        points(rotation, col=4, type="l")
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
      png(file.path("plot/tandem/", paste0(df.temp$name[1], ".png")))
      par(pin=c(6,3))
      plot(ind.dis, ylim=c(0,15), type="l", main=ind.names[i])
      points(sl.follower, type="l", col=alpha("blue", 0.5))
      points(sl.leader, type="l", col=alpha("red", 0.5))
      points(tandem*2, type="l", col="green")
      dev.off()
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
    
    ## tandem duration
    {
      if(sum(scheme=="t")>0){
        tan.timing <- which(tandem)[c(T, diff(which(tandem))>1)]
        tan.end <- which(tandem)[c(diff(which(tandem))>1,T)]
        cens = rep(1, length(tan.end))
        cens[tan.timing[j]==1 || tan.end[j]==length(tandem)] = 0
        dftemp <- data.frame(
          Video = df.temp$name[1],
          Treat = df.temp$treat[1],
          Colony = df.temp$colony[1],
          Duration = (tan.end-tan.timing+1) / 5,
          Cens = cens
        )
        df.tandem.duration = rbind(df.tandem.duration, dftemp)
      }
    }

    ## separation duration
    {
      if(sum(scheme=="s")>0){
        sep.timing <- which(scheme=="s")[c(T, diff(which(scheme=="s"))>1)]
        sep.end <- which(scheme=="s")[c(diff(which(scheme=="s"))>1,T)]
        cens = rep(1, length(sep.end))
        cens[sep.timing[j]==1 || sep.end[j]==length(tandem)] = 0
        dftemp <- data.frame(
          Video = df.temp$name[1],
          Treat = df.temp$treat[1],
          Colony = df.temp$colony[1],
          Duration = (sep.end-sep.timing+1) / 5,
          Cens = cens
        )
        df.separate.duration = rbind(df.separate.duration, dftemp)
      }
    }
      
      
    ## all data
    df.tandem = rbind(df.tandem,
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
                        rotation = c(rotation,NA),
                        competition = competition
                      )
    )
  }
  
  
  f.name <- paste0(rda.place, "/df_tandem_sum.rda")
  save(df.tandem.sum, file = f.name)
  
  f.name <- paste0(rda.place, "/df_tandem.rda")
  save(df.tandem, file = f.name)
  
  f.name <- paste0(rda.place, "/df_tandem_duration.rda")
  save(df.tandem.duration, file = f.name)
  
  f.name <- paste0(rda.place, "/df_separation_duration.rda")
  save(df.separate.duration, file = f.name)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function reads AllData-scaled.rda, then
# 1) Obtain moved distance for every second as a proxy of movement speed 
# 2) Get mean movement speed across individuals for each role and sex
# 3) Get mean speed for last 1min for tandem runs and first 1min after separation
# 4) Fitting turning patterns to wrapped Cauchy distribution for each individual
#------------------------------------------------------------------------------#
Movement.analysis <- function(){
  load(paste0(rda.place, "/AllData-scaled.rda"))
  
  ## 1. Obtain moved distance for every second as a proxy of movement speed (df.sec)
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
      print(paste(i, "/", length(ind.names), "->", ind.names[i]))
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
  
  ## 2. Get mean movement speed across individuals for each role and sex
  ## for plotting purpose (df.sec.sum)
  {
    treat.names <- unique(df.sec$treat)
    scheme.names <- unique(df.sec$scheme)
    df.sec.sum <- data.frame()
    for(i in 1:2){
      for(j in 1:3){
        df.sec.temp <- df.sec[df.sec$scheme ==scheme.names[i] &
                                df.sec$treat == treat.names[j], ]
        ymean <- tapply(df.sec.temp$speed, df.sec.temp[,c("time","role")], mean)
        yse <- tapply(df.sec.temp$speed, df.sec.temp[,c("time","role")], se)
        df.temp1 <- data.frame(treat=treat.names[j], scheme = scheme.names[i],
                               role = "Follower", time = 1:300,
                               speed.mean = ymean[,1], speed.se = yse[,1],
                               speed.mean.m.se = ymean[,1] - yse[,1],
                               speed.mean.p.se = ymean[,1] + yse[,1])
        df.temp2 <- data.frame(treat=treat.names[j], scheme = scheme.names[i],
                               role = "Leader", time = 1:300,
                               speed.mean = ymean[,2], speed.se = yse[,2],
                               speed.mean.m.se = ymean[,2] - yse[,2],
                               speed.mean.p.se = ymean[,2] + yse[,2])
        df.sec.sum <- rbind(df.sec.sum, df.temp1)
        df.sec.sum <- rbind(df.sec.sum, df.temp2)
      }
    }
    df.sec.sum$scheme <- factor(df.sec.sum$scheme, levels=c("tan", "sep"))
    df.sec.sum$treat <- factor(df.sec.sum$treat, levels=c("FM", "FF", "MM"))
    df.sec.sum[df.sec.sum$treat=="FM" & df.sec.sum$role == "Follower",]$role <- "M-Follower"
    df.sec.sum[df.sec.sum$treat=="FM" & df.sec.sum$role == "Leader",]$role <- "F-Leader"
    df.sec.sum[df.sec.sum$treat=="FF",]$role <- paste0("F-", df.sec.sum[df.sec.sum$treat=="FF",]$role)
    df.sec.sum[df.sec.sum$treat=="MM",]$role <- paste0("M-", df.sec.sum[df.sec.sum$treat=="MM",]$role)
    df.sec.sum$role <- factor(df.sec.sum$role, levels=c("F-Leader", "F-Follower",
                                                        "M-Leader", "M-Follower"))
    f.name <- paste0(rda.place, "/df_sec_sum.rda")
    save(df.sec.sum, file = f.name)
  }
  
  ## 3. Focus on the last 1min for tandem runs and first 1min after separation
  ## for statistical analysis (df.sep.sec.sum, df.tan.sec.sum)
  {
    # first 1 min after separation
    {
      df.sec.sep <- df.sec[df.sec$time < 61 & df.sec$scheme =="sep",]
      df.sep.sec.sum <- 
        data.frame(
          df.sec.sep[df.sec.sep$time == 1, 1:7],
          speed.mean = tapply(df.sec.sep$speed, df.sec.sep$name, mean)[df.sec.sep[df.sec.sep$time == 1,]$name]
        )
      
      df.sep.sec.sum$treat <- factor(df.sep.sec.sum$treat, levels = c("FM", "FF", "MM"))
    }
    
    
    # last 1 min during tandem
    {
      df.sec.tan <- df.sec[df.sec$time > 240 & df.sec$scheme =="tan",]
      df.tan.sec.sum <- 
        data.frame(
          df.sec.tan[df.sec.tan$time == 241, 1:7],
          speed.mean = tapply(df.sec.tan$speed, df.sec.tan$name, mean)[df.sec.tan[df.sec.tan$time == 241,]$name]
        )
      
      df.rem.partner <- df.tan.sec.sum$role ==  str_extract(df.tan.sec.sum$name, "Leader")
      df.rem.partner[is.na(df.rem.partner)] <- na.omit(df.tan.sec.sum$role ==  str_extract(df.tan.sec.sum$name, "Follower"))
      df.tan.sec.sum <- df.tan.sec.sum[df.rem.partner,]
      
      df.tan.sec.sum$treat <- factor(df.tan.sec.sum$treat, levels = c("FM", "FF", "MM"))
    }
    
    save(df.sep.sec.sum, file = paste0(rda.place, "/df_sec_sum_sep.rda"))
    save(df.tan.sec.sum, file = paste0(rda.place, "/df_sec_sum_tan.rda"))
  }
  
  ## 4. Fitting turning patterns to wrapped Cauchy distribution for each individual
  ## save df.sep.turn
  {
    df$angle <- c(NA, angle_cal(df$x, df$y, dim(df)[1]), NA)
    df$angle[df$time == 0] <- NA
    df$angle[df$time == 300] <- NA
    
    df.sep.focus <- df[df$scheme=="sep" & df$time <= 60,]
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
    save(df.sep.turn, file = paste0(rda.place, "/df_sep_turn.rda"))
  }  
  
}
#------------------------------------------------------------------------------#

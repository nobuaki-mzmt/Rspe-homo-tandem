## R. speratus homosexual tandem analysis
## N. Mizumoto
## 7/12/2021 ~
## script with R projects

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
}

# Computation --------------------------------------------------------------
{
  ## calcurate step length
  df <- na.omit(df)
  df$sl <- c(NA, sqrt( diff(df$x)^2 + diff(df$y)^2))
  df$sl[df$time == 0] <- NA
  
  ind.names <- unique(df$name) 
  scheme.names <- unique(df$scheme) 
  
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
      df.temp.out <- df.temp[df.temp$time %% 1==0,1:7]
      df.temp.out$speed <- df.temp$sl[(df.temp$time*10) %% 10 == 2] + 
                            df.temp$sl[(df.temp$time*10) %% 10 == 4] +
                            df.temp$sl[(df.temp$time*10) %% 10 == 4] +
                            df.temp$sl[(df.temp$time*10) %% 10 == 8] +
                            df.temp$sl[(df.temp$time*10) %% 10 == 0]
      df.sec <- rbind(df.sec, df.temp.out)
    }
  }
  df.sec$scheme[df.sec$scheme != "sep"] <- "tan"
  
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
  ggplot(df.sec.sum, aes(x=time, y=speed.mean, color=role, fill=role)) +
    geom_ribbon(aes(ymin=speed.mean.m.se, ymax=speed.mean.p.se),alpha=0.5) +
    facet_grid(treat~scheme) +
    scale_color_viridis(discrete = T, end=0.5) +
    scale_fill_viridis(discrete = T, end=0.5) +
    theme_bw() + theme(aspect.ratio = 0.75) +
    scale_y_continuous(limits = c(0,20)) +
    xlab("Time (sec)") + ylab("Speed (mm/sec)")
}

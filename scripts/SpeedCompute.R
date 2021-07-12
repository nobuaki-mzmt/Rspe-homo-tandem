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
  
  today <- Sys.Date()
  
  ## Packages
  {
  }
  
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
  
  ggplot(df.sec, aes(x=time, y=speed, col=role)) +
    #geom_path() +
    stat_smooth() +
    facet_grid(scheme~treat)
  
  ind.names <- unique(df$name) 
  df.temp <- df.sep[df.sep$name == ind.names[v],]
  df.temp
}

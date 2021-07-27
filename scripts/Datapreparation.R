## R. speratus homosexual tandem analysis
## N. Mizumoto
## 7/9/2021 ~
## script with R projects

## Data preparation

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
  
  ## parameters
  fps <- 5
  arena.size <- 85 # mm
  
  ## file location
  dir.create(file.path(PROJHOME, "/data/", "rda"), showWarnings = FALSE)
  rda.place <- file.path(PROJHOME, "/data/", "rda")
}

## Raw data use
# reads all csv files and 1) output raw trajectories in .png file
# 2) convert them in 5FPS and R data.frame saved in .rda files
{
  ## data
  dir.create(file.path(PROJHOME, "/plot"), showWarnings = FALSE)
  dir.create(file.path(PROJHOME, "/plot/raw-trajectories"), showWarnings = FALSE)
  rawdata <- list.files(file.path(PROJHOME, "data/raw"), full.names = TRUE, pattern = ".csv")
  dataname <- list.files(file.path(PROJHOME, "data/raw"), full.names = F, pattern = ".csv")
  
  ## options
  options(warn = 0)
  Plot <- F      # See overall trajectories
  Dataframe <- T # Integrate all in one dataframe
  
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
      png(file.path(PROJHOME, "plot/raw-trajectories/", str_replace(dataname[v], "csv", "png")))
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

## Scaling
# assuming that termites moved throughout the arena, scale the movement range as arena size
# scaling is performed for each observed arena (pool tandem and separation for each individual)
# also check speed time development to detect tracking error
# save scaled data as .rda
{
  Plot <- F
  dir.create(file.path(PROJHOME, "/plot/scale-trajectories"), showWarnings = FALSE)
  
  f.name <- paste0(rda.place, "/AllData.rda")
  load(f.name)
  
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
      ggsave(file.path(PROJHOME, "plot/scale-trajectories/", paste0(df.temp$name[1], ".png")), 
             plot = grid.arrange(p1, p2, nrow = 1), width = 6, height = 3)
    }
    
    # data.frame
    df[df$name == ind.names[v],] <- df.temp
  }  
  f.name <- paste0(rda.place, "/AllData-scaled.rda")
  save(df, file = f.name)
}

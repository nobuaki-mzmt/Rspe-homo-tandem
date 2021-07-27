## R. speratus homosexual tandem analysis
## N. Mizumoto
## 7/27/2021 ~
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
  data.place <- file.path(PROJHOME, "/data/sim")
  f.name <- paste0(data.place, "/SimRes_rep100000_sec60.csv")
}

# Plot --------------------------------------------------------------------
{
  d <- data.frame(fread(f.name, header=T))
  matplot(seq(0.2,60,0.2), t(d[,3:dim(d)[2]]), type="l")
  
  df <- data.frame(
    Treat = paste(d$Leader, d$Follower),
    prop60 = d$X300/100000
  )
  
  ggplot(df, aes(x=Treat, y=prop60)) +
    geom_bar(stat="identity")
}

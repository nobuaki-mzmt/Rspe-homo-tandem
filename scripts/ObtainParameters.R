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
  df$angle <- c(NA, angle_cal(df$x, df$y, dim(df)[1]), NA)
  df$angle[df$time == 0] <- NA
  df$angle[df$time == 300] <- NA
  
  df.sep.focus <- df[df$scheme=="sep" & df$time <= 60,]
  
  df.param <- df.sep.focus[df.sep.focus$time==0, 1:6]
  speed = tapply(df.sep.focus$sl, df.sep.focus$name, mean, na.rm=T)
  df.param <- data.frame(df.param, speed = speed[df.param$name])
  
  tapply(df.param$speed, df.param[,c("treat", "role")], mean)
  tapply(df.param$speed, df.param[,c("treat", "role")], se)
  as.vector(tapply(df.param$speed, df.param[,c("treat", "role")], mean))
  
  wilcox.exact(speed ~ role, data=df.param[df.param$treat == "FF",], paired=F)
  wilcox.exact(speed ~ role, data=df.param[df.param$treat == "FM",], paired=F)
  wilcox.exact(speed ~ role, data=df.param[df.param$treat == "MM",], paired=F)
  kruskal.test(speed ~ treat, data = df.param[df.param$role == "Follower",])
  kruskal.test(speed ~ treat, data = df.param[df.param$role == "Leader",])
  
  
  ind.names <- unique(df$name) 
  mu = rho = c()
  for(i in 1:length(ind.names)){
    df.temp <- df.sep.focus[df.sep.focus$name == ind.names[i],]
    mle_wrpcauchy <- wrpcauchy.ml(na.omit(df.temp$angle), 0, 0, acc=1e-015)
    if(F){
      theta <- seq(-pi, pi, length=1000)
      probabl <- dwrpcauchy(theta, as.numeric(mle_wrpcauchy[1]), as.numeric(mle_wrpcauchy[2]))
      Range <- ceiling(max(probabl,na.rm=T))
      truehist(df.temp$angle, prob=T, breaks=seq(-pi, pi, length=100),
               ylim=c(0,Range), xlim=c(-pi,pi), 
               xlab="angle", ylab="density", axes=F, main=ind.names[i])
      points(theta, probabl, type="l")  
      axis(1, at=seq(-pi,pi,length=3), label=c("-pi", "0", "pi"))
      axis(2, las=1)
    }
    mu <- c(mu, as.numeric(mle_wrpcauchy[1]))
    rho <- c(rho, as.numeric(mle_wrpcauchy[2]))
  }
  df.param$mu = mu
  df.param$rho = rho
  
  tapply(df.param$mu, df.param[,c("treat", "role")], mean)
  tapply(df.param$rho, df.param[,c("treat", "role")], mean)
  as.vector(tapply(df.param$rho, df.param[,c("treat", "role")], mean))
  
  df.temp <- data.frame(
    treat = rep(c("FF", "FM", "MM"), 2),
    role = rep( c("F", "L"), each=3),
    speed = round(as.vector(tapply(df.param$speed, df.param[,c("treat", "role")], mean))*5, 3),
    sinuosity =round( as.vector(tapply(df.param$rho, df.param[,c("treat", "role")], mean)), 3)
  )
  
  write.csv((df.temp), file= file.path(PROJHOME, "/data/", "Parameters.csv"), row.names = F, col.names = F)
  
}

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
  df <- data.frame(
    Treat = c("FF-real", "FF-simulate", "FM-real", "MM-real", "MM-simulate"),
    prop60 = d$X300/100000
  )
  df$Treat <- factor(df$Treat, level = rev(df$Treat[c(3,1,2,4,5)]))
  
  ggplot(df, aes(y=Treat, x=prop60)) +
    geom_bar(stat="identity") +
    theme_bw() +
    theme(legend.position = "none", aspect.ratio = 1)+ 
    xlim(c(0,0.5))+
    xlab("Reencounter rate after 60 seconds")
  ggsave(filename = file.path(PROJHOME, "/plot/", paste0(today, "-SimulationRes.pdf")),
         width=4, height = 4, family="PT Sans")
  
}

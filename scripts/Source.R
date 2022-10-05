## R. speratus homosexual tandem analysis
## N. Mizumoto

#------------------------------------------------------------------------------#
# This file is for loading packages and Functions for other scripts.
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Packages
#------------------------------------------------------------------------------#
{
  library(data.table)
  library(stringr)
  
  library(viridis)
  library(ggplot2)
  require(grid)
  require(gridExtra)
  library(ggridges)
  
  library(exactRankTests)
  library(CircStats)
  
  library(survival)
  library("survminer")
  require(coxme)
  
  library(lme4)
  library(car)
  library(multcomp)
  library(exactRankTests)
  
  require(phytools)
  library(stringr)
  library(scales)
  
  library(extrafont)
  #font_import(pattern="PT")
  loadfonts()
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
{
  dir.create("./plot/")
  dir.create("./plot/trajectories/")
  dir.create("./plot/scale-trajectories/")
  dir.create("./plot/tandem-reverse/")
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Constants
#------------------------------------------------------------------------------#
{
  # general
  today <- Sys.Date()
  fps <- 5
  arena.size <- 85 # mm
  rda.place <- file.path("data/rda")
  
  # tandem analysis
  interaction.threshold <- 7
  min.tandem.sec <- 2
  threshold.moving.dis.tandem <- 30
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Functions
#------------------------------------------------------------------------------#
{
  # Standard Error
  se  <-  function(x){
    y  <-  x[!is.na(x)]  #  remove  the  missing  values
    sqrt(var(as.vector(y))/length(y))
  }
  
  angle_cal <- function(X, Y, Length){
    Ax <- (X[3:Length-1] - X[3:Length-2])
    Bx <- (X[3:Length] - X[3:Length-1])
    Ay <- (Y[3:Length-1] - Y[3:Length-2])
    By <- (Y[3:Length] - Y[3:Length-1])
    hugo <- (Ax * By - Ay * Bx + 0.000001)/abs(Ax * By - Ay * Bx + 0.000001)
    cos <- round((Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5,14)
    return(acos(cos)*hugo)
  }
  
  # Smoothing tandem vector 
  # so that tandem shorter than min.sec is treated as non tandem
  # vec = tandem, xvec = non tandem, output = tandem
  tandem.smoothing <- function(vec, xvec, min.sec){
    if(sum(vec)>0){
      timing <- which(vec)[c(T, diff(which(vec))>1)]
      end <- which(vec)[c(diff(which(vec))>1,T)]
      for(fi in 1:length(timing)){
        if(length( vec[timing[fi]:end[fi]]) < min.sec ){
          vec[timing[fi]:end[fi]] <- F
        }
      }
    }
    return(vec)
  }
  
  # Smoothing tandem vector 
  # so that non-interaction shorter than min.sec is treated as tandem
  # vec = non interaction, xvec = separation, output = !tandem
  tandem.smoothing2 <- function(vec, xvec, min.sec){
    if(sum(vec)>0){
      timing <- which(vec)[c(T, diff(which(vec))>1)]
      end <- which(vec)[c(diff(which(vec))>1,T)]
      for(fi in 1:length(timing)){
        if(any(xvec[timing[fi]:end[fi]])){
          vec[timing[fi]:end[fi]] <- T
        } else if(length( vec[timing[fi]:end[fi]]) < min.sec ){
          vec[timing[fi]:end[fi]] <- F
        }
      }
    }
    return(vec)
  }
  
  # Calcurate effect size 
  # based on the equation (10) in Nakagawa and Cuthill 2007 Biol Rev
  effectsize_d <- function(z, df, n1, n2){
    return((z*(n1+n2)) / (sqrt(n1*n2)*sqrt(df)))
  }
  
}
#------------------------------------------------------------------------------#


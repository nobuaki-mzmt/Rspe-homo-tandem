## R. speratus homosexual tandem analysis
## N. Mizumoto

#------------------------------------------------------------------------------#
# This file is for output all results, including plots and statistical analysis
#------------------------------------------------------------------------------#
rm(list = ls())

{
  source("scripts/Source.R")
  tandem.output()
  duration.output()
  movement.parameter()
  sim.plot()
  timedist.comp()
  group.tandem.plot()
  nestmate.analysis()
}

#------------------------------------------------------------------------------#
# This function reads df_tandem.rda and df_tandem_sum.rda, then
# 1) compare duration of tandems across pair combinations
# 2) examine leader-follower interactions: correlation between acceleration and inter-individual distance.
# 3) plot distribution of metrics for competition (rotation index, distance)
#------------------------------------------------------------------------------#
tandem.output <- function(){
  load("data/rda/df_tandem_sum.rda")
  load("data/rda/df_tandem.rda")
  
  ## 1. Comparison of tandem duration
  {
    
    df.tandem.sum$Treat <- factor(df.tandem.sum$Treat, levels = c("FM", "FF", "MM"))
    ggplot( df.tandem.sum, aes(x=Treat, y= TandemComp.Proportion, fill=Treat)) +
      #geom_boxplot()
      #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, binwidth=.01)+
      geom_point(alpha = 0.4, position = position_jitter(.05), size = 2) +
      geom_boxplot(outlier.shape = NA, alpha=0.3, width =.2) +
      guides(fill = F)+
      scale_fill_viridis(discrete=T)+
      ylim(c(0,1.01))+
      theme_bw()+
      xlab("") + ylab("Proportion of time spent in tandem")
    ggsave(filename = file.path("plot/", paste0(today, "-TandemProportion.pdf")),
           width=4, height = 4, family="PT Sans")
    
    e = 0.01
    y <- df.tandem.sum$TandemComp.Proportion
    df.tandem.sum$logit.TandemComp.Proportion <- (log((y+e)/(1-y+e)))
    plot(y, (log((y+e)/(1-y+e))))
    hist(log((y+e)/(1-y+e)))
    r <- lmer(logit.TandemComp.Proportion ~ Treat + (1|Colony), data=df.tandem.sum)
    res <- Anova(r)
    
    print("1. Comparison of tandem duration")
    print("lmer(logit.TandemComp.Proportion ~ Treat + (1|Colony), data=df.tandem.sum)")
    print(res)
  }
  
  ## 2. Correalation between acceleration and inter-individual distance
  {
    df.temp <- data.frame(
      df.tandem[, c(1:6)],
      role = rep(c("leader", 'follower'), each=dim(df.tandem)[1]),
      scheme = df.tandem$scheme,
      competition = df.tandem$competition,
      acc = c(df.tandem$acc.leader, df.tandem$acc.follower),
      ind.dis = df.tandem$ind.dis)
    
    set.seed(125)
    df.temp <- df.temp[sample(1:dim(df.temp)[1], 10000, replace = F),]
    df.temp$treat <- factor(df.temp$treat, levels = c("FM", "FF", "MM"))
    ggplot(df.temp[df.temp$scheme == 't' & !df.temp$competition,], aes(x = ind.dis, y=acc, col=role))+
      geom_point(alpha=0.1, size=0.75)+
      facet_wrap(vars(treat)) + 
      ylim(c(-15,15)) +
      xlim(c(0,11.5))+
      stat_smooth(se =T, method = "lm")+
      theme_bw() + theme(aspect.ratio = 1)+
      scale_color_viridis(discrete = T, end=0.5)+
      xlab("Distance between partners (mm)") +
      ylab("Acceleration (mm/sec2)")
    ggsave(filename = file.path("plot/", paste0(today, "-TandemAcceleration.pdf")),
           width=8, height = 4, family="PT Sans")
    
    df.temp <- data.frame(
      df.tandem[, c(1:6)],
      role = rep(c("leader", 'follower'), each=dim(df.tandem)[1]),
      scheme = df.tandem$scheme,
      competition = df.tandem$competition,
      acc = c(df.tandem$acc.leader, df.tandem$acc.follower),
      ind.dis = df.tandem$ind.dis)
    df.temp <- df.temp[df.temp$scheme == 't' & !df.temp$competition,]
    
    print("2. Correalation between acceleration and inter-individual distance")
    
    print("FM follower")
    r <- lmer(acc ~ ind.dis + (1|colony),
              data=df.temp[df.temp$treat=="FM" & df.temp$role=="follower",])
    print(r)
    print(Anova(r))
    
    print("FM leader")
    r <- lmer(acc ~ ind.dis + (1|colony),
              data=df.temp[df.temp$treat=="FM" & df.temp$role=="leader",])
    print(r)
    print(Anova(r))
    
    print("FF follower")
    r <- lmer(acc ~ ind.dis + (1|colony),
              data=df.temp[df.temp$treat=="FF" & df.temp$role=="follower",])
    print(r)
    print(Anova(r))
    
    print("FF leader")
    r <- lmer(acc ~ ind.dis + (1|colony),
              data=df.temp[df.temp$treat=="FF" & df.temp$role=="leader",])
    print(r)
    print(Anova(r))
    
    print("MM follower")
    r <- lmer(acc ~ ind.dis + (1|colony),
              data=df.temp[df.temp$treat=="MM" & df.temp$role=="follower",])
    print(r)
    print(Anova(r))
    
    print("MM leader")
    r <- lmer(acc ~ ind.dis + (1|colony),
              data=df.temp[df.temp$treat=="MM" & df.temp$role=="leader",])
    print(r)
    print(Anova(r))
  }
  
  ## 3. Distribution of metrics for competition
  {
    print("pairs showing competition for more than 5 seconds")
    print( table(subset(df.tandem.sum, Competition>5)$Treat) )
  }
  
  ## 4. obtain parameter d for simulation
  {
    df.tandem2 <- df.tandem[,c(1:6,15,16)]
    sep.timing <- which(df.tandem2$scheme=="s")[c(T, diff(which(df.tandem2$scheme=="s"))>1)]
    sep.end <- which(df.tandem2$scheme=="s")[c(diff(which(df.tandem2$scheme=="s"))>1,T)]
    res <- NULL
    for (j in 1:length(sep.timing)) {
      dftemp <- df.tandem2[ sep.timing[j]:sep.end[j], ]
      res <- c(res, mean(dftemp$ind.dis))
    }
    print("parameter d for simulation is:")
    print(mean(res))
  }
  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function reads df_tandem_duration.rda and df_separation_duration.rda, then
# 1) compare duration of each tandem and separation event across pair combinations
#------------------------------------------------------------------------------#
duration.output <- function(){
  odir = "plot/"
  load("data/rda/df_tandem_duration.rda")
  load("data/rda/df_separation_duration.rda")
  
  ## 1. Comparison of tandem duration
  {
    fit1 <- survfit(Surv(Duration, Cens) ~ Treat, 
                    type = "kaplan-meier", 
                    data = df.tandem.duration)
    ggsurvplot(
      fit = fit1,
      conf.int = T,
      xlab = "Duration (sec)", 
      ylab = "Tandem Prob",
      xlim = c(0,300),
      palette = viridis(3,option = "A"),
      legend = c(0.8,0.8),
      ggtheme = theme_bw(),
      data=df.tandem.duration
    )
    
    ggsave(paste0(odir, "Tandem-duration.pdf"), width=4, height = 3)
    
    f.name = "Duration.txt"
    sink(paste0(odir, f.name))
    cat("Tandem duration----------------------------------------\n")
    cat("coxme(Surv(Tandem.time, Cens) ~ Treat + (1|Source/Video))\n")
    m <- coxme(Surv(Duration, Cens) ~ Treat + (1|Colony/Video), data = df.tandem.duration)
    print(summary(m))
    print(Anova(m))
    
    cat("\n-----------------------------------------------\n")
  }
  
  ## 2. Correalation between acceleration and inter-individual distance
  {
    fit1 <- survfit(Surv(Duration, Cens) ~ Treat, 
                    type = "kaplan-meier", 
                    data = df.separate.duration)
    ## Rs
    ggsurvplot(
      fit = fit1,
      conf.int = T,
      xlab = "Duration (sec)", 
      ylab = "Separation Prob",
      xlim = c(0,150),
      palette = viridis(3,option = "B"),
      legend = c(0.8,0.8),
      ggtheme = theme_bw(),
      data=df.separate.duration
    )
    
    ggsave(paste0(odir, "Separate-duration.pdf"), width=4, height = 3)
    
    cat("Separation duration ----------------------------------------\n")
    cat("coxme(Surv(Separation.time, Cens) ~ Treat + (1|Source/Video))\n")
    m <- coxme(Surv(Duration, Cens) ~ Treat + (1|Colony/Video), data = df.separate.duration)
    print(summary(m))
    print(Anova(m))
    
    cat("\n-----------------------------------------------\n")
    sink()
  }
  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function reads df.sec.sum, df_sec_sum_sep.rda, df_sec_sum_tan.rda, df_sep_turn.rda then
# 1) plot time development of speed
# 2) perform comparison of speed
# 3) perform comparison of turning pattern
# 4) output paramters for simulations
#------------------------------------------------------------------------------#
movement.parameter <- function(){
  ## 1. time development of speed
  {
    load(paste0(rda.place, "/df_sec_sum.rda"))
    ## plot time development of speed
    col1 <- c(viridis(2, end = 0.5), viridis(2, begin=0.2, end = 0.6))
    col1 <- col1[c(2,4,3,1)]
    ggplot(df.sec.sum, aes(x=time, y=speed.mean, color=role, fill=role)) +
      geom_ribbon(aes(ymin=speed.mean.m.se, ymax=speed.mean.p.se), alpha=0.6) +
      facet_grid(treat~scheme) +
      scale_color_manual(values=col1) +
      scale_fill_manual(values=col1) +
      #scale_color_viridis(discrete = T, end=0.5) +
      #scale_fill_viridis(discrete = T, end=0.5) +
      theme_bw() + theme(aspect.ratio = 0.75) +
      scale_y_continuous(limits = c(0,17.5)) +
      xlab("Time (sec)") + ylab("Speed (mm/sec)")
    ggsave(filename = file.path("plot/", paste0(today, "-SpeedTimeDevelopment.pdf")),
           width=6, height = 8, family="PT Sans")
  }
  
  ## 2. Comparison of mean movement speed during the last 1min for tandem runs and first 1min after separation
  {
    ## After separation
    {
      load(paste0(rda.place, "/df_sec_sum_sep.rda"))
      ggplot(df.sep.sec.sum, aes(x=treat, y = speed.mean, fill=role, col=role)) + 
        geom_point(alpha = 0.4, position = position_dodge(width = 0.75) , size = 2) +
        geom_boxplot(alpha = 0.2) +
        scale_fill_viridis(discrete = T, end=0.5)+
        scale_color_viridis(discrete=T, end=0.5)+
        theme_bw() + theme(aspect.ratio = 1) +
        ylab("Speed mean (mm/sec)") + xlab("")
      ggsave(filename = file.path("plot/", paste0(today, "-SpeedAfterSeparation.pdf")),
             width=4, height = 4, family="PT Sans")
      
      print("After separation")
      print("Comparison between roles in FF")
      df_temp = df.sep.sec.sum[df.sep.sec.sum$treat == "FF",]
      r <- lmer(speed.mean ~ role + (1|colony), data=df_temp)
      print(Anova(r))
      res = summary(r)
      res_coef = res$coefficients
      res_z = res_coef["roleLeader","t value"]
      effectsize_d(z = res_z, df = dim(df_temp)[1]-2, n1 = table(df_temp$role)[1], n2 = table(df_temp$role)[1])
      
      print("Comparison between roles in FM")
      df_temp = df.sep.sec.sum[df.sep.sec.sum$treat == "FM",]
      r <- lmer(speed.mean ~ role + (1|colony), data=df_temp)
      print(Anova(r))
      res = summary(r)
      res_coef = res$coefficients
      res_z = res_coef["roleLeader","t value"]
      effectsize_d(z = res_z, df = dim(df_temp)[1]-2, n1 = table(df_temp$role)[1], n2 = table(df_temp$role)[1])

      print("Comparison between roles in MM")
      df_temp = df.sep.sec.sum[df.sep.sec.sum$treat == "MM",]
      r <- lmer(speed.mean ~ role + (1|colony), data=df_temp)
      print(Anova(r))
      res = summary(r)
      res_coef = res$coefficients
      res_z = res_coef["roleLeader","t value"]
      effectsize_d(z = res_z, df = dim(df_temp)[1]-2, n1 = table(df_temp$role)[1], n2 = table(df_temp$role)[1])

      
      
      print("Comparison among pair combination in follower")
      df_temp = df.sep.sec.sum[df.sep.sec.sum$role == "Follower" & df.sep.sec.sum$treat != "FM",]
      r <- lmer(speed.mean ~ treat + (1|colony), data=df_temp)
      print(Anova(r))
      res = summary(r)
      res_coef = res$coefficients
      res_z = res_coef["treatMM","t value"]
      effectsize_d(z = res_z, df = dim(df_temp)[1]-2, n1 = table(df_temp$role)[1], n2 = table(df_temp$role)[1])
      
      print("Comparison among pair combination in leader")
      df_temp = df.sep.sec.sum[df.sep.sec.sum$role == "Leader" & df.sep.sec.sum$treat != "FM",]
      r <- lmer(speed.mean ~ treat + (1|colony), data=df_temp)
      print(Anova(r))
      res = summary(r)
      res_coef = res$coefficients
      res_z = res_coef["treatMM","t value"]
      effectsize_d(z = res_z, df = dim(df_temp)[1]-2, n1 = table(df_temp$role)[1], n2 = table(df_temp$role)[1])
      
    }
      
    ## Tandem run
    {
      load(paste0(rda.place, "/df_sec_sum_tan.rda"))
      ggplot(df.tan.sec.sum, aes(x=treat, y = speed.mean, fill=role, col=role)) + 
        geom_point(alpha = 0.4, position = position_dodge(width = 0.75) , size = 2) +
        geom_boxplot(alpha = 0.2) +
        scale_fill_viridis(discrete = T, end=0.5)+
        scale_color_viridis(discrete=T, end=0.5)+
        theme_bw() + theme(aspect.ratio = 1) +
        ylab("Speed mean (mm/sec)") + xlab("")
      ggsave(filename = file.path("plot/", paste0(today, "-SpeedDuringTandem.pdf")),
             width=4, height = 4, family="PT Sans")
      
      print("Duaring tandem")
      print("Comparison between roles in FF")
      r <- lmer(speed.mean ~ role + (1|colony), data=df.tan.sec.sum[df.tan.sec.sum$treat == "FF",])
      print(Anova(r))
      
      print("Comparison between roles in FM")
      r <- lmer(speed.mean ~ role + (1|colony), data=df.tan.sec.sum[df.tan.sec.sum$treat == "FM",])
      print(Anova(r))
      
      print("Comparison between roles in MM")
      r <- lmer(speed.mean ~ role + (1|colony), data=df.tan.sec.sum[df.tan.sec.sum$treat == "MM",])
      print(Anova(r))
      
      print("Comparison among pair combination in follower")
      r <- lmer(speed.mean ~ treat + (1|colony), data=df.tan.sec.sum[df.tan.sec.sum$role == "Follower",])
      print(Anova(r))
      multicomparison<-glht(r,linfct=mcp(treat="Tukey"))
      print(summary(multicomparison))
      
      print("Comparison among pair combination in leader")
      r <- lmer(speed.mean ~ treat + (1|colony), data=df.tan.sec.sum[df.tan.sec.sum$role == "Leader",])
      print(Anova(r))
      multicomparison<-glht(r,linfct=mcp(treat="Tukey"))
      print(summary(multicomparison))
    }
  }
  
  ## 3. Comparison of movement sinuousity (turning pattern)
  {
    load(paste0(rda.place, "/df_sep_turn.rda"))
    ggplot(df.sep.turn, aes(x=treat, y = rho, fill=role, col=role)) + 
      geom_point(alpha = 0.4, position = position_dodge(width = 0.75) , size = 2) +
      geom_boxplot(alpha = 0.2) +
      scale_fill_viridis(discrete = T, end=0.5)+
      scale_color_viridis(discrete=T, end=0.5)+
      theme_bw() + theme(aspect.ratio = 1) +
      ylab("Speed mean (mm/sec)") + xlab("")
    
    print("After separation")
    print("Comparison between roles in FF")
    df_temp = df.sep.turn[df.sep.turn$treat == "FF",]
    r <- lmer(rho ~ role + (1|colony), data=df_temp)
    print(Anova(r))
    res = summary(r)
    res_coef = res$coefficients
    res_z = res_coef["roleLeader","t value"]
    effectsize_d(z = res_z, df = dim(df_temp)[1]-2, n1 = table(df_temp$role)[1], n2 = table(df_temp$role)[1])
    
    print("Comparison between roles in FM")
    df_temp = df.sep.turn[df.sep.turn$treat == "FM",]
    r <- lmer(rho ~ role + (1|colony), data=df_temp)
    print(Anova(r))
    res = summary(r)
    res_coef = res$coefficients
    res_z = res_coef["roleLeader","t value"]
    effectsize_d(z = res_z, df = dim(df_temp)[1]-2, n1 = table(df_temp$role)[1], n2 = table(df_temp$role)[1])
    
    
    print("Comparison between roles in MM")
    df_temp = df.sep.turn[df.sep.turn$treat == "MM",]
    r <- lmer(rho ~ role + (1|colony), data=df_temp)
    print(Anova(r))
    res = summary(r)
    res_coef = res$coefficients
    res_z = res_coef["roleLeader","t value"]
    effectsize_d(z = res_z, df = dim(df_temp)[1]-2, n1 = table(df_temp$role)[1], n2 = table(df_temp$role)[1])

  }
  
  ## 4. Obtain parameters for simulations
  {
    df.temp <- data.frame(
      treat = rep(c("FF", "FM", "MM"), 2),
      role = rep( c("F", "L"), each=3),
      speed = round(as.vector(tapply(df.sep.sec.sum$speed.mean, df.sep.sec.sum[,c("treat", "role")], mean)), 3),
      sinuosity =round( as.vector(tapply(df.sep.turn$rho, df.sep.turn[,c("treat", "role")], mean)), 3)
    )
    write.csv((df.temp), file= file.path("data/sim/", "Parameters.csv"), row.names = F, col.names = F)
  }
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function reads simulation results (SimRes_rep100000_sec60.csv), and then
# 1) plot the results
#------------------------------------------------------------------------------#
sim.plot <- function(){
  d <- data.frame(fread("data/sim/SimRes_sec60.csv", header=T))
  
  ## time development of encounter rates
  df.sim <- NULL
  
  for(i in 1: (dim(d)[1]) ){
    df <- data.frame(
      comb     = (i-1) %% 5,
      ind1     = d[i,1], 
      ind2     = d[i,2],
      sim_iter = d[i,3],
      time     = 1:300,
      rate     = as.numeric(d[i,4:303])
    )
    df.sim = rbind(df.sim, df)
  }

  df.sim.sum = data.frame(
    df.sim[df.sim$sim_iter==0, c(1:3,5)],
    rate_mean = as.vector(tapply(df.sim$rate/100000, df.sim[,c("time", "comb")], mean)),
    rate_sd = as.vector(tapply(df.sim$rate/100000, df.sim[,c("time", "comb")], sd)),
    count = as.vector(tapply(df.sim$rate, df.sim[,c("time", "comb")], sum))
  )
  
  df.sim.sum$comb[df.sim.sum$comb==0] = "FF_observed"
  df.sim.sum$comb[df.sim.sum$comb==1] = "FF_simulated"
  df.sim.sum$comb[df.sim.sum$comb==2] = "FM_observed"
  df.sim.sum$comb[df.sim.sum$comb==3] = "MM_observed"
  df.sim.sum$comb[df.sim.sum$comb==4] = "MM_simulated"
  
  ggplot(df.sim.sum, 
         aes(x=time/5, y=rate_mean,
             group=comb, fill=as.factor(comb))) +
    geom_ribbon(aes(ymin=rate_mean-rate_sd*3, ymax=rate_mean+rate_sd*3), alpha = 0.2) +
    geom_line(aes(col=as.factor(comb)), size = 0.5) +
    scale_color_viridis(discrete = T) + 
    scale_fill_viridis(discrete = T) +
    ylim(c(0,0.5)) +
    theme_bw() +
    theme(aspect.ratio = 0.75) +
    xlab("Time (sec)") + ylab("Probability to encounter")
  
  ggsave(filename = file.path("plot/", paste0("SimulationEncounterRatio.pdf")),
         width=4, height = 4, family="PT Sans")

  ## histogram for the time required for encounter
  {
    ### XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  }
  
  ## sensitivity analysis
  {
    SA.path.list = list.files("data/sim/SA", full.names = T)
    SA.file.list = list.files("data/sim/SA", full.names = F)
    for(i_SA in 1:length(SA.file.list)){
      d <- data.frame(fread(SA.path.list[i_SA], header=T))
      
      df.sim <- NULL
      
      for(i in 1: (dim(d)[1]) ){
        df <- data.frame(
          comb     = (i-1) %% 5,
          ind1     = d[i,1], 
          ind2     = d[i,2],
          sim_iter = d[i,3],
          time     = 1:300,
          rate     = as.numeric(d[i,4:303])
        )
        df.sim = rbind(df.sim, df)
      }
      
      df.sim.sum = data.frame(
        df.sim[df.sim$sim_iter==0, c(1:3,5)],
        rate_mean = as.vector(tapply(df.sim$rate/100000, df.sim[,c("time", "comb")], mean)),
        rate_sd = as.vector(tapply(df.sim$rate/100000, df.sim[,c("time", "comb")], sd))
      )
      
      df.sim.sum$comb[df.sim.sum$comb==0] = "FF_observed"
      df.sim.sum$comb[df.sim.sum$comb==1] = "FF_simulated"
      df.sim.sum$comb[df.sim.sum$comb==2] = "FM_observed"
      df.sim.sum$comb[df.sim.sum$comb==3] = "MM_observed"
      df.sim.sum$comb[df.sim.sum$comb==4] = "MM_simulated"
      
      ggplot(df.sim.sum, 
             aes(x=time/5, y=rate_mean,
                 group=comb, fill=as.factor(comb))) +
        geom_ribbon(aes(ymin=rate_mean-rate_sd*3, ymax=rate_mean+rate_sd*3), alpha = 0.2) +
        geom_line(aes(col=as.factor(comb)), size = 0.5) +
        scale_color_viridis(discrete = T) + 
        scale_fill_viridis(discrete = T) +
        ylim(c(0,0.5)) +
        theme_bw() +
        theme(aspect.ratio = 0.75) +
        xlab("Time (sec)") + ylab("Probability to encounter") +
        ggtitle(str_remove(SA.file.list[i_SA], ".csv"))
      
      ggsave(filename = file.path("plot/", str_replace(SA.file.list[i_SA], "csv", "pdf")),
             width=4, height = 4, family="PT Sans")
      
    }
  }

}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Analysis for comparison of the time required for encounter (Fig. S6)
#------------------------------------------------------------------------------#
timedist.comp <- function(){
  
  # data format
  load("data/rda/df_separation_duration.rda")
  
  d <- data.frame(fread("data/sim/SimRes_sec60.csv", header=T))
  df.sim <- NULL
  for(i in 1: (dim(d)[1]) ){
    df <- data.frame(
      comb     = (i-1) %% 5,
      ind1     = d[i,1], 
      ind2     = d[i,2],
      sim_iter = d[i,3],
      time     = 1:300,
      rate     = as.numeric(d[i,4:303])
    )
    df.sim = rbind(df.sim, df)
  }
  
  df.sim.sum = data.frame(
    df.sim[df.sim$sim_iter==0, c(1:3,5)],
    count = as.vector(tapply(df.sim$rate, df.sim[,c("time", "comb")], sum))
  )
  df.sim.sum$comb[df.sim.sum$comb==0] = "FF_observed"
  df.sim.sum$comb[df.sim.sum$comb==1] = "FF_simulated"
  df.sim.sum$comb[df.sim.sum$comb==2] = "FM_observed"
  df.sim.sum$comb[df.sim.sum$comb==3] = "MM_observed"
  df.sim.sum$comb[df.sim.sum$comb==4] = "MM_simulated"
  df.sim.sum$count = c(0,diff(df.sim.sum$count))
  df.sim.sum[df.sim.sum$time == 1,]$count = 0
  
  # plot
  sim_label = c("FF_observed", "FF_simulated", "FM_observed", "MM_observed", "MM_simulated")
  sim_cols = c("#d01b1bff", "#964242ff", "black", "#47abd8ff", "#9596ecff")
  for(i in 1:5){
    ggplot(data = df.sim.sum[df.sim.sum$comb == sim_label[i],]) +
      geom_histogram(aes(x=time/5, y=..count..,weight=count),
                     binwidth = 0.2, alpha = 0.4, fill = sim_cols[i])+
      xlim(c(0,60)) +
      ylab("Count") +
      xlab("Time until encounter (sec)") +
      theme_bw() + 
      theme(legend.position = "none", aspect.ratio = 0.75) +
      ggtitle(paste("Sim", sim_label[i]))
    ggsave(paste0("plot/", sim_label[i], "-time_dstr.pdf"),
           width=4, height = 4, family="PT Sans")
  }
  
  emp_label = c("FM", "FF", "MM")
  emp_cols = viridis(3, option = "B")
  for(i in 1:3){
    ggplot(data = df.separate.duration[df.separate.duration$Treat==emp_label[i],]) +
      geom_histogram(aes(x=Duration),
                     binwidth = 0.2, alpha = 0.4, fill = emp_cols[i])+
      xlim(c(0,60)) +
      ylab("Count") +
      xlab("Time until encounter (sec)") +
      theme_bw() + 
      theme(legend.position = "none", aspect.ratio = 0.75) +
      ggtitle(paste("Emp", emp_label[i]))
    ggsave(paste0("plot/", emp_label[i], "-time_dstr.pdf"),
           width=4, height = 4, family="PT Sans")
  }
  
  # ks test
  df.sim$comb[df.sim$comb==0] = "FF_observed"
  df.sim$comb[df.sim$comb==1] = "FF_simulated"
  df.sim$comb[df.sim$comb==2] = "FM_observed"
  df.sim$comb[df.sim$comb==3] = "MM_observed"
  df.sim$comb[df.sim$comb==4] = "MM_simulated"
  
  df_ks = NULL
  for(i_label in 1:5){
    for(i_iter in 0:9){
      df.temp = df.sim[df.sim$comb == sim_label[i_label] & df.sim$sim_iter == i_iter,]
      sim.vec = rep( df.temp$time, c(0,diff(df.temp$rate)))/5
      if(i_label < 3){
        emp.vec = df.separate.duration[df.separate.duration$Treat=="FF",]$Duration
      } else if (i_label == 3){
        emp.vec = df.separate.duration[df.separate.duration$Treat=="FM",]$Duration
      } else {
        emp.vec = df.separate.duration[df.separate.duration$Treat=="FF",]$Duration
      }
      ks_res <- ks.test(sim.vec, emp.vec)
      df_ks.temp = data.frame(label = sim_label[i_label],
                              iter = i_iter,
                              ks_d = ks_res$statistic
      )
      df_ks <- rbind(df_ks, df_ks.temp)
    }
  }
  df_ks$label <- factor(df_ks$label, levels = sim_label[c(4,5,1,2,3)])
  ggplot(df_ks, aes(x=label, y=ks_d)) +
    geom_jitter(width = 0.25, cex=0.5)+
    ylab("KS D statistics") +
    xlab("") +
    theme_bw() + 
    theme(legend.position = "none", aspect.ratio = 2.5) 
  ggsave(paste0("plot/KS_D.pdf"),
         width=4, height = 5, family="PT Sans")
  
  print(wilcox.exact(ks_d ~ label, paired=F, 
               data=df_ks[df_ks$label == "FF_observed" | df_ks$label == "FF_simulated",]))
  print(wilcox.exact(ks_d ~ label, paired=F, 
               data=df_ks[df_ks$label == "MM_observed" | df_ks$label == "MM_simulated",]))
  
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# This function reads Grouptandem.rda (the results of group observation), and then
# 1) compare num pf tandem engaging individuals among 5F5M, 10F and 10M
#------------------------------------------------------------------------------#
group.tandem.plot <- function(){
  print("Group tandem")
  
  load("data/rda/Grouptandem.rda")
  r <- glmer(cbind(IndEngageTandem, NonEng) ~ treat + (1|colony), family="binomial", data=d)
  print(Anova(r))
  multicomparison<-glht(r,linfct=mcp(treat="Tukey"))
  print(summary(multicomparison))
  
  effectsize_d(z = 2.912, df = 16, n1 = 9, n2 = 9)
  effectsize_d(z = 3.740, df = 16, n1 = 9, n2 = 9)
  effectsize_d(z = 0.935, df = 16, n1 = 9, n2 = 9)
  
  ggplot(data=d, aes(x=treat, y=IndEngageTandem/10)) +
    geom_dotplot(binaxis = "y", stackdir='center', dotsize=0.1, binwidth = 0.1)+
    ylim(c(0,1))+
    stat_summary(fun.data=mean_se, 
                 geom="pointrange", color="red") +
    theme_bw() + ylab("Prop of individuals engaging in tandems") +
    theme(legend.position = "none", aspect.ratio = 1) 
  ggsave(filename = file.path("plot/", paste0(today, "-GroupTandemProportion.pdf")),
         width=4, height = 4, family="PT Sans")
}  
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Analysis for nestmate experiments (Text S1)
#------------------------------------------------------------------------------#
nestmate.analysis <- function(){
  load("data/rda/df_nestmate.rda")
  
  print("heterosexual tandem")
  print("comp num of runs observed")
  r <- glm(count ~ type, data = d[d$type == "h_nest" | d$type == "h_nonnest", ], family="poisson")
  print(Anova(r))
  summary(r)
  z = 1.085
  dis = (z*(6+6)) / (sqrt(6*6)*sqrt(11))
  print("effect size ")
  print(z/sqrt(z^2+11))
  
  print("male-male tandem")
  print("comp num of runs observed")
  r <- glm(count ~ type, data = d[d$type == "m_nest" | d$type == "m_nonnest", ], family="poisson")
  print(Anova(r))
  summary(r)
  z = 3.375
  dis = (z*(6+6)) / (sqrt(6*6)*sqrt(11))
  print("effect size ")
  print(z/sqrt(z^2+11))
  
  print("num of observation is")
  print( tapply(d$count, d$type, sum) )
  
}
#------------------------------------------------------------------------------#






















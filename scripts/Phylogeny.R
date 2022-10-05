## Tandem run phylogenetic comparative analysis
## N. Mizumoto

#------------------------------------------------------------------------------#
# This file is for phylogenetic comparative analysis
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Setup
#------------------------------------------------------------------------------#
{
  rm(list = ls())
  source("scripts/Source.R")
  
  ## Data
  # d = read.table("clipboard",header =T)
  # save(d, file = "data/rda/TermiteTandemInformation.rda")
  load("data/rda/TermiteTandemInformation.rda")
  
  ## Phylogeny
  tree_no3rd <- read.nexus("data/phylogeny/Termitephylogeny.tre")
  
  # drop outgroup and fossil species
  label <- tree_no3rd$tip.label
  label[label=="Coptotermes_priscus"] <-  "Coptotermes_priscus_fossil"
  tree_no3rd$tip.label <- label
  length(label)
  
  out.group <- is.na(str_locate(label, "Ter")[,1] | str_locate(label, "ter")[,1] | str_locate(label, "Cryptocercus")[,1])
  termite.tree <- drop.tip(tree_no3rd, label[out.group])
  termite.tree <- ladderize(termite.tree)
  
  label<- termite.tree$tip.label
  fossil <- !is.na(str_locate(label, "fossil")[,1])
  termite.tree <- drop.tip(termite.tree, label[fossil])
  label <- termite.tree$tip.label
  
  label <- str_replace(label, "_.*", "")
  label <- str_replace(label, "Homalotermes", "Homallotermes")
  termite.tree$tip.label <- label
  
  plot(termite.tree, cex=0.5)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Plot the results
#------------------------------------------------------------------------------#
plot_result <- function(data, termite.tree, filename){
  
  # Matching phylogeny and tandem data
  label <- termite.tree$tip.label
  df <- data[!is.na(data$Leader) | !is.na(data$Tandem),]
  Res <- NULL
  remove.label <- NULL
  for(i in 1:length(label)){
    Genus <- label[i]
    if(sum(df$Genus == Genus) == 0){
      #print(Genus)
      remove.label <- c(remove.label, Genus)
    } else {
      Res <- rbind(Res, df[df$Genus == Genus,])
    }
  }
  row.names(Res) <- 1:dim(Res)[1]

  termite.tree <- drop.tip(termite.tree, remove.label)
  termite.tree <- ladderize(termite.tree)
  
  ## get tiles
  is_tip <- termite.tree$edge[,2] <= length(termite.tree$tip.label)
  ordered_tips <- termite.tree$edge[is_tip, 2]
  ordered_tips.label <- termite.tree$tip.label[ordered_tips]
  
  tandem <- Res$Tandem
  leader <- Res$Leader
  leader.male <- (Res$Leader=="male" | Res$Leader=="both")*1
  leader.female <- (Res$Leader=="female" | Res$Leader=="both")*1
  
  names(tandem) = names(leader) = names(leader.male) = names(leader.female) = Res$Genus
  ordered_tamdem <- as.numeric(rev(tandem[ordered_tips.label]))
  ordered_leader <- rev(leader[ordered_tips.label])
  ordered_leader.male <- rev(leader.male[ordered_tips.label])
  ordered_leader.female <- rev(leader.female[ordered_tips.label])
  
  cols = viridis(2, end=0.5, alpha=0.5)
  col.show_tandem <- alpha(c(cols[1], "white"), 0.5)
  col.show_female <- alpha(c("pink", "white"), 0.5)
  col.show_male <- alpha(c("light blue", "white"), 0.5)
  data.tile <- c(matrix(c(col.show_tandem[2-ordered_tamdem],
                          col.show_female[2-ordered_leader.female],
                          col.show_male[2-ordered_leader.male]), 3, byrow=T))
  pdf(file.path("plot/", paste0(filename, "-PhylogenyPanel.pdf")), width=4, height=6, family = "PT Serif", paper = "a4")
  show_col(data.tile, ncol = 3, labels =F)
  dev.off()
  
  col_tandem = c("white", "#44015480")
  col_female = c("white", "pink")
  col_male   = c("white", "light blue")
  
  pdf(file.path("plot/", paste0(filename, "-Phylogeny.pdf")), width=4, height=6, family = "PT Serif")
  
  plotTree(termite.tree,fsize=0.7,ftype="i",lwd=1)
  node.size = 0.5
  
  fitER_tandem <- ace(tandem, termite.tree, model="ER", type="discrete")
  nodelabels(node=1:termite.tree$Nnode+Ntip(termite.tree),
             pie=fitER_tandem$lik.anc,piecol=col_tandem,cex=node.size, adj = -4)
  
  fitER_female <- NULL
  if( sum(leader.female == 0, na.rm=T) > 0 ){
    fitER_female <- ace(leader.female, termite.tree, model="ER", type="discrete")
    nodelabels(node=1:termite.tree$Nnode+Ntip(termite.tree),
               pie=fitER_female$lik.anc,piecol=col_female,cex=node.size)
  }

  fitER_male <- ace(leader.male, termite.tree, model="ER", type="discrete")
  nodelabels(node=1:termite.tree$Nnode+Ntip(termite.tree),
             pie=fitER_male$lik.anc,piecol=col_male,cex=node.size, adj = 5)
  dev.off()
  
  plotTree(termite.tree,fsize=0.7,ftype="i",lwd=1)
  nodelabels(text=1:termite.tree$Nnode,node=1:termite.tree$Nnode+Ntip(termite.tree))
  
  
  return(list(fitER_tandem, fitER_female, fitER_male))
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# main calcuration
#------------------------------------------------------------------------------#
# Cryptocercus lacks tandem run
dftemp <- d
dftemp[dftemp$Genus == "Cryptocercus", c(2,3)] = c(0, NA)
dftemp[dftemp$Genus == "Mastotermes", c(2,3)]  = c(1, "both")
dftemp[dftemp$Genus == "Hodotermes", c(2,3)]   = c(1, "male")
fitres = plot_result(data = dftemp, termite.tree = termite.tree, filename = "Crylack")
fitres[[1]]$lik.anc
fitres[[2]]$lik.anc
fitres[[3]]$lik.anc

# Cryptocercus has tandem run
dftemp <- d
dftemp[dftemp$Genus == "Cryptocercus", c(2,3)] = c(1, "both")
dftemp[dftemp$Genus == "Mastotermes", c(2,3)]  = c(1, "both")
dftemp[dftemp$Genus == "Hodotermes", c(2,3)]   = c(1, "male")
fitres = plot_result(data = dftemp, termite.tree = termite.tree, filename = "CryTan")
fitres[[1]]$lik.anc
fitres[[2]]$lik.anc
fitres[[3]]$lik.anc


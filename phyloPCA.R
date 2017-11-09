# Script for Paniw et al. XXXXXX - Appendix S2 - PCA analyses

# This script is divided into two parts:
# PART A (lines 38-441): Calculate life-history traits and relate the demographic data (needed for PCA) to phylogenies obtained from open sources
# PART B (lines 442-644): Perform a phylogenetically-informed PCA on the demographic data

# Author: Maria Paniw & Roberto Salguero-Gomez
# Created: 15 May 2016

########## 
#Clean memory
rm(list=ls(all=TRUE))

# load neccessary packages 

library(plyr)
library(ggplot2)
library(mgcv)
library(stringr)
library(fpc)
library("Cairo")
library(car)
library(ggrepel)
library(ape)
library(caper)
library(phytools)
library(agricolae)
library(popbio)
library(MASS)

# Set working directory 

setwd("/Users/mariapaniw/Dropbox/TempAutoProject/SuppMat")

#####################################################
## PART A - LINK SPECIES TO PHYLOGENY
########################################################

# Read in the data containing vital rates

load("matsMean")
################# Life history traits: necessary functions (written by R. Salguero-Gomez)

# where the arguments of the function above are taken from use 
# lifeTimeReprodEvents and makeLifeTable in https://github.com/jonesor/compadreDB/tree/master/Functions
# lx: age-specific survivorship
# x: ages
# LaF: mean age at sexual maturity
# LaC: mean age at clonal maturity - this may not be useful to you, as I think you discarded clonal species?
# QSD1: quasi-stable stage distribution to cut of mortality plateaus


lifeTimeRepEvents <- function(matU, matF, startLife ){
  #Function to determine probability of reaching reproduction, age at maturity and reproductive lifespan (Code adapted from H. Caswell's matlab code):
  
  uDim = dim(matU)[1]
  surv = colSums(matU)
  repLifeStages = colSums(matF)
  repLifeStages[which(repLifeStages>0)] = 1
  
  if(missing(matF) | missing(matU)){stop('matU or matF missing')}
  if(sum(matF,na.rm=T)==0){stop('matF contains only 0 values')}
  
  #Probability of survival to first reprod event
  Uprime = matU
  Uprime[,which(repLifeStages==1)] = 0
  Mprime = matrix(0,2,uDim)
  for (p in 1:uDim[1]) {
    if (repLifeStages[p]==1) Mprime[2,p] = 1 else
      Mprime[1,p] = 1-surv[p]
  }
  Bprime = Mprime%*%(ginv(diag(uDim)-Uprime))
  pRep = Bprime[2,startLife]
  
  out = data.frame(pRep = pRep)
  
  #Age at first reproduction (La; Caswell 2001, p 124)
  D = diag(c(Bprime[2,]))
  Uprimecond = D%*%Uprime%*%ginv(D)
  expTimeReprod = colSums(ginv(diag(uDim)-Uprimecond))
  La = expTimeReprod[startLife]
  
  out$La = La
 
  return(out)
}

# get lx and mx
makeLifeTable<-function(matU, matF = NULL, matC = NULL, startLife, nSteps = 1000){
  
  matDim = ncol(matU)
  
  #Age-specific survivorship (lx) (See top function on page 120 in Caswell 2001):
  matUtemp = matU
  survivorship = array(NA, dim = c(nSteps, matDim))
  for (o in 1:nSteps){
    survivorship[o, ] = colSums(matUtemp %*% matU)
    matUtemp = matUtemp %*% matU
  }
  
  lx = survivorship[, startLife]
  lx = c(1, lx[1:(length(lx) - 1)])
  
  #Make room for dx and qx under assumption of 0.5 in age gap distributions
  
  #Start to assemble output object
  out = data.frame(x = 0:(length(lx)-1),lx = lx)
  
  if(!missing(matF)){
    if(sum(matF,na.rm=T)==0){
      warning("matF contains only 0 values")
    }
    #Age-specific fertility (mx, Caswell 2001, p. 120)
    ageFertility = array(0, dim = c(nSteps, matDim))
    fertMatrix = array(0, dim = c(nSteps, matDim))
    matUtemp2 = matU
    e = matrix(rep(1, matDim))
    for (q in 1:nSteps) {
      fertMatrix = matF %*% matUtemp2 * (as.numeric((ginv(diag(t(e) %*% matUtemp2)))))
      ageFertility[q, ] = colSums(fertMatrix)
      matUtemp2 = matUtemp2 %*% matU
    }  
    mx = ageFertility[, startLife]
    mx = c(0, mx[1:(length(mx) - 1)])
    out$mx = mx
  }
  
  if(!missing(matC)){
    if(sum(matC,na.rm=T)==0){
      warning("matC contains only 0 values")
    }
    #Age-specific clonality (cx)
    ageClonality = array(0, dim = c(nSteps, matDim))
    clonMatrix = array(0, dim = c(nSteps, matDim))
    matUtemp2 = matU
    e = matrix(rep(1, matDim))
    for (q in 1:nSteps) {
      clonMatrix = matC %*% matUtemp2 * (as.numeric((ginv(diag(t(e) %*% matUtemp2)))))
      ageClonality[q, ] = colSums(clonMatrix)
      matUtemp2 = matUtemp2 %*% matU
    }  
    cx = ageClonality[, startLife]
    cx = c(0, cx[1:(length(cx) - 1)])
    out$cx = cx
  }
  
  return(out)
}

qsdConverge <- function(matU, conv = 0.05, startLife = 1, nSteps = 1000){
  #Function to determine the cutoff age at quasi-convergence for lx and mx (Code adapted from H. Caswell's matlab code):
  
  uDim = dim(matU)
  eig = eigen.analysis(matU)
  qsd = eig$stable.stage
  qsd = as.numeric(t(matrix(qsd / sum(qsd))))
  
  #Set up a cohort
  nzero = rep(0, uDim[1]) #Set a population vector of zeros
  nzero[startLife] = 1 #Set the first stage to = 1
  n = nzero #Rename for convenience
  
  #Iterate the cohort (n= cohort population vector, p = proportional structure)
  dist = p = NULL
  survMatrix1 <- matU
  for (j in 1:nSteps){ #j represent years of iteration
    p = n / sum(n) #Get the proportional distribution
    dist[j] = 0.5 * (sum(abs(p - qsd)))
    n = survMatrix1 %*% n #Multiply the u and n matrices to iterate
  }
  #Find the ages for convergence to conv. (default = 0.05).
  #i.e. within 5% of the QSD.
  convage = min(which(dist < conv))
  return(convage) 
}


demetriusEntropy <- function(A,lx,mx,QSD1){
  r <- log(max(Re(eigen(A)$values)))
  lxmx<-lx[1:QSD1]*mx[1:QSD1]
  for (i in 1:QSD1) {
    if (lxmx[i]==0) {lxmx[i]<-1}
  }
  loglxmx<-log(lxmx)
  loglxmx[which(lxmx==0)]<-NA
  #demetrius<-exp(r*c(1:length(lxmx)))%*%(lxmx)
  demetrius<-abs(sum(lxmx*loglxmx)/sum(lxmx))
  return(demetrius)
}

###################### Calculate life-history traits for MPMs

LH_traits=array(NA,c(length(mats),5))
sp.name=rep(NA,length(mats)) # SpeciesAuthor
spA.name=rep(NA,length(mats))# SpeciesAccepted

### calculate sensitivity and elasticity 

for(x in 1:length(mats)){
  
  sp.name[x]=mats[[x]]$species
  spA.name[x]=mats[[x]]$speciesAccepted

  F.mat=mats[[x]]$matF
  U.mat=mats[[x]]$matU
  surv=colSums(U.mat)
  U.mat.g=U.mat
  for(xx in 1:ncol(U.mat.g)){
    U.mat.g[,xx]=U.mat.g[,xx]/surv[xx]
  }
 
  U.mat.g[!is.finite(U.mat.g)]=0
  
  if(colSums(U.mat)[ncol(U.mat)]==1){U.mat[,ncol(U.mat)]<-U.mat[,ncol(U.mat)]*0.995}
  
  mean.mat=U.mat+F.mat
  
  life = which(mats[[x]]$class%in%c("PR","R"))[1] # start life at 1st pre or reproductive stage
 
  # round mean mat for some species where generation time doesn't work otherwise:
  
  if(is.na(generation.time(mean.mat,F.mat))){
    
    mean.mat=round(mean.mat,4)
    # generation time
    
    LH_traits[x,1] = generation.time(mean.mat,F.mat)
  }else{
    
    LH_traits[x,1] = generation.time(mean.mat,F.mat)
  }
  
  La=tryCatch({
    
    lifeTimeRepEvents(U.mat,F.mat,startLife=life)$La
  },error=function(e) NA)
  
  # Age at sexual maturity
  LH_traits[x,2] <- La
  
  
  # Degree of iteroparity
  a=makeLifeTable(U.mat,matF=F.mat,startLife=life)
  QSD1=qsdConverge(U.mat, conv = 0.05, startLife = life)
  
  ifelse(is.infinite(QSD1),QSD1<-1,QSD1)
  
  a$mx[is.na(a$mx)]=0
  S= demetriusEntropy(mean.mat,a$lx,a$mx,QSD1)
  
  LH_traits[x,3]<- S
  
  # net reproductive rate
  LH_traits[x,4] = net.reproductive.rate(mean.mat,F.mat)
  
  
  # vital rates weighted by stage distribution 
  
  stable.stage=eigen.analysis(mean.mat)$stable.stage/sum(eigen.analysis(mean.mat)$stable.stage)

  stable.stage3=stable.stage[which(colSums(F.mat)>0)]/sum(stable.stage[which(colSums(F.mat)>0)])
  #Mean reproduction
  LH_traits[x,5]=sum(colSums(F.mat)[which(colSums(F.mat)>0)]*stable.stage3)
  
}

colnames(LH_traits)=c("T","La","S","Ro","Phi")
rownames(LH_traits)=sp.name

LH_traits[,3]=LH_traits[,3]+1

for(a in 1:ncol(LH_traits)){
  
  
  # Arcsine transform for PCA analysis
  LH_traits[,a]=log(LH_traits[,a])
  
}


El=as.data.frame(LH_traits)

El$species=spA.name
El$ID=rownames(El)

# PREPARE PHYLOGENY

El$species <- gsub("_"," ", El$species)
El$species <- gsub("[0-9]+", "", El$species)
El$species <- gsub("\\svar.+", "", El$species)#remove varieties bc many times phylogenies are only resolved at the species level


# read in phylogenetic tree 
algaeTree <- read.tree("phylogeny/AlgaePhylo.tre")
plantTree <- read.tree("phylogeny/PlantPhylo.tre")
animalTree <- read.nexus("phylogeny/AnimalPhylo.nexus")

#In case you need to drop species that are no present in the data but are present in the tree. E.g.:
plantTree <- drop.tip(plantTree, 
                      c("Gracilaria_gracilis", "Mazzaella_splendens", "Gelidium_sesquipedale"))

supertree <- bind.tree(algaeTree, plantTree, where="root")
supertree <- bind.tree(supertree, animalTree, where="root")
supertree <- compute.brlen(supertree)

plot(supertree)

# Stiching it all together, plants, algae, and animals
supertree$tip.label <- gsub("_", " ", supertree$tip.label)
supertree$tip.label <- gsub("[0-9]+", "", supertree$tip.label)
supertree$tip.label <- gsub(" ott", "", supertree$tip.label)
tree <- supertree
supertree$tip.label[which(duplicated(supertree$tip.label))] # should be 0

# What species are in the data but not tree? Due to mismatched names!

missingSpecies=El$species[!El$species%in%tree$tip.label]# 40 species missing 

# Change species names manually:

El$species[El$species=="Cottus sp."]="Cottus bairdii"
El$species[El$species=="Thalassarche melanophris"]="Thalassarche melanophrys"
El$species[El$species=="Odocoileus virginianus subsp. borealis"]="Odocoileus virginianus"
El$species[El$species=="Xenosaurus sp."]="Xenosaurus grandis"
El$species[El$species=="Fucus vesiculosus"]="Fucus serratus"
El$species[El$species=="Brassica napus"]="Brassica insularis"
El$species[El$species=="Asplenium adulterinum"]="Asplenium scolopendrium"
El$species[El$species=="Anthyllis vulneraria alpicola"]="Anthyllis vulneraria"
El$species[El$species=="Antirrhinum molle lopesianum"]="Antirrhinum subbaeticum"
El$species[El$species=="Centaurea stoebe"]="Centaurea corymbosa"
El$species[El$species=="Silene glaucifolia pseudoviscosa"]="Silene acaulis"
El$species[El$species=="Polemonium van-bruntiae"]="Polemonium vanbruntiae"
El$species[El$species=="Silene douglasii oraria"]="Silene regia"
El$species[El$species=="Geonoma pohliana weddelliana"]="Geonoma schottiana"
El$species[El$species=="Helianthemum juliae"]="Helianthemum motae"
El$species[El$species=="Sambucus sieboldiana"]="Sambucus racemosa subsp. sieboldiana"
El$species[El$species=="Vella pseudocytisus paui"]="Vella pseudocytisus"
El$species[El$species=="Escobaria robbinsorum"]="Escobaria robbinsorum"
El$species[El$species=="Mammillaria napia"]="Mammillaria napina"
El$species[El$species=="Choerospodnias axillaris"]="Choerospondias axillaris"
El$species[El$species=="Rangifer tarandus subsp. platyrhynchus"]="Rangifer tarandus"
El$species[El$species=="Cryptantha flava"]="Oreocarya flava"
El$species[El$species=="Limonium geronense"]="Limonium narbonense"
El$species[El$species=="Limonium malacitanum"]="Limonium lilacinum"
El$species[El$species=="Mimulus cardinalis"]="Erythranthe cardinalis"
El$species[El$species=="Mimulus lewisii"]="Erythranthe lewisii"
El$species[El$species=="Solidago altissima"]="Solidago canadensis"
El$species[El$species=="Viola persicifolia"]="Phoma persicifolia"
El$species[El$species=="Persoonia glaucescens"]="Perigonia glaucescens"
El$species[El$species=="Rosmarinus tomentosus"]="Coprinopsis cinerea"
El$species[El$species=="Nothofagus fusca"]="Fuscospora fusca"   
El$species[El$species=="Sapium sebiferum"]="Triadica sebifera"   
El$species[El$species=="Phalacrocorax auritus"]="Nannopterum auritus"  
El$species[El$species=="Thalassarche melanophris"]="Thalassarche melanophrys"    
El$species[El$species=="Cornu aspersa"]="Helix aspersa"
El$species[El$species=="Proclossiana eunomia"]="Boloria eunomia"                                
El$species[El$species=="Clethrionomys rufocanus"]="Myodes rufocanus"
El$species[El$species=="Podocnemis lewyana "]="Podocnemis lewyana"
El$species[El$species=="Calathea micans"]="Goeppertia micans"
El$species[El$species=="Escobaria robbinsorum"]="Escobaria robbinsiorum"
El$species[El$species=="Stenocactus crispatus"]="Ferocactus crispatus"
El$species[El$species=="Carya sinensis"]="Annamocarya sinensis"
El$species[El$species=="Falco peregrinus subsp. anatum"]="Falco peregrinus"
El$species[El$species=="Ovis canadensis subsp. sierrae"]="Ovis canadensis"
El$species[El$species=="Cottus bairdi"]="Cottus bairdii"

# remove duplicates (populations) initially - add them later as branches to the tree

El2 = El[-which(duplicated(El$species)),] 

duplicates=El[which(duplicated(El$species)),"species"] # remember duplicates to add them later as branches

# Give duplicates a unique name (same name as for branches)
El[which(duplicated(El$species)),"species"]=paste(El[which(duplicated(El$species)),"species"],"1")

# Drop all the tips in the tree that are not in the data
drop<-tree$tip.label[which(!tree$tip.label%in%El$species)]

#Should be TRUE
length(drop)+length(unique(El2$species))==length(tree$tip.label)

# Trim tree of not needed species, if any at all
smalltree<-drop.tip(tree, drop)


all(sort(El2$species)==sort(smalltree$tip.label))

#Final tree

tree=smalltree
tree$node.label<-as.character(1:length(tree$node.label))

#The tree must have branch lengths to be able to build the phylogenetic variance-covariance matrix into the function phyl.pca for the phylogenetically-informed PCA
tree <- compute.brlen(tree)

# To use the tree in for an MCMC analyses, it needs to be rooted
# By default the tree should be unrooted, but we can root it
#Using the comand root or root.tree
tree <- root(tree, tree$tip.label[1])

tree <- multi2di(tree, random=FALSE)
is.rooted(tree)

any(duplicated(tree$node.label))
tree$node.label <- unique(tree$node.label)

#Another possible issue could be that the distance between taxa
#Could be 0, and this is not enabled for some analyses (including MCMCglmm)
#Thus we can solve this by adding an 0.001 to all tips
tree$edge.length <-  tree$edge.length + 0.001

#Transform the tree into an ultrametric one
is.ultrametric(tree)
tree <- chronos(tree, lambda=0, model="correlated") 

### Add populations as extra branches
for(a in 1:length(duplicates)){
  
  name=paste(duplicates[a],"1", sep="_")
  
  tree<-bind.tip(tree, name, where=which(tree$tip.label==duplicates[a]))
  tree$tip.label=gsub("_"," ", tree$tip.label)
}

is.ultrametric(tree)
class(tree) <- "phylo"

any(duplicated(tree$node.label))
tree$node.label <- unique(tree$node.label)

# Plot the final tree

library("Cairo")

CairoPDF("plots/phylotree.pdf", width = 16, height =67,pointsize=16)          

plot(tree,cex=0.5)

dev.off()
###########################################################################
#######################################################
#### PART B PHYLOGENETICALLY - INFORMED PCA

# Phylo PCA
row.names(El)=El$species

pcaData=scale(El[,c(1:5)])

pca=phyl.pca(tree,pcaData,method="lambda",mode="corr") 

summary(pca)

pca$lambda# 0.61 
pca$logL

diag(pca$Eval)^2 # Apply Kaiser's criterion to determine how many PCA axes to keep

ncomp=2 # keep the first 2 axes

# Varimax correction on the axes

rawLoadings     <- pca$L[,1:ncomp]

# find the variance maximizing roation of loadings 
rotatedLoadings <- varimax(rawLoadings)$loadings

name1=0.344
name2=0.333

# Create the inverse of the loading to calculate new scores: data multiplied by rotatio matrix 
invLoadings     <- t(pracma::pinv(rotatedLoadings))
scores          <- pcaData %*% invLoadings


x <- list() 
x$scores <- scores[,1:2]
colnames(x$scores)=c("PC1", "PC2")
x$scores[,1:2] <- (-1)*x$scores[,1:2]
x$loadings <- rotatedLoadings[,1:2]
colnames(x$loadings)=c("PC1", "PC2")
x$loadings[,1:2] <- (-1)*x$loadings[,1:2]

#### SAVE THE PCA SCORES FOR REGRESSION ANALYSES (APPENDIX sigma_r)

write.csv(cbind(x$scores,El$ID),"PCAscores.csv",row.names=F)

### Get standard deviation of Pagel's lambda by simple bootstrap

nsim=100

Plambda=rep(NA,nsim)

for(i in 1:nsim){
  
  pcaData1=pcaData[sample(1:nrow(pcaData),nrow(pcaData),replace=T),]
  pcaData1=pcaData1[!duplicated(pcaData1), ]
  Plambda[i]=phyl.pca(tree,pcaData1,method="lambda",mode="corr")$lambda 
  
}

##############################################################

# Plot results

### With point labels (population ID)
data1 <- data.frame(x$scores,point.lab=1:nrow(x$scores))

PCbiplot <- function(PC, x="PC1", y="PC2") {

  plot <- ggplot(data1, aes_string(x=x, y=y)) 
  plot <- plot + geom_point(aes(size=3,color="grey20"),alpha=0.5)
  plot <- plot + geom_hline(aes(yintercept=0), size=.2) + geom_vline(aes(xintercept=0), size=.2)
  plot <- plot + scale_size(range = c(0, 6)) 
  datapc <- data.frame(varnames=colnames(El[1:5]), PC$loadings)
  mult <- min(
    (max(data1[,y]) - min(data1[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data1[,x]) - min(data1[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .35 * mult* (get(x)),
                      v2 = .35 * mult* (get(y))
  )
  plot<-plot + geom_text_repel(aes(label=point.lab))
  plot <- plot + geom_segment(data=datapc[1:2,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=2, color="sienna1")
  plot <- plot + geom_segment(data=datapc[3:5,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=2, color="palegreen3")
  plot <- plot+theme_bw()+theme(panel.grid = element_blank())+ylab(paste("PCA 2 (",round(name2,2)*100,"%)",sep=""))+xlab(paste("PCA 1 (",round(name1,2)*100,"%)",sep=""))
  plot<- plot+theme(axis.text = element_text(size=18))+theme(axis.title = element_text(size=20))+theme(legend.title = element_text(size=16, face="bold"),
                                                                                                       legend.text = element_text(size=16),
                                                                                                       legend.key.size = unit(1, "lines"))
  
  plot
}

PCbiplot(x)

# No points (NOTE: color based on Sv1 here not shown for simplicity)

library(fields)

PCbiplot <- function(PC, x="PC1", y="PC2") {
  # PC being a prcomp object
  
  plot <- ggplot(data1, aes_string(x=x, y=y)) 
  plot <- plot + geom_point(aes(size=3,color="grey20"),alpha=0.5)
  plot <- plot + geom_hline(aes(yintercept=0), size=.2) + geom_vline(aes(xintercept=0), size=.2)
  plot <- plot + scale_size(range = c(0, 6)) 
  datapc <- data.frame(varnames=colnames(El[1:5]), PC$loadings)
  mult <- min(
    (max(data1[,y]) - min(data1[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data1[,x]) - min(data1[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .35 * mult* (get(x)),
                      v2 = .35 * mult* (get(y))
  )
  
  plot<-plot + geom_text_repel(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5,  color="black", point.padding = unit(1, 'lines'))
  plot <- plot + geom_segment(data=datapc[1:2,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=2, color="sienna1")
  plot <- plot + geom_segment(data=datapc[3:5,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=2, color="palegreen3")
  plot <- plot+theme_bw()+theme(panel.grid = element_blank())+ylab(paste("PCA 2 (",round(name2,2)*100,"%)",sep=""))+xlab(paste("PCA 1 (",round(name1,2)*100,"%)",sep=""))
  plot<- plot+theme(axis.text = element_text(size=18))+theme(axis.title = element_text(size=20))+theme(legend.position="none")
  
  plot
}

PCbiplot(x)

################# Check that PCA axes capture most variation among Sv1

library(caper)

data1$species=El$ID
# load Sv1 (obtained from code simLambdaTAC.R)

sens=read.csv("Sv1_PCA.csv")
sens$sens=abs(sens$sens)
sens=sens[sens$sens>0,]

# Sum Sv1 across vital rates
sens.mu=aggregate(sens~species,data=sens[sens$cv==0.5&sens$f==0.65,],sum)

sens.mu$species=as.character(sens.mu$species)

sens.mu=left_join(sens.mu,data1,by="species")

sens.mu$SP_tree=left_join(sens.mu,El,by = c("species" = "ID"))$species.y

sens.mu$sens=log(sens.mu$sens)

# Prepare data for pGLS
comp_data <- comparative.data(phy = tree,
                              data =sens.mu,
                              names.col = SP_tree,
                              vcv=TRUE)

### Look at different lambdas explaining phylogenic effect 

pgls_mod0 <- pgls(TAC~PC1*PC2, data = comp_data, lambda =c(0.001))
summary(pgls_mod0)

pgls_mod1 <- pgls(TAC~PC1*PC2, data = comp_data, lambda = c(0.15))
summary(pgls_mod1)

pgls_mod2 <- pgls(TAC~PC1*PC2, data = comp_data, lambda = c(0.2))
summary(pgls_mod2)

pgls_mod3 <- pgls(TAC~PC1*PC2, data = comp_data, lambda = c(0.3))
summary(pgls_mod3)

pgls_mod4 <- pgls(TAC~PC1*PC2, data = comp_data, lambda = c(0.4))
summary(pgls_mod4)

pgls_mod5 <- pgls(TAC~PC1*PC2, data = comp_data, lambda = c(0.99))
summary(pgls_mod5)

AIC(pgls_mod0,pgls_mod1,pgls_mod2,pgls_mod3,pgls_mod4,pgls_mod5) # from model 2 on, AIC is not sufficiently decreased

# Search for lambda using ML

pgls_mod1 <- pgls(sens~PC1*PC2, data = comp_data, lambda = "ML", bounds=list(lambda=c(0.001,0.9)))
summary(pgls_mod1)

pgls_mod2 <- pgls(sens~PC1+PC2 +I(PC1^2)*I(PC2^2), data = comp_data, lambda = "ML", bounds=list(lambda=c(0.01,0.99)))
summary(pgls_mod2)

AIC(pgls_mod1,pgls_mod2)

pgls_mod3 <- pgls(sens~PC1+PC2 +I(PC1^2)*I(PC2^2)+PC1*I(PC1^2), data = comp_data, lambda = "ML", bounds=list(lambda=c(0.01,0.99)))
summary(pgls_mod3)

AIC(pgls_mod1,pgls_mod2,pgls_mod3)# model 2 is best

mod_profile <- pgls.profile(pgls_mod2) # model performance doesn't improve signficantly after a Pagel's lambda of 0.15
plot(mod_profile)

#  plot predictions using different Pagel's lambda
pgls_mod7a=pgls(sens~PC1+PC2 +I(PC1^2)*I(PC2^2), data = comp_data, lambda = c(0.001))

pgls_mod7b=pgls(sens~PC1+PC2 +I(PC1^2)*I(PC2^2), data = comp_data, lambda = c(0.15))

pgls_mod7c=pgls(sens~PC1+PC2 +I(PC1^2)*I(PC2^2), data = comp_data, lambda = c(0.3))

### plot
library(directlabels)
library(fields)
library(lattice)
library(latticeExtra)
library(reshape2)

data=NULL

sub1=comp_data$data

fit = Tps(cbind(sub1$PC1,sub1$PC2), sub1$sens)
pred = predictSurface(fit)
new=melt(pred$z)

new$PC1=rep(pred$x,pred$ny)
new$PC2=rep(pred$y,each=pred$nx)
new$Pl="0"
new$sens=predict(pgls_mod7a,newdata = new)[,1]
new$sens[is.na(new$value)]=NA

data=rbind(data,new)

### 2
fit = Tps(cbind(sub1$PC1,sub1$PC2), sub1$sens)
pred = predictSurface(fit)
new=melt(pred$z)

new$PC1=rep(pred$x,pred$ny)
new$PC2=rep(pred$y,each=pred$nx)
new$Pl="0.15"
new$sens=predict(pgls_mod7b,newdata = new)[,1]
new$sens[is.na(new$value)]=NA

data=rbind(data,new)

### 3
fit = Tps(cbind(sub1$PC1,sub1$PC2), sub1$sens)
pred = predictSurface(fit)
new=melt(pred$z)

new$PC1=rep(pred$x,pred$ny)
new$PC2=rep(pred$y,each=pred$nx)
new$Pl="0.5"
new$sens=predict(pgls_mod7c,newdata = new)[,1]
new$sens[is.na(new$value)]=NA

data=rbind(data,new)

data$Pl=factor(data$Pl)
levels(data$Pl)=c(expression(paste("Pagel's ", lambda, "= 0")),expression(paste("Pagel's ", lambda, "= 0.15")),expression(paste("Pagel's ", lambda, "= 0.5")))


sub1$sens=exp(sub1$sens)


breaks=c(-1,-4,-6.5)
ggplot(data, aes(PC1, PC2, z =sens))+
  geom_raster(aes(fill = sens),interpolate=F) +
  geom_point(data=sub1,aes(size=sens),shape=1)+
  facet_wrap( ~ Pl,ncol=3,scales="fixed", labeller = label_parsed)+ 
  scale_fill_gradientn(limits=c(-6.5,-1),colours=tim.colors(128),na.value="white",breaks=breaks)+
  guides(color=F,size=F,fill=F)+
  xlab("PCA 1")+ylab("PCA 2")+
  scale_x_continuous( expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size=18,colour="black"))+
  theme(axis.title = element_text(size=24))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.9),
        plot.title = element_text(size=21,lineheight=.8, face="bold"))+
  geom_hline(aes(yintercept=0), size=.2) + 
  geom_vline(aes(xintercept=0), size=.2)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(panel.spacing = unit(0.7, "lines"))+
  theme(strip.background =element_blank(),
        strip.text = element_text(size=20))

# No significant differences in results when Pagel's lambda increases
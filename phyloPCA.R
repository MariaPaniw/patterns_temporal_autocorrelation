# Script for Paniw et al. XXXXXX - Appendix S2 - PCA analyses

# This script is divided into two parts:
# PART A (lines 38-305): Calculate deterministic vital-rate elasticities and relate the demographic data (needed for PCA) to phylogenies obtained from open sources
# PART B (lines 307-353): Perform a phylogenetically-informed PCA on the demographic data

# Author: Maria Paniw & Rob Salguero-Gomez
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


# Set working directory 

setwd("C:/Users/Maria/Dropbox/TempAutoProject")

#####################################################
## PART A - LINK SPECIES TO PHYLOGENY
########################################################

# Read in the data containing vital rates

raw=read.csv("vitalRates.csv")

# subset the data frame to have only the vital rates

vr=raw[,c("s1","s2","s3","s4","g21","g31","g32","g43","g34","f13","f23","f33")]

# Turn NAs into 0 (to facilitate further calculations) 

vr[is.na(vr)]=0

### Calculate elasticities (and sensitivities since the PCA can also be performed on sensitivities)

sens.elast=function(vitalRates){
  
 
  sens=array(NA,c(nrow(vitalRates),ncol(vitalRates)))
  
  elast=array(NA,c(nrow(vitalRates),ncol(vitalRates)))
  
 
  ### calculate sensitivity and elasticity 
  
  for(x in 1:nrow(vitalRates)){
    
    # Create matrix from the mean vital rates per population
    mean.mat=matrix(c(vitalRates[x,"s1"]*(1-vitalRates[x,"g21"]-vitalRates[x,"g31"]-0),vitalRates[x,"s1"]*vitalRates[x,"g21"],vitalRates[x,"s1"]*vitalRates[x,"g31"],0,
                      0,vitalRates[x,"s2"]*(1-0-vitalRates[x,"g32"]-0),vitalRates[x,"s2"]*vitalRates[x,"g32"],0,
                      vitalRates[x,"f13"],vitalRates[x,"f23"],vitalRates[x,"f33"]+vitalRates[x,"s3"]*(1-vitalRates[x,"g43"]),vitalRates[x,"s3"]*vitalRates[x,"g43"],
                      0,0,vitalRates[x,"s4"]*vitalRates[x,"g34"],vitalRates[x,"s4"]*(1-vitalRates[x,"g34"])),
                    4,4,byrow=F)
    
    # get sensitivities and elasticities
    # for the popbio vitalsens() functions, two parameters are required: a list with vital rate values and names
    # and a description (elements) of how the vital rates are put together in a transition matrix 
    
    vital.rates=lapply(seq_len(ncol(vr[x,])),function(i) vitalRates[x,i])
    
    names(vital.rates)=c("s1","s2","s3","s4","g21","g31","g32","g43","g34","f13","f23","f33")
    
    elements=expression(s1*(1-g21-g31),s1*g21,s1*g31,0,
                        0,s2*(1-g32),s2*g32,0,
                        f13,f23,f33+s3*(1-g43),s3*g43,
                        0,0,s4*g34,s4*(1-g34))
    
    sens[x,]<-vitalsens(elements, vital.rates)$sensitivity
    
    elast[x,]<-vitalsens(elements, vital.rates)$elasticity
    
    
  }
  
  colnames(sens)=colnames(elast)=rownames(vitalsens(elements, vital.rates))
  rownames(sens)=rownames(elast)=raw$SpeciesAuthor
  
  return(list(sens=sens,elast=elast))
}

sensElastList=sens.elast(vr)

# Exctract elasticities and calculate relative elasticities

elast=abs(sensElastList$elast)# we are interested in the magnitude, not direction of change, so we take absolute values

elast.sum=as.numeric(apply(elast,1,sum))# sum of elasticities across vital rates

rel.elast=apply(elast,2,function(x) x/elast.sum)

rel.elast[is.na(rel.elast)]=0

colnames(rel.elast)=c("s1","s2","s3","s4","g21","g31","g32","g43","g34","f13","f23","f33")


### use only vital rates that are present in > 10 % of populations

1-apply(rel.elast,c(2),function(x) length(x[x==0])/length(x)) # for what % populations are vital-rates non-zero

El=cbind(rel.elast[,-c(4,6,8,9,12)])

colnames(El)=c(
  as.character(interaction("REL",colnames(vr)[-c(4,6,8,9,12)],sep="x")))


for(a in 1:7){
  
  
  # Arcsine transform for PCA analysis
  El[,a]=asin((El[,a]))
  
}


El[is.na(El)]=0
El[is.infinite(El)]=0

El=as.data.frame(El)

El$species=as.character(raw$SpeciesAccepted)
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

### Order El

El=El[order(match(El$ID, raw$SpeciesAuthor)),]
all(raw$SpeciesAuthor==El$ID)

# Phylo PCA
row.names(El)=El$species

pcaData=scale(El[,c(1:7)])

pca=phyl.pca(tree,pcaData,method="lambda",mode="corr") 

summary(pca)

pca$lambda# 0.185 
pca$logL

name1=(diag(pca$Eval)/sum(pca$Eval))[1]
name2=(diag(pca$Eval)/sum(pca$Eval))[2]

diag(pca$Eval)^2 # Apply Kaiser's criterion to determine how many PCA axes to keep

ncomp=2 # keep the first 2 axes

# Varimax correction on the axes

rawLoadings     <- pca$L[,1:ncomp]

# find the variance maximizing roation of loadings 
rotatedLoadings <- varimax(rawLoadings)$loadings

# Create the inverse of the loading to calculate new scores: data multiplied by rotatio matrix 
invLoadings     <- t(pracma::pinv(rotatedLoadings))
scores          <- pcaData %*% invLoadings

x <- list() 
x$scores <- scores[,1:2]
colnames(x$scores)=c("PC1", "PC2")
x$loadings <- rotatedLoadings[,1:2]
colnames(x$loadings)=c("PC1", "PC2")

#### SAVE THE PCA SCORES FOR REGRESSION ANALYSES (APPENDIX S3)

write.csv(x$scores,"PCAscores.csv",row.names=F)

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
  datapc <- data.frame(varnames=colnames(El[2:8]), PC$loadings)
  mult <- min(
    (max(data1[,y]) - min(data1[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data1[,x]) - min(data1[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .35 * mult* (get(x)),
                      v2 = .35 * mult* (get(y))
  )
  plot<-plot + geom_text_repel(aes(label=point.lab))
  plot <- plot + geom_segment(data=datapc[1:3,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=2, color="sienna1")
  plot <- plot + geom_segment(data=datapc[4:5,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=2, color="skyblue2")
  plot <- plot + geom_segment(data=datapc[6:7,], aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), size=2, color="palegreen3")
  
  plot <- plot+theme_bw()+theme(panel.grid = element_blank())+ylab(paste("PCA 2 (",round(name2,2)*100,"%)",sep=""))+xlab(paste("PCA 1 (",round(name1,2)*100,"%)",sep=""))
  plot<- plot+theme(axis.text = element_text(size=18))+theme(axis.title = element_text(size=20))+theme(legend.title = element_text(size=16, face="bold"),
                                                                                                       legend.text = element_text(size=16),
                                                                                                       legend.key.size = unit(1, "lines"))
  
  plot
}
library("Cairo")

CairoPDF("plots/PCA_points_phylo.pdf", width = 24, height =21,pointsize=16)          


PCbiplot(x)

dev.off()


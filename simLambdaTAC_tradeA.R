# Script for Paniw et al. XXXXXX - Appendix S4 - Simulations of temporal autocorrelation in vital rates assuming tradeoff A

# This script simulates temporal autocorrelation assuming a tradeoff structure between survival and reproduction
# and calculate sensitivity of the stochastic growht rate to changes in correlation structure

# Author: Maria Paniw

########## 
# Clean memory
rm(list=ls(all=TRUE))

### load libraries

library(ggplot2)
library(lme4)
library(reshape2)
library(MASS)
library(Cairo)
library(plyr)
library(pbkrtest)

# set working directory 

setwd("C:/Users/Maria/Dropbox/TempAutoProject")

#####################################################################################
######################################################################
############# PART A - SIMULATIONS

raw=read.csv("vitalRates.csv")

# subset the data frame to have only the vital rates

vr=raw[,c("s1","s2","s3","s4","g21","g31","g32","g43","g34","f13","f23","f33")]

# Turn NAs into 0

vr[is.na(vr)]=0

# another data frame with vital rates that are correlated 
vr.cor=vr[,c("s3","s4","g43","g34","f13","f23","f33")]


# Create some parameters necessary for the simulation (loop)

# Autocorrelation:

v1=c(-0.3,0.3)

# Coefficient of variation

cv=c(0.2,0.5,0.8)

### Define how long simulations run and what the characteristics of the environments are

# Define simulation time (to run faster simulations, the user may decrease both tr and ts)
tr=10000 #the discard time
ts=90000 # how many years to keep

# run each simulation for ts+tr years 
trun=ts+tr

# The different environments 
env=c("good","bad")
env.n=1:2

# long-term frequency of good condition = first value of right eigenvector describing 
# stationary distribution of environmental states (see Tuljapurkar & Haridas, Ecol Lett, 2006 for more detail)

f=c(0.35,0.65)

#### Define trade-offs

# the covariances are given by the trade-off structures

a.cor=c(0.2,0.6,0.8)

#### Functions for truncated distributions (to new averages of vital rates in good and bad environments)

qtrunc <- function(p, spec, a = -Inf, b = Inf, ...)
{
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
  return(tt)
}

rtrunc <- function(n, spec, a = -Inf, b = Inf, ...)
{
  x <- u <- runif(n, min = 0, max = 1)
  x <- qtrunc(u, spec, a = a, b = b,...)
  return(x)
}

### Function to get positive definite matrix 
pos.def = function(origMat){ 
  
  cholStatus <- try(u <- chol(origMat), silent = T) 
  
  cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
  
  # fix the correl matrix
  
  newMat <- origMat
  
  iter <- 0
  
  while (cholError) {
    
    iter <- iter + 1
    
    cat("iteration ", iter, "\n")
    
    # replace -ve eigen values with small +ve number
    
    newEig <- eigen(newMat)
    
    newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)
    
    # create modified matrix eqn 5 from Brissette et al 2007, inv = transp for
    
    # eig vectors
    
    newMat <- newEig$vectors %*% diag(newEig2) %*% t(newEig$vectors)
    
    # normalize modified matrix eqn 6 from Brissette et al 2007
    
    newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
    
    # try chol again
    
    cholStatus <- try(u <- chol(newMat), silent = TRUE)
    
    cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
    
  }
  
  newMat}


# The output array will have 6 dimensions: correlation x species x vital rates x v1 x cv x f

sp.sub=1:nrow(raw) # the user may wish to work with a smaller subset of species to speed up simulations

lambda.s=array(NA,c(length(a.cor),length(sp.sub),ncol(vr.cor),length(v1),length(cv),length(f)))

### ACTUAL SIMULATIONS

# First main loop: vital-rate correlation

for(sc in 1:length(a.cor)){
  # Define vital-rate correaltion matrix 
  sigma <- matrix(c(1,0,0,0,0,0,0,0,0,0,0,0,
                    0,1,0,0,0,0,0,0,0,0,0,0,
                    0,0,1,0,0,0,0,0,0,-a.cor[sc],-a.cor[sc],-a.cor[sc],
                    0,0,0,1,0,0,0,0,0,-a.cor[sc],-a.cor[sc],-a.cor[sc],
                    0,0,0,0,1,0,0,0,0,0,0,0,
                    0,0,0,0,0,1,0,0,0,0,0,0,
                    0,0,0,0,0,0,1,0,0,0,0,0,
                    0,0,0,0,0,0,0,1,0,a.cor[sc],a.cor[sc],a.cor[sc],
                    0,0,0,0,0,0,0,0,1,-a.cor[sc],-a.cor[sc],-a.cor[sc],
                    0,0,-a.cor[sc],-a.cor[sc],0,0,0,a.cor[sc],-a.cor[sc],1,0,0,
                    0,0,-a.cor[sc],-a.cor[sc],0,0,0,a.cor[sc],-a.cor[sc],0,1,0,
                    0,0,-a.cor[sc],-a.cor[sc],0,0,0,a.cor[sc],-a.cor[sc],0,0,1),nrow=12,ncol=12,byrow=T)
  
  colnames(sigma)=rownames(sigma)=c("s1","s2","s3","s4","g21","g31","g32","g43","g34","f13","f23","f33")
  
  sigma=sigma[-c(1,2,5,6,7),-c(1,2,5,6,7)]# keep only the relevant non-zero correlations
  
  # Second main loop: frequency of good 
  for(z in 1:length(f)){
    
    # Third main loop: coef of var 
    for(i in 1:length(cv)){
      
      # Fourth main loop: autocorrelation
      for(j in 1:length(v1)){
        
        # Fifth main loop: species
        
        for(x in 1:length(sp.sub)){
          
          ### Focus on the vital rates that have non-zero entries 
          if(any(vr.cor[sp.sub[x],]>0)){
            sigma.sub=pos.def(sigma[which(vr.cor[sp.sub[x],]>0),which(vr.cor[sp.sub[x],]>0)])
            colnames(sigma.sub)=rownames(sigma.sub)=colnames(sigma[which(vr.cor[sp.sub[x],]>0),which(vr.cor[sp.sub[x],]>0)])
            if(!any(sigma.sub!=0&sigma.sub!=1)) sigma.sub<-NULL
            
          }else{sigma.sub=NULL}
          
          
          # Sixth main loop: vital rates
          
          for(y in 1:ncol(vr.cor)){
            
            #############################################
            ###############################
            # perturb the vr data frame to create bad and good condition mean vital rates
            
            bad.mean=vr.cor[sp.sub[x],]
            
            good.mean=vr.cor[sp.sub[x],]
            
            vr.mean=vr[sp.sub[x],]
            
            dat= vr.cor[sp.sub[x],]
            
            good.cor=bad.cor=NULL
            
            ########## BETWEEN-STATE VARIATION 
            
            ### START WITH FECUNDITIES
            if(y>=5){ # For fecundity, keep regular CV to define variation
              
              if(bad.mean[,y]>0&dat[,"s3"]>0&!is.null(sigma.sub)){ #good/bad environment matrices only for non-zero fecundities
                
                # simulate from a gamma distribution (FENCUDITY)
                mu=dat[,y]
                var=(cv[i]*dat[,y])^2
                
                # convert mean and variance to gamma parameters
                shape= mu^2/var
                scale=var/mu
                
                # new means for good and bad are the averages of two truncated gamma distributions
                # both with a mean given by the vital rate at its average across environments 
                # and variance given by the CV 
                
                bad.mean[,y]=mean(rtrunc(10000,spec="gamma",a=0,b=mu,shape=shape,scale=scale)) # upper bound truncated at the mean 
                good.mean[,y]=mean(rtrunc(10000,spec="gamma",a=mu,b=Inf,shape=shape,scale=scale)) # lower bound truncated at the mean
                
                # Define correlation structure 
                cor=sigma.sub[rownames(sigma.sub)[rownames(sigma.sub)!=colnames(vr.cor)[y]],colnames(vr.cor)[y]]
                
                names(cor)=rownames(sigma.sub)[rownames(sigma.sub)!=colnames(vr.cor)[y]]
                
                cor=cor[cor>(abs(min(sigma.sub))-0.1)|cor<(min(sigma.sub)+0.1)]
                
                good.cor=bad.cor=matrix(as.numeric(dat[,1:4]),1,4)
                
                colnames(good.cor)=colnames(bad.cor)=colnames(vr.cor)[1:4]
                
                
                for(cr in 1:length(cor)){
                  
                  # simulate from a beta distribution 
                  # (for the vital rates that are correlated with fecundity)
                  # by making their variation depend on the correlation coefficient cor
                  
                  mu2=dat[,names(cor)[cr]]
                  maxCV2=sqrt(mu2*(1-mu2))/mu2
                  var2=abs(cor[cr])*(cv[i]*maxCV2*mu2)^2
                  
                  # convert mean and variance to beta parameters
                  a2= mu2*((mu2*(1-mu2))/var2-1)
                  b2=(1-mu2)*((mu2*(1-mu2))/var2-1)
                  
                  if(cor[cr]<0){
                    bad.cor[1,names(cor)[cr]]=mean(rtrunc(10000,spec="beta",a=mu2,b=1,shape1=a2,shape2=b2))
                    good.cor[1,names(cor)[cr]]=mean(rtrunc(10000,spec="beta",a=0,b=mu2,shape1=a2,shape2=b2))
                    
                  }else{
                    bad.cor[1,names(cor)[cr]]=mean(rtrunc(10000,spec="beta",a=0,b=mu2,shape1=a2,shape2=b2))
                    good.cor[1,names(cor)[cr]]=mean(rtrunc(10000,spec="beta",a=mu2,b=1,shape1=a2,shape2=b2))
                    
                    
                  }
                  
                  
                }
                
                bad.cor[is.na(bad.cor)]=1
                good.cor[is.na(good.cor)]=1
                
                # put new values for average good and bad conditions into a vector
                mu.vec.good=c(good.cor[1,"s3"],good.cor[1,"s4"],good.cor[1,"g43"],good.cor[1,"g34"],good.mean[,"f13"],good.mean[,"f23"],good.mean[,"f33"]) 
                mu.vec.bad=c(bad.cor[1,"s3"],bad.cor[1,"s4"],bad.cor[1,"g43"],bad.cor[1,"g34"],bad.mean[,"f13"],bad.mean[,"f23"],bad.mean[,"f33"])
                
                }else{
                
                  # If fecundity is 0, fill the vectors simply with average values  
                  mu.vec.good=mu.vec.bad=c(vr.mean[,"s3"],vr.mean[,"s4"],vr.mean[,"g43"],vr.mean[,"g34"],vr.mean[,"f13"],vr.mean[,"f23"],vr.mean[,"f33"]) 
                  
              }
              
              
              
            }else{ 
              
              # For survival/growth, make variation dependent on maximum CV of a binomial variable
              # and use a beta distribution to simulate values of vital rates 
              
              maxCV=sqrt(dat[,y]*(1-dat[,y]))/dat[,y]
              
              if(is.na(maxCV)) maxCV = 0
              
              if(maxCV>0&!is.null(sigma.sub)){
                
                # simulate from a beta distribution (FOR SURVIVAL)
                mu=dat[,y]
                var=(cv[i]*maxCV*dat[,y])^2
                
                # convert mean and variance to beta parameters
                shape= mu*((mu*(1-mu))/var-1)
                scale=(1-mu)*((mu*(1-mu))/var-1)
                
                bad.mean[,y]=mean(rtrunc(10000,spec="beta",a=0,b=mu,shape1=shape,shape2=scale)) # for the vital rates that are 0, this will remain 0
                good.mean[,y]=mean(rtrunc(10000,spec="beta",a=mu,b=1,shape1=shape,shape2=scale))  # for the vital rates that are 0, this will remain 0
                
                ### Correlation
                cor=sigma.sub[rownames(sigma.sub)[rownames(sigma.sub)!=colnames(vr.cor)[y]],colnames(vr.cor)[y]]
                
                names(cor)=rownames(sigma.sub)[rownames(sigma.sub)!=colnames(vr.cor)[y]]
                cor=cor[cor>(abs(min(sigma.sub))-0.1)|cor<(min(sigma.sub)+0.1)]
                
                good.cor=bad.cor=matrix(as.numeric(dat[,5:7]),1,3)
                
                colnames(good.cor)=colnames(bad.cor)=colnames(vr.cor)[5:7]
                
                for(cr in 1:length(cor)){
                  
                  # simulate from a gamma distribution (FOR VITAL RATES CORRELATED WITH SURVIVAL)
                  mu2=dat[,names(cor)[cr]]
                  var2=abs(cor[cr])*(cv[i]*mu2)^2
                  
                  # convert mean and variance to gamma parameters
                  a2= mu2^2/var2
                  b2=var2/mu2
                  
                  
                  if(cor[cr]<0){
                    bad.cor[1,names(cor)[cr]]=mean(rtrunc(10000,spec="gamma",a=mu2,b=Inf,shape=a2,scale=b2))
                    good.cor[1,names(cor)[cr]]=mean(rtrunc(10000,spec="gamma",a=0,b=mu2,shape=a2,scale=b2))
                    
                  }else{
                    bad.cor[1,names(cor)[cr]]=mean(rtrunc(10000,spec="gamma",a=0,b=mu2,shape=a2,scale=b2))
                    good.cor[1,names(cor)[cr]]=mean(rtrunc(10000,spec="gamma",a=mu2,b=Inf,shape=a2,scale=b2))
                    
                    
                  }
                  
                  
                }
                
                bad.cor[is.na(bad.cor)]=0.001
                good.cor[is.na(good.cor)]=0.001
                
                mu.vec.good=c(good.mean[,"s3"],good.mean[,"s4"],good.mean[,"g43"],good.mean[,"g34"],good.cor[1,"f13"],good.cor[1,"f23"],good.cor[1,"f33"]) 
                mu.vec.bad=c(bad.mean[,"s3"],bad.mean[,"s4"],bad.mean[,"g43"],bad.mean[,"g34"],bad.cor[1,"f13"],bad.cor[1,"f23"],bad.cor[1,"f33"]) 

              }else{
                
                # if survival/growth = 0
                mu.vec.good=mu.vec.bad=c(vr.mean[,"s3"],vr.mean[,"s4"],vr.mean[,"g43"],vr.mean[,"g34"],vr.mean[,"f13"],vr.mean[,"f23"],vr.mean[,"f33"])  
                
              }
              
              
            } 
            
            ## SIMULATE WITHIN-STATE VARIANCE (USING COPULAS)
            
            if(!is.null(nrow(sigma.sub))){
              z.m <- mvrnorm(100,mu=rep(0,nrow(sigma.sub)),Sigma=sigma.sub)
              #Apply the Normal CDF function to Z to obtain data that is uniform 
              # on the interval [0,1], but still correlated.
              
              u.m <- pnorm(z.m)
              
              if(length(grep("s3",colnames(sigma.sub)))>0){
                
                if(dat[,"s3"]==1){
                  mu.g=0.999
                  mu.b=0.999}else{
                    mu.g=mu.vec.good[1]
                    mu.b=mu.vec.bad[1]}
                # good
                maxCV.g=sqrt(mu.g*(1-mu.g))/mu.g
                var.g=(0.01*maxCV.g*mu.g)^2
                a.g= mu.g*((mu.g*(1-mu.g)/var.g)-1)
                b.g=(1-mu.g)*((mu.g*(1-mu.g))/var.g-1)
                
                s3.g=qbeta(u.m[,which(colnames(sigma.sub)=="s3")],shape1=a.g,
                           shape2=b.g)
                #bad
                maxCV.b=sqrt(mu.b*(1-mu.b))/mu.b
                var.b=(0.01*maxCV.b*mu.b)^2
                a.b= mu.b*((mu.b*(1-mu.b)/var.b)-1)
                b.b=(1-mu.b)*((mu.b*(1-mu.b))/var.b-1)
                
                s3.b=qbeta(u.m[,which(colnames(sigma.sub)=="s3")],shape1=a.b,
                           shape2=b.b)
              }else{
                s3.g=0
                s3.b=0}
              
              ### S4
              if(length(grep("s4",colnames(sigma.sub)))>0){
                
                if(dat[,"s4"]==1){
                  mu.g=0.999
                  mu.b=0.999}else{
                    mu.g=mu.vec.good[2]
                    mu.b=mu.vec.bad[2]}
                # good
                maxCV.g=sqrt(mu.g*(1-mu.g))/mu.g
                var.g=(0.01*maxCV.g*mu.g)^2
                a.g= mu.g*((mu.g*(1-mu.g)/var.g)-1)
                b.g=(1-mu.g)*((mu.g*(1-mu.g))/var.g-1)
                
                s4.g=qbeta(u.m[,which(colnames(sigma.sub)=="s4")],shape1=a.g,
                           shape2=b.g)
                #bad
                maxCV.b=sqrt(mu.b*(1-mu.b))/mu.b
                var.b=(0.01*maxCV.b*mu.b)^2
                a.b= mu.b*((mu.b*(1-mu.b)/var.b)-1)
                b.b=(1-mu.b)*((mu.b*(1-mu.b))/var.b-1)
                
                s4.b=qbeta(u.m[,which(colnames(sigma.sub)=="s4")],shape1=a.b,
                           shape2=b.b)
              }else{
                s4.g=0
                s4.b=0}
              
              ### g43
              if(length(grep("g43",colnames(sigma.sub)))>0){
                
                if(dat[,"g43"]==1){
                  mu.g=0.999
                  mu.b=0.999}else{
                    mu.g=mu.vec.good[3]
                    mu.b=mu.vec.bad[3]}
                # good
                maxCV.g=sqrt(mu.g*(1-mu.g))/mu.g
                var.g=(0.01*maxCV.g*mu.g)^2
                a.g= mu.g*((mu.g*(1-mu.g)/var.g)-1)
                b.g=(1-mu.g)*((mu.g*(1-mu.g))/var.g-1)
                
                g43.g=qbeta(u.m[,which(colnames(sigma.sub)=="g43")],shape1=a.g,
                            shape2=b.g)
                #bad
                maxCV.b=sqrt(mu.b*(1-mu.b))/mu.b
                var.b=(0.01*maxCV.b*mu.b)^2
                a.b= mu.b*((mu.b*(1-mu.b)/var.b)-1)
                b.b=(1-mu.b)*((mu.b*(1-mu.b))/var.b-1)
                
                g43.b=qbeta(u.m[,which(colnames(sigma.sub)=="g43")],shape1=a.b,
                            shape2=b.b)
              }else{
                g43.g=0
                g43.b=0}
              
              ### g34
              if(length(grep("g34",colnames(sigma.sub)))>0){
                
                if(dat[,"g34"]==1){
                  mu.g=0.999
                  mu.b=0.999}else{
                    mu.g=mu.vec.good[4]
                    mu.b=mu.vec.bad[4]}
                # good
                maxCV.g=sqrt(mu.g*(1-mu.g))/mu.g
                var.g=(0.01*maxCV.g*mu.g)^2
                a.g= mu.g*((mu.g*(1-mu.g)/var.g)-1)
                b.g=(1-mu.g)*((mu.g*(1-mu.g))/var.g-1)
                
                g34.g=qbeta(u.m[,which(colnames(sigma.sub)=="g34")],shape1=a.g,
                            shape2=b.g)
                #bad
                maxCV.b=sqrt(mu.b*(1-mu.b))/mu.b
                var.b=(0.01*maxCV.b*mu.b)^2
                a.b= mu.b*((mu.b*(1-mu.b)/var.b)-1)
                b.b=(1-mu.b)*((mu.b*(1-mu.b))/var.b-1)
                
                g34.b=qbeta(u.m[,which(colnames(sigma.sub)=="g34")],shape1=a.b,
                            shape2=b.b)
              }else{
                g34.g=0
                g34.b=0}
              
              ### f13
              if(length(grep("f13",colnames(sigma.sub)))>0){
                
                mu.g=mu.vec.good[5]
                mu.b=mu.vec.bad[5]
                
                # good
                var.gg=(0.01*mu.g)^2
                a.g= mu.g^2/var.gg
                b.g=var.gg/mu.g
                
                f13.g=qgamma(u.m[,which(colnames(sigma.sub)=="f13")],shape=a.g,
                             scale=b.g)
                #bad
                var.bb=(0.01*mu.b)^2
                a.b= mu.b^2/var.bb
                b.b=var.bb/mu.b
                
                f13.b=qgamma(u.m[,which(colnames(sigma.sub)=="f13")],shape=a.b,
                             scale=b.b)
              }else{
                f13.g=0
                f13.b=0}
              
              ### f23
              if(length(grep("f23",colnames(sigma.sub)))>0){
                
                mu.g=mu.vec.good[6]
                mu.b=mu.vec.bad[6]
                # good
                var.gg=(0.01*mu.g)^2
                a.g= mu.g^2/var.gg
                b.g=var.gg/mu.g
                
                f23.g=qgamma(u.m[,which(colnames(sigma.sub)=="f23")],shape=a.g,
                             scale=b.g)
                #bad
                var.bb=(0.01*mu.b)^2
                a.b= mu.b^2/var.bb
                b.b=var.bb/mu.b
                
                f23.b=qgamma(u.m[,which(colnames(sigma.sub)=="f23")],shape=a.b,
                             scale=b.b)
              }else{
                f23.g=0
                f23.b=0}
              
              ### f33
              if(length(grep("f33",colnames(sigma.sub)))>0){
                
                mu.g=mu.vec.good[7]
                mu.b=mu.vec.bad[7]
                # good
                var.gg=(0.01*mu.g)^2
                a.g= mu.g^2/var.gg
                b.g=var.gg/mu.g
                
                f33.g=qgamma(u.m[,which(colnames(sigma.sub)=="f33")],shape=a.g,
                             scale=b.g)
                #bad
                var.bb=(0.01*mu.b)^2
                a.b= mu.b^2/var.bb
                b.b=var.bb/mu.b
                
                f33.b=qgamma(u.m[,which(colnames(sigma.sub)=="f33")],shape=a.b,
                             scale=b.b)
              }else{
                f33.g=0
                f33.b=0}
              
            }else{
              s3.g=s4.g=g43.g=g34.g=f13.g=f23.g=f33.g=rep(0,100) 
              s3.b=s4.b=g43.b=g34.b=f13.b=f23.b=f33.b=rep(0,100) 
              
            }
            
            
            #### CREATE AN ARRAY OF GOOD AND BAD ENVIRONMENTAL STATES
            ##### PUT ALL TOGETHER
            
            dat.vec.g=cbind(vr.mean[,c("s1","s2","g21","g31","g32")],s3.g,s4.g,g43.g,g34.g,f13.g,f23.g,f33.g)
            dat.vec.b=cbind(vr.mean[,c("s1","s2","g21","g31","g32")],s3.b,s4.b,g43.b,g34.b,f13.b,f23.b,f33.b)
            
            colnames(dat.vec.g)=colnames(dat.vec.b)=c("s1","s2","g21","g31","g32","s3","s4","g43","g34","f13","f23","f33")
            
            good.bad.array=array(0,c(4,4,100,2))
            
            for(dist in 1:100){
              
              good.bad.array[,,dist,1]=matrix(c(dat.vec.g[dist,"s1"]*(1-dat.vec.g[dist,"g21"]-dat.vec.g[dist,"g31"]-0),dat.vec.g[dist,"s1"]*dat.vec.g[dist,"g21"],dat.vec.g[dist,"s1"]*dat.vec.g[dist,"g31"],0,
                                                0,dat.vec.g[dist,"s2"]*(1-0-dat.vec.g[dist,"g32"]-0),dat.vec.g[dist,"s2"]*dat.vec.g[dist,"g32"],0,
                                                dat.vec.g[dist,"f13"],dat.vec.g[dist,"f23"],dat.vec.g[dist,"f33"]+dat.vec.g[dist,"s3"]*(1-dat.vec.g[dist,"g43"]),dat.vec.g[dist,"s3"]*dat.vec.g[dist,"g43"],
                                                0,0,dat.vec.g[dist,"s4"]*dat.vec.g[dist,"g34"],dat.vec.g[dist,"s4"]*(1-dat.vec.g[dist,"g34"])),
                                              4,4,byrow=F)
              
              good.bad.array[,,dist,2]=matrix(c(dat.vec.b[dist,"s1"]*(1-dat.vec.b[dist,"g21"]-dat.vec.b[dist,"g31"]-0),dat.vec.b[dist,"s1"]*dat.vec.b[dist,"g21"],dat.vec.b[dist,"s1"]*dat.vec.b[dist,"g31"],0,
                                                0,dat.vec.b[dist,"s2"]*(1-0-dat.vec.b[dist,"g32"]-0),dat.vec.b[dist,"s2"]*dat.vec.b[dist,"g32"],0,
                                                dat.vec.b[dist,"f13"],dat.vec.b[dist,"f23"],dat.vec.b[dist,"f33"]+dat.vec.b[dist,"s3"]*(1-dat.vec.b[dist,"g43"]),dat.vec.b[dist,"s3"]*dat.vec.b[dist,"g43"],
                                                0,0,dat.vec.b[dist,"s4"]*dat.vec.b[dist,"g34"],dat.vec.b[dist,"s4"]*(1-dat.vec.b[dist,"g34"])),
                                              4,4,byrow=F)
            }
            
            # Create environmental transition matrix P (working with f and v1)
            
            q=f[z]*(1-v1[j])
            
            p=v1[j]+q
            
            
            P=matrix(c(p,1-p,
                       q,1-q),nrow=2,ncol=2,byrow=F)
            colnames(P) <- c("1","2")
            row.names(P) <- c("1","2")
            
            # vector of environments
            simul.n = numeric(trun+1)
            
            # get sequence of environments
            # start with god environment 
            env_at_t=1
            simul.n[1]=1
            
            for(sim_t in 2:trun)
              simul.n[sim_t] = env_at_t =sample(ncol(P),1,pr=P[,env_at_t])
            
            
            # Calculate stochastic lambda
            ########################################################
            nstage <- dim(good.bad.array)[1] 
            states <- simul.n
            states.within <- sample(1:100,trun+1,replace=T)
            growth <- array(0,ts)   
            
            # Initial population vector
            vec1=rep(1,nstage)
            vec1 <- vec1/sum(vec1)
            vec1 <- t(vec1) 
            
            # ITERATION TO CALCULATE LAMBDA FOR EACH TIME STEP
            for (a  in 1:trun){
              i2<-states.within[a]
              i3 <- states[a]
              mat1 <-  good.bad.array[,,i2,i3]
              vec1 <- mat1%*%as.numeric(vec1)
              growth1 <- sum(vec1) # population growth at one time step 
              vec1 <- vec1/growth1
              if( a > tr){ # after the burn-in, save the growth rate for each time step      
                i1 <- a - tr
                growth[i1] <- growth1
              }
              
            }
            
            a1=sum(log(growth[1:ts]))
            
            lambda.s[sc,x,y,j,i,z]= a1/ts  
            
            
          }
          
        }
      }
    }
  }
}




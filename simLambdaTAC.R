# Script for Paniw et al. XXXXXX - Appendix S3 - Simulations of temporal autocorrelation in vital rates 

# This script is divided into two parts:
# PART A (lines 33-542): Simulate temporal autocorrelation and calculate sensitivity of the stochastic growth rate to changes in correlation structure
# PART B (lines 543-843): Regress Sv1 (sensitivity to temporal autocorrelation) against PCA scores, vital-rate classes, and habitat types

# Author: Maria Paniw

########## 
# Clean memory
rm(list=ls(all=TRUE))

### load libraries

library(ggplot2)
library(reshape2)
library(MASS)
library(Cairo)
library(plyr)
library(pbkrtest)
library(dplyr)
# set working directory 

setwd("/Users/mariapaniw/Dropbox/TempAutoProject/SuppMat")

#####################################################################################
######################################################################
############# PART A - SIMULATIONS

load("matsMean")

# Create some parameters necessary for the simulation (loop)

# Autocorrelation:

v1=c(-0.3,0,0.3)

# Coefficient of variation

cv=c(0.2,0.5,0.8)

### Define how long simulations run and what the characteristics of the environments are

# Define simulation time (to run faster simulations, the user may decrease both tr and ts)
tr=20000 #the discard time
ts=5000 # how many years to keep

# run each simulation for ts+tr years 
trun=ts+tr

# The different environments 
env=c("good","bad")
env.n=1:2

# long-term frequency of good condition = first value of right eigenvector describing 
# stationary distribution of environmental states (see Tuljapurkar & Haridas, Ecol Lett, 2006 for more detail)

f=c(0.35,0.65)

rep.simul=5
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


# The output array will have 5 dimensions: species x vital rates x v1 x cv x f

sp.sub=1:length(mats) # the user may wish to work with a smaller subset of species to speed up simulations

lambda.s=vector("list", length(sp.sub))

### ACTUAL SIMULATIONS

# First main loop: ppopulations

for(x in 1:length(sp.sub)){
  
  vr=mats[[sp.sub[x]]]$vr[mats[[sp.sub[x]]]$vr>0]
  if(length(which(vr[-grep("f",names(vr))]==1))>0) vr=vr[-which(vr[-grep("f",names(vr))]==1)]
  
  Umat=mats[[sp.sub[x]]]$matU
  Fmat=mats[[sp.sub[x]]]$matF
  lambda.s.temp=array(NA,c(length(vr),length(v1),length(cv),length(f),rep.simul))
  var.temp=array(NA,c(length(vr),length(v1),length(cv),length(f),rep.simul))
  
  # Second main loop: Replicate simulations 
  
  for(sm in 1:rep.simul){
    
    number.sim=seq(100,149)[sm]
    
    # Third main loop: frequency of good environment  
    for(z in 1:length(f)){
      
      # Fourth main loop: coefficient of variation 
      for(i in 1:length(cv)){
        
        # Fifth main loop: autocorelation coefficient
        
        for(j in 1:length(v1)){
          
          # Sixth main loop: vital rates
          
          
          for(y in 1:length(vr)){
            
            #############################################
            ###############################
            # perturb the vr data frame to create bad and good condition mean vital rates
            
            bad.mean=vr
            
            good.mean=vr
            
            vr.mean=vr
            
            good.vec=bad.vec=numeric(100)
            if(y>=min(grep("f",names(vr)))){ # For fecundity, keep regular CV to define variation
              
              
              # simulate from a gamma distribution
              mu=vr.mean[y]
              var=(cv[i]*vr.mean[y])^2
              
              # convert mean and variance to gamma parameters
              shape= mu^2/var
              scale=var/mu
              
              # new means for good and bad are the averages of two truncated gamma distributions
              # both with a mean given by the vital rate at its average across environments 
              # and variance given by the CV 
              set.seed(number.sim)
              bad.mean[y]=mean(rtrunc(10000,spec="gamma",a=0,b=mu,shape=shape,scale=scale)) # upper bound truncated at the mean 
              
              set.seed(number.sim)
              good.mean[y]=mean(rtrunc(10000,spec="gamma",a=mu,b=Inf,shape=shape,scale=scale)) # lower bound truncated at the mean
              
              # good and bad vectors:
              
              # within each good and bad environment,
              # 100 random values of a vital rate are produced 
              # from a gamma distribution using the average values for new means of bad/good environments
              # and very low coefficient of variation = 0.01
              var.g=(0.01*good.mean[y])^2
              var.b=(0.01*bad.mean[y])^2
              
              set.seed(number.sim)
              good.vec=rgamma(100,shape=good.mean[y]^2/var.g,scale=var.g/good.mean[y])
              
              set.seed(number.sim)
              bad.vec=rgamma(100,shape=bad.mean[y]^2/var.b,scale=var.b/bad.mean[y])
              
              
              
              
              
            }else if(y<min(grep("f",names(vr)))){ 
              
              # For survival/growth, make variation dependent on maximum CV of a binomial variable
              # and use a beta distribution to simulate values of vital rates 
              
              maxCV=sqrt(vr.mean[y]*(1-vr.mean[y]))/vr.mean[y]
              mu=vr.mean[y]
              var=(cv[i]*maxCV*vr.mean[y])^2
              # convert mean and variance to beta parameters
              shape= mu*((mu*(1-mu))/var-1)
              scale=(1-mu)*((mu*(1-mu))/var-1)
              
              set.seed(number.sim)
              bad.mean[y]=mean(rtrunc(10000,spec="beta",a=0,b=mu,shape1=shape,shape2=scale)) # for the vital rates that are 0, this will remain 0
              set.seed(number.sim)
              good.mean[y]=mean(rtrunc(10000,spec="beta",a=mu,b=1,shape1=shape,shape2=scale))  # for the vital rates that are 0, this will remain 0
              
              # good and bad vectors:
              
              # define a new maximum variance
              maxCV.g=sqrt(good.mean[y]*(1-good.mean[y]))/good.mean[y]
              maxCV.b=sqrt(bad.mean[y]*(1-bad.mean[y]))/bad.mean[y]
              
              # new coefficient of variation
              var.g=(0.01*maxCV.g*good.mean[y])^2
              var.b=(0.01*maxCV.b*bad.mean[y])^2
              
              shape.g= good.mean[y]*((good.mean[y]*(1-good.mean[y]))/var.g-1)
              scale.g=(1-good.mean[y])*((good.mean[y]*(1-good.mean[y]))/var.g-1)
              
              shape.b= bad.mean[y]*((bad.mean[y]*(1-bad.mean[y]))/var.b-1)
              scale.b=(1-bad.mean[y])*((bad.mean[y]*(1-bad.mean[y]))/var.b-1)
              
              set.seed(number.sim)
              good.vec=rbeta(100,shape1=shape.g,shape2=scale.g)
              good.vec[!is.finite(good.vec)]=0
              
              set.seed(number.sim)
              bad.vec=rbeta(100,shape1=shape.b,shape2=scale.b)
              bad.vec[!is.finite(bad.vec)]=0
              
              
            } 
            
            # Constarin survival and growth to not go below 0 or above 1 (just in case)
            
            if(y<min(grep("f",names(vr)))){
              
              good.vec[good.vec>1]=1
              
            }
            
            bad.vec[bad.vec<0]=0
            
            #############################################
            ###############################
            
            # Create good and bad condition matrices from the perturbed and non-perturbed vital rates
            
            ### Create a new matrix that has 100 rows (simulated values for a give vital rate) and 12 columns(vr)
            
            ## Array with 1 matrix per row in dataframe
            good.bad.array=array(0,c(dim(Umat)[1],dim(Umat)[1],100,2))
            mean.mat=Umat+Fmat
            
            for(dist in 1:100){
              
              # survival
              if(y<=max(grep("s",names(vr)))){
                if(good.vec[dist]>0|bad.vec[dist]>0){
                  
                  surv=surv.g=surv.b=colSums(Umat)
                  surv.g[y]=good.vec[dist]
                  surv.b[y]=bad.vec[dist]
                  
                  #take out survival from suvival-trnasition matrix
                  U.mat.g=U.mat.b=Umat
                  for(xx in 1:ncol(U.mat.g)){
                    
                    U.mat.g[,xx]=Umat[,xx]/surv[xx]
                    U.mat.g[,xx]=U.mat.g[,xx]*surv.g[xx]
                    U.mat.b[,xx]=Umat[,xx]/surv[xx]
                    U.mat.b[,xx]=U.mat.b[,xx]*surv.b[xx]
                  }
                  
                  U.mat.g[!is.finite(U.mat.g)]=0
                  U.mat.b[!is.finite(U.mat.b)]=0
                  
                  
                  good.bad.array[,,dist,1]=U.mat.g+Fmat
                  good.bad.array[,,dist,2]=U.mat.b+Fmat
                  
                }else{
                  
                  good.bad.array[,,dist,1]=Umat+Fmat
                  good.bad.array[,,dist,2]=Umat+Fmat
                }
                
                # Transitions (progression & retrogression)
              }else if(y>max(grep("s",names(vr)))&y<min(grep("f",names(vr)))){
                
                # Control: after changing values for growth and retrogresion, columns Umat excluding survival must sum to 1
                # This is then technically a multinomial distribution
                
                if(good.vec[dist]>0&bad.vec[dist]>0){
                  
                  all=which((substring(names(vr), 1, 1)=="g"|substring(names(vr), 1, 1)=="r")&substring(names(vr), 4, nchar(names(vr)))==substring(names(vr[y]), 4, nchar(names(vr[y]))))
                  
                  col.now=as.numeric(substring(names(vr)[y], 4, nchar(names(vr[y]))))
                  
                  
                  vr.all.good=vr[all[!all%in%y]]
                  vr.all.bad=vr[all[!all%in%y]]
                  vr.all2.g=vr[all]
                  vr.all2.b=vr[all]
                  vr.all2.g[names(vr.all2.g)%in%names(vr)[y]]=good.vec[dist]
                  vr.all2.b[names(vr.all2.b)%in%names(vr)[y]]=bad.vec[dist]
                  
                  # If no stasis:
                  if(diag(mean.mat)[col.now]==0){
                    stasis.g=0
                    stasis.b=0
                    
                    if((1-(good.vec[dist]+sum(vr.all.good)))<0){
                      
                      vr.all.good[vr.all.good-abs(1-(good.vec[dist]+sum(vr.all.good)))>0]=vr.all.good[vr.all.good-abs(1-(good.vec[dist]+sum(vr.all.good)))>0]-abs(1-(good.vec[dist]+sum(vr.all.good)))/length(vr.all.good[vr.all.good-abs(1-(good.vec[dist]+sum(vr.all.good)))>0])
                    }else if((1-(bad.vec[dist]+sum(vr.all.bad)))>0){
                      
                      vr.all.bad=vr.all.bad+abs(1-(bad.vec[dist]+sum(vr.all.bad)))/length(vr.all.bad)
                    }
                  }else{
                    
                    # If stasis
                    stasis.b=1-(bad.vec[dist]+sum(vr.all.bad))
                    stasis.g=1-(good.vec[dist]+sum(vr.all.good))
                    if(stasis.g<0){
                      vr.all.good[vr.all.good-abs(stasis.g)>0]=vr.all.good[vr.all.good-abs(stasis.g)>0]-abs(stasis.g)/length(vr.all.good[vr.all.good-abs(stasis.g)>0])
                      stasis.g=1-(good.vec[dist]+sum(vr.all.good))
                    }
                  }
                  
                  
                  vr.all2.g[names(vr.all2.g)%in%names(vr.all.good)]=vr.all.good
                  vr.all2.b[names(vr.all2.b)%in%names(vr.all.bad)]=vr.all.bad
                  
                  surv=colSums(Umat)
                  U.mat.gr.g=U.mat.gr.b=Umat
                  for(xx in 1:ncol(U.mat.gr.g)){
                    U.mat.gr.g[,xx]=Umat[,xx]/surv[xx]
                    U.mat.gr.b[,xx]=Umat[,xx]/surv[xx]
                    
                  }
                  
                  U.mat.gr.g[!is.finite(U.mat.gr.g)]=0
                  U.mat.gr.b[!is.finite(U.mat.gr.b)]=0
                  
                  # Fill transition matrix with new values
                  diag(U.mat.gr.g)[col.now] <- stasis.g
                  diag(U.mat.gr.b)[col.now] <- stasis.b
                  
                  U.mat.gr.g[,col.now][lower.tri(U.mat.gr.g)[,col.now]][U.mat.gr.g[,col.now][lower.tri(U.mat.gr.g)[,col.now]]>0]=vr.all2.g[grep("g",names(vr.all2.g))]
                  U.mat.gr.g[,col.now][upper.tri(U.mat.gr.g)[,col.now]][U.mat.gr.g[,col.now][upper.tri(U.mat.gr.g)[,col.now]]>0]=vr.all2.g[grep("r",names(vr.all2.g))] 
                  
                  U.mat.gr.b[,col.now][lower.tri(U.mat.gr.b)[,col.now]][U.mat.gr.b[,col.now][lower.tri(U.mat.gr.b)[,col.now]]>0]=vr.all2.b[grep("g",names(vr.all2.b))]
                  U.mat.gr.b[,col.now][upper.tri(U.mat.gr.b)[,col.now]][U.mat.gr.b[,col.now][upper.tri(U.mat.gr.b)[,col.now]]>0]=vr.all2.b[grep("r",names(vr.all2.b))] 
                  
                  # Make sure suvival-dependent transition matrix sums to 1
                  for(xx in 1:ncol(U.mat.gr.b)){
                    
                    U.mat.gr.g[,xx][ U.mat.gr.g[,xx]<0]=0
                    U.mat.gr.g[,xx][ U.mat.gr.g[,xx]>0]=U.mat.gr.g[,xx][ U.mat.gr.g[,xx]>0]/sum(U.mat.gr.g[,xx][ U.mat.gr.g[,xx]>0])
                    U.mat.gr.b[,xx][ U.mat.gr.b[,xx]<0]=0
                    U.mat.gr.b[,xx][ U.mat.gr.b[,xx]>0]=U.mat.gr.b[,xx][ U.mat.gr.b[,xx]>0]/sum(U.mat.gr.b[,xx][ U.mat.gr.b[,xx]>0])
                    
                    U.mat.gr.g[,xx]=U.mat.gr.g[,xx]*surv[xx]
                    
                    U.mat.gr.b[,xx]=U.mat.gr.b[,xx]*surv[xx]
                    
                  }
                  
                  good.bad.array[,,dist,1]=U.mat.gr.g+Fmat
                  good.bad.array[,,dist,2]=U.mat.gr.b+Fmat
                  
                }else{
                  
                  good.bad.array[,,dist,1]=Umat+Fmat
                  good.bad.array[,,dist,2]=Umat+Fmat
                }
                
                # Reproduction 
              }else if(y>=min(grep("f",names(vr)))){
                
                F.mat.g=F.mat.b=Fmat
                F.mat.g[as.numeric(substring(names(vr)[y], 2, 3)),as.numeric(substring(names(vr)[y], 4, nchar(names(vr[y]))))]=good.vec[dist]
                F.mat.b[as.numeric(substring(names(vr)[y], 2, 3)),as.numeric(substring(names(vr)[y], 4, nchar(names(vr[y]))))]=bad.vec[dist]
                
                good.bad.array[,,dist,1]=Umat+F.mat.g
                good.bad.array[,,dist,2]=Umat+F.mat.b
              }
            }
            
            # Create environmental transition matrix P (working with f and v1)
            
            q=f[z]*(1-v1[j])
            
            p=v1[j]+q
            
            
            P=matrix(c(p,1-p,
                       q,1-q),nrow=2,ncol=2,byrow=F)
            colnames(P) <- c("1","2")
            row.names(P) <- c("1","2")
            
            #simulations for variance
            vr.sim.result=rep(0,50)
            for(vr.sim in 1:length(vr.sim.result)){
              
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
              
              vr.sim.result[vr.sim]=sum(log(growth[1:ts])) 
            }
            
            # mean of stochastic lambda across simulations
            a.hat=sum(vr.sim.result)/(length(vr.sim.result)*ts)
            
            # variance of stochastic lambda
            var.hat=sum((vr.sim.result-a.hat*ts)^2)/(length(vr.sim.result)*ts)
            lambda.s.temp[y,j,i,z,sm]= a.hat
            var.temp[y,j,i,z,sm]= var.hat
            
          }
          
        }
      }
    }
  }
  
  
  
  lambda.s[[x]]$lambda.s=lambda.s.temp
  lambda.s[[x]]$lambda.s.var=var.temp
  lambda.s[[x]]$names.vr=names(vr)
}


##############################################################################
################################# Put results into data frame
sens=NULL
lambda.final=NULL
for(x in 1:length(mats)){
  
  lambda=lambda.s[[x]]$lambda.s
  
  ### SV1 (sensitivities) output
  
  #lambda dims are for vr, v1 (here 3),cv, f, and sim runs 
  sens.temp=array(0,c(dim(lambda)[1],dim(lambda)[3],dim(lambda)[4],dim(lambda)[5]))
  
  for(xx in 1:dim(lambda)[5]){
    
    for(z in 1:dim(lambda)[4]){
      
      for(a in 1:dim(lambda)[3]){
        
        sens.temp[,a,z,xx]=(lambda[,3,a,z,xx]-lambda[,2,a,z,xx])/(abs(v1[3]-v1[2]))
        
      }
    }
  }
  
  sens.temp2=apply(sens.temp,c(1,2,3),mean) # average Sv1 across simualtion runs
  
  data=adply(sens.temp2,c(1,2,3))
  colnames(data)=c("vr","cv","f","sens")
  
  data$var.sens=adply(apply(sens.temp,c(1,2,3),var),c(1,2,3))[,4]
  
  data$species=mats[[x]]$species
  data$matDim=dim(mats[[x]]$matU)[1]
  levels(data$vr)=lambda.s[[x]]$names.vr
  levels(data$cv)=cv
  levels(data$f)=f
  
  class=mats[[x]]$class
  vr.name=lambda.s[[x]]$names.vr
  
  data$class=NA
  for(xx in 1:length(vr.name)){
    if(length(grep("s",vr.name[xx]))>0){
      
      data$class[data$vr==vr.name[xx]]=paste("s",class[as.numeric(substring(vr.name[xx], 2, 2))],sep="")
      
    }else if(length(grep("g",vr.name[xx]))>0|length(grep("r",vr.name[xx]))>0){
      
      data$class[data$vr==vr.name[xx]]=paste(substring(vr.name[xx], 1, 1),class[as.numeric(substring(vr.name[xx], 4, 5))],sep="")
    }else{
      
      data$class[data$vr==vr.name[xx]]=paste(substring(vr.name[xx], 1, 1),class[as.numeric(substring(vr.name[xx], 2, 3))],sep="")
      
    }
    
  }
  sens=rbind(sens,data)
  
  ### LAMBDA OUTPUT
  
  # mean lambda across simulation runs
  lambda=apply(lambda,c(1,2,3,4),mean)
  lambda.var=lambda.s[[x]]$lambda.s.var
  
  # mean variance across simulation runs
  lambda.var=apply(lambda.var,c(1,2,3,4),mean)
  
  lambda.sub=adply(lambda,c(1,2,3,4))
  colnames(lambda.sub)=c("vr","v1","cv","f","lambda")
  lambda.sub$var=adply(lambda.var,c(1,2,3,4))[,5]
  
  lambda.sub$species=mats[[x]]$species
  lambda.sub$matDim=dim(mats[[x]]$matU)[1]
  levels(lambda.sub$vr)=lambda.s[[x]]$names.vr
  levels(lambda.sub$cv)=cv
  levels(lambda.sub$f)=f
  levels(lambda.sub$v1)=v1
  
  lambda.final=rbind(lambda.final,lambda.sub)
  
}
### add PCA scores

pca.scores=read.csv("PCAscores.csv")

sens$PC1=left_join(sens,PCA_scores,by="species")$PC1
sens$PC2=left_join(sens,PCA_scores,by="species")$PC2

##################################################################################
################################# 
############# # PART B - GAMs (sensitivity to temporal autocorrelation) against PCA scores

library(mgcv)

##### Running the simulations on all populations is time-comsuming
# but some models may not converge if simulations were done on a smaller
# dataset that potentially excluded vital rates

# We therefore provide the data file Sv1_PCA that includes sensitivities 
# for all the vital rates and 455 populations 

data=read.csv("Sv1_PCA.csv")

data$sens=abs(data$sens)
data=data[data$sens>0,]

data$cv=as.factor(data$cv)
data$f=as.factor(data$f)

data$class=as.character(data$class)
data$class=gsub("dorm","NR",data$class)

data$class=gsub("prop","P",data$class)

data$class=factor(data$class,levels=c("sP","sPR","sR","sNR","gP","gPR","gR","gNR","rPR","rR","rNR","fP","fPR","fR","fNR"))

# Put all transitions (progression & retrogression) together

levels(data$class)[c(9,10,11)]=c("gPR","gR","gNR")

data$cv=factor(data$cv)
data$f=factor(data$f)

########## STEP 1: Sv1 ~ CV + f
dat=data[data$sens>0,]

dat=aggregate(sens~cv+f+class+species,data=dat,sum)

dat$sens=log(dat$sens)
dat=dat[is.finite(dat$sens),]
dat=left_join(dat,data[,c("species","PC1","PC2","habitat","matDim")],by=c("species"))
dat=dat[-which(duplicated(dat)),]

dat=dat[!is.na(dat$PC1),]

dat$species=factor(dat$species)
dat$matDim=factor(dat$matDim)
levels(dat$matDim)[8:28]="10up"
levels(dat$matDim)[1:2]="4low"

################### STEP 1 - SENS ~ CV + F
mod1=gam(sens ~ s(matDim,bs="re",k=4)+s(species,bs="re",k=4),data=dat,gamma=1.4)
mod2=gam(sens ~ s(matDim,bs="re",k=4)+s(species,bs="re",k=4) +cv,data=dat,gamma=1.4)
mod3=gam(sens ~ s(matDim,bs="re",k=4)+s(species,bs="re",k=4) +cv+f,data=dat,gamma=1.4)
mod4=gam(sens ~ s(matDim,bs="re",k=4)+s(species,bs="re",k=4) +cv*f,data=dat,gamma=1.4)


AIC(mod1,mod2,mod3,mod4)

summary(mod3)
### Repeat without outlier
outlier=function(x,k=2) {
  aa<-logical(length(x))
  aa[which(x>mean(x,na.rm=T)+k*sd(x,na.rm=T)|x<mean(x,na.rm=T)-k*sd(x,na.rm=T))]<-T
  aa}

# exclude extreme values of Sv1
mod3a=gam(sens ~ s(matDim,bs="re",k=4)+s(species,bs="re",k=4)+cv+f,data=dat[!outlier(dat$sens),],method="ML")

summary(mod3a)

#exclude extreme values of PCA scores
mod3b=gam(sens ~ s(matDim,bs="re",k=4)+s(species,bs="re",k=4)+cv+f,data=dat[!outlier(dat$PC1)&!outlier(dat$PC2),],method="ML")

summary(mod3b)

##################################

########## STEP 2: Sv1 ~ vital rate + PCA 1 + PCA 2 (example at CV = 0.5 and f = 0.65)


dat1=droplevels(dat[dat$cv=="0.5"&dat$f=="0.65",])
table(dat1$matDim)

dat1=droplevels(dat1[dat1$class!="gNR"&dat1$class!="gP"&dat1$class!="sNR"&dat1$class!="fR"&dat1$class!="fNR",])

dat1$species=factor(dat1$species)

names(which(rowSums(table(dat1$species,dat1$class))==1))

dat1=droplevels(dat1[!dat1$species%in%names(which(rowSums(table(dat1$species,dat1$class))==1)),])

mod0=gam(sens ~ 1,data=dat1,gamma=1.4)
mod1=gam(sens ~ s(matDim,bs="re",k=4)+s(species,bs="re",k=4),data=dat1,gamma=1.4)
mod2= gam(sens ~ s(matDim,bs="re",k=4)+s(species,bs="re",k=4)+te(PC1,k=4) ,data=dat1,gamma=1.4)
mod3= gam(sens ~ s(matDim,bs="re",k=4)+s(species,bs="re",k=4)+te(PC1,k=4) + te(PC2,k=4)  ,data=dat1,gamma=1.4)
mod4= gam(sens ~ s(matDim,bs="re",k=4)+s(species,bs="re",k=4)+te(PC1,PC2,k=4) ,data=dat1,gamma=1.4)

modF= gam(sens ~ s(matDim,bs="re",k=4)+s(species,bs="re",k=4)+te(PC1,PC2,k=4)+ class,data=dat1,gamma=1.4)
modF2=gam(sens ~ s(matDim,bs="re",k=4)+s(species,bs="re",k=4)+te(PC1,PC2,k=4)+ class +te(PC1,by=class,k=4),data=dat1,gamma=1.4)
modF3=gam(sens ~ s(matDim,bs="re",k=4)+s(species,bs="re",k=4)+te(PC1,PC2,k=4)+ class +te(PC1,by=class,k=4)+te(PC2,by=class,k=4),data=dat1,gamma=1.4)

dum=rep(1,nrow(dat1))
modF4=gam(sens ~ te(PC1,PC2,by=class,k=4)+ class+s(matDim,by=dum,bs="re",k=4)+s(species,by=dum,bs="re",k=4),data=dat1,gamma=1.4)

AIC(mod1,mod2,mod3,mod4,modF,modF2,modF3,modF4)

# no outlier in SV1
modF4a=gam(sens ~ te(PC1,PC2,by=class,k=4)+ class+s(matDim,by=dum,bs="re",k=4)+s(species,by=dum,bs="re",k=4),data=dat1[!outlier(dat1$sens),])

summary(modF4a)

# no outlier in PCA scores
modF4b=gam(sens ~ te(PC1,PC2,by=class,k=4)+ class+s(matDim,by=dum,bs="re",k=4)+s(species,by=dum,bs="re",k=4),data=dat1[!outlier(dat1$PC1)&!outlier(dat1$PC2),])

summary(modF4b)

## PLOT

### plot
library(directlabels)
library(fields)
library(lattice)
library(latticeExtra)

dat1$class=factor(dat1$class,levels=c("sP","sPR","sR","gPR","gR","fP","fPR"))

new.data=dat1[dat1$PC2>(-4)&dat1$PC2<6,]

new.data$dum=0

pred=predict(modF4, newdata = new.data,se.fit=T)

new.data$sens=exp(pred$fit)
new.data$up=exp(pred$fit+2*pred$se.fit)
new.data$low=exp(pred$fit-2*pred$se.fit)

data2=NULL

for(i in 1:length(unique(new.data$class))){
  
  sub1=droplevels(new.data[new.data$class==levels(new.data$class)[i],])
  
  fit = Tps(cbind(sub1$PC1,sub1$PC2), sub1$sens)
  pred = predictSurface(fit)
  
  new=melt(pred$z)
  
  new$PC1=rep(pred$x,pred$ny)
  new$PC2=rep(pred$y,each=pred$nx)
  new$class=levels(new.data$class)[i]
  new$species=dat1$species[1]
  new$matDim=dat1$matDim[1]
  new$dum=0
  new$sens=predict(modF4, newdata = new,se.fit=T)$fit
  new$sens[is.na(new$value)]=NA
  
  
  data2=rbind(data2,new)
}
data2$class=factor(data2$class,levels=c("sP","sPR","sR","gPR","gR","fP","fPR"))

breaks <- c(-2.6, -5,-10,-16)

ggplot(data2, aes(PC1, PC2, z =sens))+
  geom_raster(aes(fill = sens),interpolate=F) +
  geom_point(data=dat1, aes(size=exp(sens)),shape=1)+
  facet_wrap( ~ class,ncol=3,scales="fixed")+  
  scale_fill_gradientn(limits=c(-16,-2.6),colours=tim.colors(200),na.value="white",breaks=breaks)+
  guides(color=F,size=F,fill=F)+
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
  theme(strip.text.x = element_blank(),
        strip.background =element_blank())



##################################

########## STEP 3: mean Sv1 ~ PCA 1 + PCA 2 + habitat (example at CV = 0.5 and f = 0.65)

dat=data[data$sens>0,]

levels(dat$class)=c("survival","survival","survival","survival","growth","growth","growth","growth","fec","fec","fec","fec")
dat=aggregate(sens~cv+f+class+species,data=dat,sum)

dat$sens=log(dat$sens)
dat=dat[is.finite(dat$sens),]
dat=left_join(dat,data[,c("species","PC1","PC2","habitat","matDim")],by=c("species"))
dat=dat[-which(duplicated(dat)),]

dat=dat[!is.na(dat$PC1),]

dat$species=factor(dat$species)
dat$matDim=factor(dat$matDim)
levels(dat$matDim)[8:28]="10up"
levels(dat$matDim)[1:2]="4low"

dat1=droplevels(dat[dat$cv=="0.5"&dat$f=="0.65",])

xtabs(~class+habitat,dat1)


mod1=gam(sens ~ te(PC1,PC2,by=class,k=4)+class,data=dat1,gamma=1.4)
mod2=gam(sens ~ te(PC1,PC2,by=class,k=4)+class+habitat,data=dat1,gamma=1.4)
mod3=gam(sens ~ te(PC1,PC2,by=class,k=4)+class+habitat +te(PC1,by=habitat,k=4),data=dat1,gamma=1.4)
mod4=gam(sens ~ te(PC1,PC2,by=class,k=4)+class+habitat +te(PC1,by=habitat,k=4)+te(PC2,by=habitat,k=4),data=dat1,gamma=1.4)
mod5=gam(sens ~ te(PC1,PC2,by=class,k=4)+class+habitat+ te(PC1,PC2,by=habitat,k=4),data=dat1,gamma=1.4)
dum=rep(1,nrow(dat1))
mod6=gam(sens ~ te(PC1,PC2,by=class,k=4)+class+habitat+ te(PC1,PC2,by=habitat,k=4)+s(matDim,by=dum,bs="re",k=4),data=dat1,gamma=1.4)
mod7=gam(sens ~ te(PC1,PC2,by=class,k=4)+class+habitat+ te(PC1,PC2,by=habitat,k=4)+s(matDim,bs="re",k=4)+s(species,bs="re",k=4),data=dat1,gamma=1.4)
mod8=gam(sens ~ te(PC1,PC2,by=class,k=4)+class+habitat+ te(PC1,PC2,by=habitat,k=4)+s(species,bs="re",k=4),data=dat1,gamma=1.4)

AIC(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8)

#### Plot 

new.data=dat1[dat1$PC2>(-4)&dat1$PC2<6,]

new.data$dum=0
pred=predict(mod6, newdata = new.data,se.fit=T)

new.data$sens=exp(pred$fit)
new.data$up=exp(pred$fit+2*pred$se.fit)
new.data$low=exp(pred$fit-2*pred$se.fit)

### plot
library(directlabels)
library(fields)
library(lattice)
library(latticeExtra)

data=NULL

for(i in 1:length(levels(dat1$class))){
  
  for(j in 1:length(levels(dat1$habitat))){
    
    sub1=droplevels(new.data[new.data$class==as.character(levels(dat1$class)[i])&new.data$habitat==as.character(levels(dat1$habitat)[j]),])
    
    fit = Tps(cbind(sub1$PC1,sub1$PC2), sub1$sens)
    pred = predictSurface(fit)
    
    new=melt(pred$z)
    
    new$PC1=rep(pred$x,pred$ny)
    new$PC2=rep(pred$y,each=pred$nx)
    new$class=as.character(levels(dat1$class)[i])
    new$habitat=as.character(levels(dat1$habitat)[j])
    new$matDim=dat1$matDim[150]
    new$dum=0
    new$sens=predict(mod6, newdata = new,se.fit=T)$fit
    new$sens[is.na(new$value)]=NA
    
    
    data=rbind(data,new)
  }
}

dat1$sens=exp(dat1$sens)
dat1$sens[dat1$sens>0.1]=0.1

data$class=factor(data$class,levels=c("survival","growth","fec"))
dat1$class=factor(dat1$class,levels=c("survival","growth","fec"))

data$habitat=factor(data$habitat,levels=c("Alpine & Arctic","Temperate","Tropical & Subtropical","Arid","Aquatic"))
dat1$habitat=factor(dat1$habitat,levels=c("Alpine & Arctic","Temperate","Tropical & Subtropical","Arid","Aquatic"))

breaks <- c(-1, -5, -10,-16)
ggplot(data, aes(PC1, PC2, z =sens))+
  geom_raster(aes(fill = sens),interpolate=F) +
  geom_point(data=dat1, aes(size=exp(sens)),shape=1,alpha=0.3)+
  facet_wrap( class ~ habitat,ncol=5,scales="fixed")+  
  scale_fill_gradientn(limits=c(-16,-1),colours=tim.colors(200),na.value="white",breaks=breaks)+
  guides(color=F,size=F,fill=F)+
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
  theme(strip.text= element_blank(),
        strip.background =element_blank())

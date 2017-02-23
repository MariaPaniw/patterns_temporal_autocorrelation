# Script for Paniw et al. XXXXXX - Appendix S3 - Simulations of temporal autocorrelation in vital rates 

# This script is divided into two parts:
# PART A (lines 33-495): Simulate temporal autocorrelation and calculate sensitivity of the stochastic growht rate to changes in correlation structure
# PART B (lines 399-739): Regress Sv1 (sensitivity to temporal autocorrelation) against PCA scores

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

sp.sub=1:nrow(raw) # the user may wish to work with a smaller subset of species to speed up simulations

lambda.s=array(NA,c(length(sp.sub),ncol(vr),length(v1),length(cv),length(f)))

### ACTUAL SIMULATIONS

# First main loop: frequency of good 

for(z in 1:length(f)){
  
  # Second main loop: coef of var 
  for(i in 1:length(cv)){
    
    # third main loop: autocorrelation
    for(j in 1:length(v1)){
      
      # fourth main loop: species
      
      for(x in 1:length(sp.sub)){
        
        # fifth main loop: vital rates
        
        for(y in 1:ncol(vr)){
          
          #############################################
          ###############################
          # perturb the vr data frame to create bad and good condition mean vital rates
          
          bad.mean=vr[sp.sub[x],]
          
          good.mean=vr[sp.sub[x],]
          
          vr.mean=vr[sp.sub[x],]

          if(y>=10){ # For fecundity, keep regular CV to define variation
            
            if(bad.mean[,y]>0){ #good/bad environment matrices only for non-zero fecundities
              
              # simulate from a gamma distribution
              mu=vr.mean[,y]
              var=(cv[i]*vr.mean[,y])^2
              
              # convert mean and variance to gamma parameters
              shape= mu^2/var
              scale=var/mu
              
              # new means for good and bad are the averages of two truncated gamma distributions
              # both with a mean given by the vital rate at its average across environments 
              # and variance given by the CV 
              
              bad.mean[,y]=mean(rtrunc(10000,spec="gamma",a=0,b=mu,shape=shape,scale=scale)) # upper bound truncated at the mean 
              good.mean[,y]=mean(rtrunc(10000,spec="gamma",a=mu,b=Inf,shape=shape,scale=scale)) # lower bound truncated at the mean
              
              # good and bad vectors:
              
              # within each good and bad environment,
              # 100 random values of a vital rate are produced 
              # from a gamma distribution using the average values for new means of bad/good environments
              # and very low coefficient of variation = 0.01
              var.g=(0.01*good.mean[,y])^2
              var.b=(0.01*bad.mean[,y])^2
              good.vec=rgamma(100,shape=good.mean[,y]^2/var.g,scale=var.g/good.mean[,y])
              bad.vec=rgamma(100,shape=bad.mean[,y]^2/var.b,scale=var.b/bad.mean[,y])
            }else{
              
              # if fecundity vital rate = 0
              good.vec=rep(0,100)
              bad.vec=rep(0,100)
            }
            
            
            
          }else{ 
            
            # For survival/growth, make variation dependent on maximum CV of a binomial variable
            # and use a beta distribution to simulate values of vital rates 
            
            maxCV=sqrt(vr.mean[,y]*(1-vr.mean[,y]))/vr.mean[,y]
            
            if(is.na(maxCV)) maxCV = 0
            
            if(maxCV>0){
              
              # simulate from a beta distribution
              mu=vr.mean[,y]
              var=(cv[i]*maxCV*vr.mean[,y])^2
              
              # convert mean and variance to beta parameters
              shape= mu*((mu*(1-mu))/var-1)
              scale=(1-mu)*((mu*(1-mu))/var-1)
              
              bad.mean[,y]=mean(rtrunc(10000,spec="beta",a=0,b=mu,shape1=shape,shape2=scale)) # for the vital rates that are 0, this will remain 0
              good.mean[,y]=mean(rtrunc(10000,spec="beta",a=mu,b=1,shape1=shape,shape2=scale))  # for the vital rates that are 0, this will remain 0
              
              # good and bad vectors:
              
              # define a new maximum variance
              maxCV.g=sqrt(good.mean[,y]*(1-good.mean[,y]))/good.mean[,y]
              maxCV.b=sqrt(bad.mean[,y]*(1-bad.mean[,y]))/bad.mean[,y]
              
              # new coefficient of variation
              var.g=(0.01*maxCV.g*good.mean[,y])^2
              var.b=(0.01*maxCV.b*bad.mean[,y])^2
              
              shape.g= good.mean[,y]*((good.mean[,y]*(1-good.mean[,y]))/var.g-1)
              scale.g=(1-good.mean[,y])*((good.mean[,y]*(1-good.mean[,y]))/var.g-1)
              
              shape.b= bad.mean[,y]*((bad.mean[,y]*(1-bad.mean[,y]))/var.b-1)
              scale.b=(1-bad.mean[,y])*((bad.mean[,y]*(1-bad.mean[,y]))/var.b-1)
              
              good.vec=rbeta(100,shape1=shape.g,shape2=scale.g)
              bad.vec=rbeta(100,shape1=shape.b,shape2=scale.b)
            }else{
              
              # if survival/growth = 0
              good.vec=rep(0,100)
              bad.vec=rep(0,100)
            }
            
            
          } 
          
          # Constarin survival and growth to not go below 0 or above 1 (just in case)
          
          if(y<10){
            
            good.mean[,y][good.mean[,y]>1]=1
            
          }
          
          bad.mean[,y][bad.mean[,y]<0]=0
          
          # Control for the fact that we have two growth values in first column (growth21 and growth31)
          # This is then technically a multinomial distribution
          
          ## LOOP OVER GOOD AND BAD VECTORS
          
          if(y==5){ # if we are modeling g21, constrain g31 to be 1-g21 if necessary
            
            good.gr31=rep(good.mean[,"g31"],100)
            
            for(gr in 1:100){
              
              if(1-good.vec[gr]-good.mean[,6]<0){
                
                good.gr31[gr] = 1-good.vec[gr]
              }
            }
            
          }else if(y==6){# if we are modeling g31, constrain g21 to be 1-g31 if necessary
            
            good.gr21=rep(good.mean[,"g21"],100)
            
            for(gr in 1:100){
              if((1-good.vec[gr]-good.mean[,5]<0)){
                
                good.gr21[gr] = 1-good.vec[gr]
              }  
              
            }
            
          }
          
          #############################################
          ###############################
          # Create good and bad condition matrices from the perturbed and non-perturbed vital rates
          
          ### Create a new matrix that has 100 rows (simulated values for a give vital rate) and 12 columns(vr)
          
          dat.vec=matrix(rep(as.numeric(vr.mean),each=100),100,ncol(vr),byrow=F)
          
          
          colnames(dat.vec)=c("s1","s2","s3","s4","g21","g31","g32","g43","g34","f13","f23","f33")
          
          dat.vec.g=dat.vec.b=dat.vec
          
          dat.vec.g[,colnames(dat.vec)[y]]=good.vec
          dat.vec.b[,colnames(dat.vec)[y]]=bad.vec
          
          if(y==5){
            
            dat.vec.g[,colnames(dat.vec)[y+1]]=good.gr31
            
          }else if(y==6){
            
            dat.vec.g[,colnames(dat.vec)[y-1]]=good.gr21
          }
          
          ## Array with 1 matrix per row in dataframe
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
          
          lambda.s[x,y,j,i,z]= a1/ts  
          
          
        }
        
      }
    }
  }
}


##############################################################################
#################################
# Calculate sensitivity to temporal autocorrelation

sens=array(0,c(dim(lambda.s)[1] , dim(lambda.s)[2] , dim(lambda.s)[4],dim(lambda.s)[5]))

for(z in 1:dim(lambda.s)[5]){
  
  for(a in 1:dim(lambda.s)[4]){
    
    sens[,,a,z]=(lambda.s[,,1,a,z]-lambda.s[,,2,a,z])/(abs(v1[1]-v1[2]))
    
    # make 0 all the entries where the vital rates were 0
    
    sub=as.matrix(vr[sp.sub,])
    sub[sub==0&!is.na(sub)]=NA
    sub[sub>0&!is.na(sub)]=1
    
    sens[,,a,z]=sens[,,a,z]*sub
    
  }
}


data=adply(sens,c(1,2,3,4))
colnames(data)=c("species","vr","cv","f","sens")

data$species=as.numeric(data$species)
levels(data$vr)=c("s1","s2","s3","s4","g21","g31","g32","g43","g34","f13","f23","f33")
levels(data$cv)=cv
levels(data$f)=f

data$habitat=rep(as.character(raw$habitat[sp.sub]),nrow(data)/length(raw$habitat[sp.sub]))

##################################################################################
################################# 
############# # PART B - Regress Sv1 (sensitivity to temporal autocorrelation) against PCA scores

# Function to remove outliers
outlier=function(x,k=2) {
  aa<-logical(length(x))
  aa[which(x>mean(x,na.rm=T)+k*sd(x,na.rm=T)|x<mean(x,na.rm=T)-k*sd(x,na.rm=T))]<-T
  aa}


##### Running the simulations on all species is time-comsuming
# but some models may not converge if simulations were done on a smaller
# dataset that potentially excluded vital rates

# We therefore provide the data file Sv1_PCA that includes sensitivities 
# for all the vital rates and 428 populations 

data2=read.csv("Sv1_PCA.csv")

data2$cv=as.factor(data2$cv)
data2$f=as.factor(data2$f)
data2$vr=factor(data2$vr,levels=c("s1","s2","s3","s4","g21","g31","g32","g43","g34","f13","f23","f33"))

### CORRELATION DAMPING AND SV1
#### DAMPING RATIOS

damping.r=function(vitalRates){
  
  damping.ratio=rep(NA,nrow(vitalRates))
  
  for(x in 1:nrow(vitalRates)){
  # Create matrix from the mean vital rates per population
  mean.mat=matrix(c(vitalRates[x,"s1"]*(1-vitalRates[x,"g21"]-vitalRates[x,"g31"]-0),vitalRates[x,"s1"]*vitalRates[x,"g21"],vitalRates[x,"s1"]*vitalRates[x,"g31"],0,
                    0,vitalRates[x,"s2"]*(1-0-vitalRates[x,"g32"]-0),vitalRates[x,"s2"]*vitalRates[x,"g32"],0,
                    vitalRates[x,"f13"],vitalRates[x,"f23"],vitalRates[x,"f33"]+vitalRates[x,"s3"]*(1-vitalRates[x,"g43"]),vitalRates[x,"s3"]*vitalRates[x,"g43"],
                    0,0,vitalRates[x,"s4"]*vitalRates[x,"g34"],vitalRates[x,"s4"]*(1-vitalRates[x,"g34"])),
                  4,4,byrow=F)
  
  # calculate damping ratios
  
  values <- eigen(mean.mat) 
  dr <- rle(round(Mod(values$values), 5))$values
  v.dom  <- dr[1]
  
  v.sub <-  dr[2]
  
  damping.ratio[x]=v.sub/v.dom
  
  
}

return(damping.ratio)
}


### Test for correlation at CV=0.5 and f=0.65

sub=data2[data2$cv=="0.5"&data2$f=="0.65",]
sub.mu=aggregate(sens~species,data=sub,mean,na.rm=T)

### add the PCA scores 

sub.mu$PC1=data2$PC1[1:nrow(sub.mu)]
sub.mu$PC2=data2$PC2[1:nrow(sub.mu)]

sub.mu$damping=damping.r(vr)

# do not restrict by PCA

mod=lm(sens~damping,data=sub.mu)
summary(mod) # not significant

# restrict by PCA (only long-lived)

mod=lm(sens~damping,data=sub.mu[sub.mu$PC1<(-1)|sub.mu$PC2>(1),])
summary(mod) # not significant

sub.mu$pred[sub.mu$PC1<(-1)|sub.mu$PC2>(1)]=predict(mod)

ggplot(data=sub.mu[sub.mu$PC1<(-1)|sub.mu$PC2>(1),],aes(damping,sens))+
  geom_point(size=3,col="grey",alpha=0.5)+
  geom_line(aes(damping,pred),color="black",size=1.5)+
  theme_bw()+
  
  xlab("Damping ratio") +
  ylab(expression(bar(S^v1)))+
  theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size=20))+
  theme(axis.title = element_text(size=26))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.9))



####################### REGRESSIONS

dat=data2[!is.na(data2$sens),]

# Drop vital rates with very few samples 
dat=droplevels(dat[dat$vr!="f33"&dat$vr!="g43"&dat$vr!="g31"&dat$vr!="s4"&dat$vr!="g34",])

###################################################
#### A - Sv1 as a function of vital rates, CV, and F

## All candidate models are described in Appendix S3. 
## Here only some models are presented to showcase the model selectino procedure

mod8=lmer(sens~(PC1+I(PC1^2)+PC2+I(PC2^2))*cv+f+vr+(1|species),dat=dat[!outlier(dat$sens),])
mod9=lmer(sens~(PC1+I(PC1^2)+PC2+I(PC2^2))*cv*vr+f+(1|species),dat=dat[!outlier(dat$sens),])

PBmodcomp(mod9,mod8) # model 9 is significantly better

##### Plot predictions of model 
## (here only shown at CV=0.5 and f=0.65; remaining plots can be seen in Appendix S3)

# Predictions for PCA 1 (at average of PCA 2)
new.data=expand.grid(PC1=data2$PC1[1:425],
                     PC2=mean(data2$PC2[1:425]),
                     vr=levels(dat$vr),
                     cv=levels(dat$cv),
                     f=levels(dat$f),
                     sens=0)
mm <- model.matrix(terms(mod9),new.data)


new.data$sens <- mm %*% as.matrix(fixef(mod9))

new.data$PC="PCA 1"
new.data$x=new.data$PC1

# Predictions for PCA 2 (at average of PCA 1)
new.data2=expand.grid(PC1=mean(data2$PC1[1:425]),
                      PC2=data2$PC2[1:425],
                      vr=levels(dat$vr),
                      cv=levels(dat$cv),
                      f=levels(dat$f),
                      sens=0)
mm <- model.matrix(terms(mod9),new.data2)


new.data2$sens <- mm %*% as.matrix(fixef(mod9))
new.data2$PC="PCA 2"
new.data2$x=new.data2$PC2

dataPred=rbind(new.data,new.data2)


dat1=dat
dat1$PC="PCA 1"
dat1$x=dat$PC1
dat2=dat
dat2$PC="PCA 2"
dat2$x=dat$PC2

dataPoints=rbind(dat1,dat2)


### Plot
sub1=droplevels(dataPred[dataPred$cv=="0.5"&dataPred$f=="0.65",])
sub2=droplevels(dataPoints[dataPoints$cv=="0.5"&dataPoints$f=="0.65",])

ggplot(data=sub1,aes(x,sens,col=vr))+
  geom_point(data=sub2,aes(x,sens,col=vr),size=3,alpha=0.5)+
  ylim(-0.05,0.1)+
  facet_wrap( ~ PC, ncol=2,scales="free_x")+
  geom_line(aes(x,sens,col=vr),size=1.5)+
  scale_color_manual(name="",values=c("red4","red","orange","purple","royalblue1",
                                      "greenyellow","green4"),
                     labels=list(bquote(sigma[P]),bquote(sigma[J]),bquote(sigma[R]),bquote(gamma[JP]),bquote(gamma[RJ]),
                                 bquote(phi[PR]),bquote(phi[JR])))+
  theme_bw()+
  theme(legend.title = element_text(size=26, face="bold"),
        legend.text = element_text(size=18),
        legend.key.size = unit(2, "lines"))+
  
  xlab("") +
  ylab(expression(S^v1))+
  ggtitle("f = 0.65; CV = 0.5")+
  theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size=20))+
  theme(axis.title = element_text(size=26))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.9),
        plot.title = element_text(size=22,lineheight=.8, face="bold"))+
  theme(strip.text.x = element_text(size = 18),
        strip.background =element_blank())


###################################################
#### A - Mean Sv1 as a function of  CV,  f, and ecoregion 

data.mu=aggregate(sens~species+cv+f,data=data2,mean,na.rm=T)

head(data.mu)

data.mu$habitat=rep(data2$habitat[1:length(unique(data.mu$species))],nrow(data.mu)/length(unique(data.mu$species)))
data.mu$PC1=rep(data2$PC1[1:length(unique(data.mu$species))],nrow(data.mu)/length(unique(data.mu$species)))
data.mu$PC2=rep(data2$PC2[1:length(unique(data.mu$species))],nrow(data.mu)/length(unique(data.mu$species)))


mod8=lm(sens~(PC1+I(PC1^2)+PC2+I(PC2^2))+f*cv,data=data.mu[!outlier(data.mu$sens),])

mod9=lm(sens~(PC1+I(PC1^2)+PC2+I(PC2^2))*f*cv,data=data.mu[!outlier(data.mu$sens),])

PBmodcomp(mod9,mod8)

##### Plot predictions of model 
## (here only shown at CV=0.5 and f=0.65; remaining plots can be seen in Appendix S3)

# Predictions for PCA 1 (at average of PCA 2)
new.data=expand.grid(PC1=data2$PC1[1:425],
                     PC2=mean(data2$PC2[1:425]),
                     cv=levels(dat$cv),
                     f=levels(dat$f))

new.data$sens=predict(mod9,newdata=new.data,se.fit=T)$fit

new.data$PC="PCA 1"
new.data$x=new.data$PC1

# Predictions for PCA 2 (at average of PCA 1)
new.data2=expand.grid(PC1=mean(data2$PC1[1:425]),
                      PC2=data2$PC2[1:425],
                      cv=levels(dat$cv),
                      f=levels(dat$f))
new.data2$sens=predict(mod9,newdata=new.data2,se.fit=T)$fit
new.data2$PC="PCA 2"
new.data2$x=new.data2$PC2

dataPred=rbind(new.data,new.data2)


dat1=data.mu
dat1$PC="PCA 1"
dat1$x=data.mu$PC1
dat2=data.mu
dat2$PC="PCA 2"
dat2$x=data.mu$PC2

dataPoints=rbind(dat1,dat2)


### Plot
sub1=droplevels(dataPred[dataPred$cv=="0.5"&dataPred$f=="0.65",])
sub2=droplevels(dataPoints[dataPoints$cv=="0.5"&dataPoints$f=="0.65",])

ggplot(data=sub1,aes(x,sens))+
  geom_point(data=sub2,aes(x,sens),col="grey",size=3,alpha=0.5)+
  ylim(-0.05,0.1)+
  facet_wrap( ~ PC, ncol=2,scales="free_x")+
  geom_line(aes(x,sens),col="black",size=1.5)+
  theme_bw()+
  theme(legend.title = element_text(size=26, face="bold"),
        legend.text = element_text(size=18),
        legend.key.size = unit(2, "lines"))+
  guides(color=F) +
  xlab("") +
  ylab(expression(bar(S^v1)))+
  ggtitle("f = 0.65; CV = 0.5")+
  theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size=20))+
  theme(axis.title = element_text(size=26))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.9),
        plot.title = element_text(size=22,lineheight=.8, face="bold"))+
  theme(strip.text.x = element_text(size = 18),
        strip.background =element_blank())


####### INCLUDE HABITAT TYPE


mod2=lm(sens~(PC1+I(PC1^2)+PC2+I(PC2^2))*f*cv+habitat,data=data.mu[!outlier(data.mu$sens),])

mod3=lm(sens~(PC1+I(PC1^2))*habitat+(PC1+I(PC1^2)+PC2+I(PC2^2))*f*cv,data=data.mu[!outlier(data.mu$sens),])

PBmodcomp(mod3,mod2)

### Plot predictions 
new.data=expand.grid(PC1=data2$PC1[1:425],
                     PC2=mean(data2$PC2[1:425]),
                     habitat=levels(data.mu$habitat),
                     cv=levels(data.mu$cv),
                     f=levels(data.mu$f))

new.data$sens=predict(mod3,newdata=new.data,se.fit=T)$fit
new.data$PC="PCA 1"
new.data$x=new.data$PC1

new.data2=expand.grid(PC1=mean(data2$PC1[1:425]),
                      PC2=data2$PC2[1:425],
                      habitat=levels(data.mu$habitat),
                      cv=levels(data.mu$cv),
                      f=levels(data.mu$f))
new.data2$sens=predict(mod3,newdata=new.data2,se.fit=T)$fit
new.data2$PC="PCA 2"
new.data2$x=new.data2$PC2

dataPred=rbind(new.data,new.data2)
dataPred$habitat=factor(dataPred$habitat,levels=c("Alpine & Arctic","Temperate",
                                    "Tropical & Subtropical","Arid","Aquatic"))


dat1=data.mu
dat1$PC="PCA 1"
dat1$x=data.mu$PC1
dat2=data.mu
dat2$PC="PCA 2"
dat2$x=data.mu$PC2

dataPoints=rbind(dat1,dat2)
dataPoints$habitat=factor(dataPoints$habitat,levels=c("Alpine & Arctic","Temperate",
                                                  "Tropical & Subtropical","Arid","Aquatic"))




sub1=droplevels(dataPred[dataPred$cv=="0.5"&dataPred$f=="0.65"&dataPred$PC=="PCA 1",])
sub2=droplevels(dataPoints[dataPoints$cv=="0.5"&dataPoints$f=="0.65"&dataPoints$PC=="PCA 1",])

ggplot(data=sub1,aes(x,sens,col=habitat))+
  geom_point(data=sub2,aes(x,sens,col=habitat),size=3,alpha=0.5)+
  ylim(-0.05,0.1)+
  geom_line(aes(x,sens,col=habitat),size=1.5)+
  scale_color_manual(name="",values=c("yellowgreen","green2","darkgreen","orange",
                                      "blue"))+
  theme_bw()+
  theme(legend.title = element_text(size=26, face="bold"),
        legend.text = element_text(size=18),
        legend.key.size = unit(2, "lines"))+
  
  xlab("") +
  ylab(expression(bar(S^v1)))+
  ggtitle("f = 0.65; CV = 0.5")+
  theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size=20))+
  theme(axis.title = element_text(size=26))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.9),
        plot.title = element_text(size=22,lineheight=.8, face="bold"))+
  theme(strip.text.x = element_text(size = 18),
        strip.background =element_blank())


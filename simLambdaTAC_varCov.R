# Script for Paniw et al. XXXXXX - Appendix S3 - Simulations of temporal autocorrelation in vital rates using observed vital-rate covariation 

# This script is divided into two parts:
# PART A (lines 33-440): Simulate temporal autocorrelation and calculate sensitivity of the stochastic growth rate to changes in correlation structure
# PART B (lines 441-579): Regress Sv1 (sensitivity to temporal autocorrelation) against PCA scores and habitat types

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

load("matsVarCov")

# Create some parameters necessary for the simulation (loop)

# Autocorrelation:

v1=c(-0.3,0,0.3)

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

rep.simul=5

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


# The output array will have 5 dimensions: species x vital rates x v1 x cv x f

sp.sub=1:length(matsVarCov) # the user may wish to work with a smaller subset of species to speed up simulations

lambda.s=vector("list", length(sp.sub))

### ACTUAL SIMULATIONS

# First main loop: ppopulations

for(x in 1:length(sp.sub)){
  
  vr.all=matsVarCov[[sp.sub[x]]]$vr.mu
  
  vr=matsVarCov[[sp.sub[x]]]$vr.mu[matsVarCov[[sp.sub[x]]]$vr.mu>0]
  Umat=matsVarCov[[sp.sub[x]]]$matU
  Fmat=matsVarCov[[sp.sub[x]]]$matF
  
  if(length(which(vr[-grep("f",names(vr))]==1))>0) vr=vr[-which(vr[-grep("f",names(vr))]==1)]
  
  
  sigma=matsVarCov[[sp.sub[x]]]$corr
  if(x==53|x==55) colnames(sigma)=names(vr)
  
  
  if(any(vr>0)){
    if(sp.sub[x]==34|sp.sub[x]==38|sp.sub[x]==55){
      
      sigma[lower.tri(sigma)]=t(sigma)[lower.tri(sigma)]
      
      if(sp.sub[x]==34){sigma=sigma[1:5,1:5]}
      sigma.sub=pos.def(sigma[colnames(sigma)%in%names(vr[vr>0]),colnames(sigma)%in%names(vr[vr>0])])
      
    }else{
      sigma.sub=sigma[colnames(sigma)%in%names(vr[vr>0]),colnames(sigma)%in%names(vr[vr>0])]
    }
    
    colnames(sigma.sub)=rownames(sigma.sub)=colnames(sigma[colnames(sigma)%in%names(vr[vr>0]),colnames(sigma)%in%names(vr[vr>0])])
    vr.cor=vr[names(vr)%in%rownames(sigma.sub)]
    
    min.f=ifelse(is.finite(min(grep("f",names(vr.cor)))),min(grep("f",names(vr.cor))),length(vr.cor)+20)
    # real variances:
    
    real.var=matsVarCov[[sp.sub[x]]]$var[names(matsVarCov[[sp.sub[x]]]$var)%in%names(vr.cor)]
    
    if(!any(sigma.sub!=0)) sigma.sub<-NULL
    
  }else{sigma.sub=NULL}
  
  lambda.s.temp=array(NA,c(length(v1),length(f),rep.simul))
  var.temp=array(NA,c(length(v1),length(f),rep.simul))
  
  # Second main loop: Replicate simulations 
  
  for(sm in 1:rep.simul){
    number.sim=seq(100,149)[sm]
    
    for(z in 1:length(f)){
      
      
      for(j in 1:length(v1)){
        
        ## SIMULATE WITHIN-STATE VARIANCE (USING COPULAS)
        corr.values=matrix(rep(vr,each=1000),1000,length(vr),byrow=F)
        colnames(corr.values)=names(vr)
        
        if(!is.null(nrow(sigma.sub))){
          
          set.seed(number.sim)
          z.m <- mvrnorm(1000,mu=rep(0,nrow(sigma.sub)),Sigma=sigma.sub)
          #Apply the Normal CDF function to Z to obtain data that is uniform 
          # on the interval [0,1], but still correlated.
          
          u.m <- pnorm(z.m)
          
          for(mm in 1:length(names(vr.cor))){
            
            # create marginal beta distribution for survival/transitions
            if(mm<min.f){
              
              mu=min(0.99999,vr[names(vr.cor)[mm]])
              var.real=real.var[names(vr.cor)[mm]]
              CV.real=sqrt(var.real)/mu
              
              maxCV=sqrt(mu*(1-mu))/mu
              
              maxCV <- ifelse(CV.real>maxCV,0.99*maxCV,CV.real)
              
              var=(maxCV*mu)^2
              
              
              a= mu*((mu*(1-mu)/var)-1)
              b=(1-mu)*((mu*(1-mu))/var-1)
              
              set.seed(number.sim)
              # add marginal beta to Z
              corr.values[,names(vr.cor)[mm]]=qbeta(u.m[,names(vr.cor)[mm]],shape1=a,
                                                    shape2=b)
              
              # create marginal gamma distribution for reproduction 
            }else if(mm>=min.f){
              
              mu=vr[names(vr.cor)[mm]]
              
              var=real.var[names(vr.cor)[mm]]
              a= mu^2/var
              b=var/mu
              
              set.seed(number.sim)
              
              corr.values[,names(vr.cor)[mm]]=qgamma(u.m[,names(vr.cor)[mm]],shape=a,
                                                     scale=b)
              
            }
          }
          
          
        }else{
          
          corr.values=corr.values
          
          
        }
        
        
        #### CREATE AN ARRAY OF 1000 matrices
        ##### PUT ALL VITAL RATES TOGETHER
        
        mat.array=array(0,c(dim(Umat)[1],dim(Umat)[1],1000))
        mean.mat=Umat+Fmat
        
        good.bad=rep(NA,1000)
        lambda.good=NULL
        lambda.bad=NULL
        lambda.mu=as.numeric(max(abs(eigen(mean.mat)$values)))
        
        for(dist in 1:1000){
          vr2=c(corr.values[dist,],vr.all[!names(vr.all)%in%colnames(corr.values)])
          vr2=vr2[names(vr.all)]
          
          F.mat.new=Fmat
          
          surv=colSums(Umat)
          U.mat.gr=Umat
          for(xx in 1:ncol(U.mat.gr)){
            U.mat.gr[,xx]=Umat[,xx]/surv[xx]
            
          }
          
          U.mat.gr[!is.finite(U.mat.gr)]=0
          
          for(y in 1:ncol(U.mat.gr)){
            
            if(y<10){all=which((substring(names(vr2), 1, 1)=="g"|substring(names(vr2), 1, 1)=="r")&substring(names(vr2), 4, nchar(names(vr2)))==paste(0,y,sep=""))}else if(y>=10){
              
              all=which((substring(names(vr2), 1, 1)=="g"|substring(names(vr2), 1, 1)=="r")&substring(names(vr2), 4, nchar(names(vr2)))==paste(y,sep=""))
            }
            
            vr.sub=vr2[all]
            
            if(diag(mean.mat)[y]==0){
              
              stasis=0
              
              if((1-sum(vr.sub))<0){
                
                vr.sub[vr.sub>0&(vr.sub-abs(1-sum(vr.sub)))>0]=vr.sub[vr.sub>0&(vr.sub-abs(1-sum(vr.sub)))>0]-abs(1-sum(vr.sub))/length(vr.sub[vr.sub>0&(vr.sub-abs(1-sum(vr.sub)))>0])
              }else if((1-sum(vr.sub))>0){
                
                vr.sub[vr.sub>0]=vr.sub[vr.sub>0]+abs(1-sum(vr.sub))/length(vr.sub[vr.sub>0])
              }
              
            }else{
              
              stasis=1-sum(vr.sub)
              if(stasis<0){
                vr.sub[vr.sub>0&(vr.sub-abs(stasis))>0]=vr.sub[vr.sub>0&(vr.sub-abs(stasis))>0]-abs(stasis)/length(vr.sub[vr.sub>0&(vr.sub-abs(stasis))>0])
                stasis=1-sum(vr.sub)
              }
              
            }
            
            diag(U.mat.gr)[y] <- stasis
            U.mat.gr[,y][lower.tri(U.mat.gr)[,y]]=vr.sub[grep("g",names(vr.sub))]
            U.mat.gr[,y][upper.tri(U.mat.gr)[,y]]=vr.sub[grep("r",names(vr.sub))]
            
            U.mat.gr[,y]=U.mat.gr[,y]*vr2[grep("s",names(vr2))][y]
          }
          
          for(yy in 1:length(vr2[grep("f",names(vr2))])){
            
            F.mat.new[as.numeric(substring(names(vr2[grep("f",names(vr2))])[yy], 2, 3)),as.numeric(substring(names(vr2[grep("f",names(vr2))])[yy], 4, nchar(names(vr2[grep("f",names(vr2))][yy]))))]= vr2[grep("f",names(vr2))][yy]
          }
          
          mat.array[,,dist]=U.mat.gr+F.mat.new
          
          lambda.temp=as.numeric(max(abs(eigen(mat.array[,,dist])$values)))
          
          if(sp.sub[x]!=102){
            if(round(lambda.temp,3)>round(lambda.mu,3)){
              good.bad[dist]<-1
              lambda.good=c(lambda.good,lambda.temp)
            }else if(round(lambda.temp,3)<round(lambda.mu,3)){
              
              good.bad[dist]<-2
              lambda.bad=c(lambda.bad,lambda.temp)
            }else{good.bad[dist]<-0}
          }else{
            
            if(lambda.temp>1){
              good.bad[dist]<-1
              lambda.good=c(lambda.good,lambda.temp)
            }else if(lambda.temp<1){
              
              good.bad[dist]<-2
              lambda.bad=c(lambda.bad,lambda.temp)
            }else{good.bad[dist]<-0}
          }
          
          
          
        }
        
        
        # Create environmental transition matrix P (working with f and v1)
        
        q=f[z]*(1-v1[j])
        
        p=v1[j]+q
        
        
        P=matrix(c(p,1-p,
                   q,1-q),nrow=2,ncol=2,byrow=F)
        colnames(P) <- c("1","2")
        row.names(P) <- c("1","2")
        
        vr.sim.result=rep(0,50)
        for(vr.sim in 1:50){
          simul.n = numeric(trun+1)
          
          # get sequence of environments
          # start with god environment 
          env_at_t=1
          simul.n[1]=1
          
          for(sim_t in 2:trun)
            simul.n[sim_t] = env_at_t =sample(ncol(P),1,pr=P[,env_at_t])
          
          
          # Calculate stochastic lambda
          ########################################################
          nstage <- dim(mat.array)[1] 
          states <- simul.n
          growth <- array(0,ts)   
          
          # Initial population vector
          vec1=rep(1,nstage)
          vec1 <- vec1/sum(vec1)
          vec1 <- t(vec1) 
          
          # ITERATION TO CALCULATE LAMBDA FOR EACH TIME STEP
          for (a  in 1:trun){
            i3 <- states[a]
            if(i3==1){i2<-sample(which(good.bad==1),1)}else{i2<-sample(which(good.bad==2),1)}
            mat1 <-  mat.array[,,i2]
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
        a.hat=sum(vr.sim.result)/(3*ts)
        var.hat=sum((vr.sim.result-a.hat*ts)^2)/(3*ts)
        lambda.s.temp[j,z,sm]= a.hat
        
        var.temp[j,z,sm]= var.hat
        
        
      }
      
    }
  }

  lambda.s[[x]]$lambda.s=lambda.s.temp
  lambda.s[[x]]$lambda.s.var=var.temp
  lambda.s[[x]]$names.vr=names(vr.cor)
}


##############################################################################
################################# Put results into data frame
sens=NULL

for(x in 1:length(matsVarCov)){
  
  lambda=lambda.s[[x]]$lambda.s
  
  # save sensitivities (Sv1) per f and simulation run
  sens.temp=array(0,c(dim(lambda)[2],dim(lambda)[3]))
  
  for(xx in 1:dim(lambda)[3]){
    # sensitivities are extracted from v1 = 0.3 - v1 = 0
    sens.temp[1,xx]=(lambda[3,1,xx]-lambda[2,1,xx])/(abs(v1[3]-v1[2])) # f = 0.3
    sens.temp[2,xx]=(lambda[3,2,xx]-lambda[2,2,xx])/(abs(v1[3]-v1[2])) # f = 0.65
  }
  
  sens.temp2=apply(sens.temp,1,mean) # average Sv1 across simualtion runs
  
  data=adply(sens.temp2,c(1))
  colnames(data)=c("f","sens")
  data$var.sens=apply(sens.temp,1,var)
  
  data$species=matsVarCov[[x]]$species
  data$matDim=dim(matsVarCov[[x]]$matU)[1]
  levels(data$f)=f
  
  sens=rbind(sens,data)
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

# We therefore provide the data file Sv1_PCA_varCov that includes sensitivities 
# for all the vital rates and 455 populations 

data=read.csv("Sv1_PCA_varCov.csv")

data$sens=abs(data$sens)

data$habitat=as.factor(data$habitat)

data$species=as.factor(data$species)
data$matDim=as.factor(data$matDim)
levels(data$matDim)[7:10]="9up"
levels(data$matDim)[1:2]="4low"


########## STEP 1: Sv1 ~ f
dat=data[data$sens>0,]
dat$sens=log(dat$sens)

mod1=gam(sens ~ 1,data=dat,gamma=1.4)
mod1a=gam(sens ~ s(matDim,bs="re"),data=dat,gamma=1.4)
mod2= gam(sens ~ f +s(matDim,bs="re"),data=dat,gamma=1.4)

AIC(mod1,mod1a,mod2)

##################################

########## STEP 2: Sv1 ~ habitat + PCA 1 + PCA 2 (example at f = 0.65)

dat1=droplevels(dat[dat$f=="0.65",])
levels(dat1$habitat)[c(1,2,3,5)]="other"

dat1$habitat=factor(dat1$habitat,levels=c("Temperate","other"))

xtabs(~habitat,dat1)

mod1=gam(sens ~ 1,data=dat1,gamma=1.4)
mod1a=gam(sens ~ s(matDim,bs="re",k=4),data=dat1,gamma=1.4) # MPM dimension not significant
mod2= gam(sens ~ te(PC1,k=4),data=dat1,gamma=1.4)
mod2a= gam(sens ~ te(PC1,k=4) + te(PC2,k=4) ,data=dat1,gamma=1.4)
mod3= gam(sens ~  te(PC1,PC2,k=4),data=dat1,gamma=1.4)

AIC(mod1,mod2,mod2a,mod3)

summary(mod3)

mod4= gam(sens ~ te(PC1,PC2,k=4) +habitat ,data=dat1,gamma=1.4)
mod5= gam(sens ~ te(PC1,PC2,k=4) + te(PC1,by=habitat,k=4) +habitat,data=dat1,gamma=1.4)
mod6= gam(sens ~ te(PC1,PC2,k=4) + te(PC2,by=habitat,k=4) +habitat,data=dat1,gamma=1.4)
mod7= gam(sens ~ te(PC1,PC2,k=4) + te(PC1,by=habitat,k=4)+ te(PC2,by=habitat,k=4) +habitat ,data=dat1,gamma=1.4)
mod8= gam(sens ~ te(PC1,PC2,by=habitat,k=4)  +habitat ,data=dat1,gamma=1.4)
AIC(mod3,mod4,mod5,mod6,mod7,mod8)

## plot 

new.data=dat1

pred=predict(mod8, newdata = new.data,se.fit=T)

new.data$sens=exp(pred$fit)
new.data$up=exp(pred$fit+2*pred$se.fit)
new.data$low=exp(pred$fit-2*pred$se.fit)


### plot
library(directlabels)
library(fields)
library(lattice)
library(latticeExtra)

data=NULL

################ Temperate
sub1=droplevels(dat1[dat1$habitat=="Temperate",])

fit = Tps(cbind(sub1$PC1,sub1$PC2), sub1$sens)
pred = predictSurface(fit)

new=melt(pred$z)

new$PC1=rep(pred$x,pred$ny)
new$PC2=rep(pred$y,each=pred$nx)
new$habitat="Temperate"
new$sens=predict(mod8, newdata = new,se.fit=T)$fit
new$sens[is.na(new$value)]=NA

data=rbind(data,new)

################ Other
sub1=droplevels(dat1[dat1$habitat=="other",])

fit = Tps(cbind(sub1$PC1,sub1$PC2), sub1$sens)
pred = predictSurface(fit)

new=melt(pred$z)

new$PC1=rep(pred$x,pred$ny)
new$PC2=rep(pred$y,each=pred$nx)
new$habitat="other"
new$sens=predict(mod8, newdata = new,se.fit=T)$fit
new$sens[is.na(new$value)]=NA

data=rbind(data,new)
data$habitat=factor(data$habitat,levels=c("Temperate","other"))

dat1$sens=exp(dat1$sens)
dat1$sens[dat1$sens>0.05]=0.05

breaks=c(-3.8,-8,-11.2)
ggplot(data, aes(PC1, PC2, z =sens))+
  geom_raster(aes(fill = sens),interpolate=F) +
  geom_point(data=dat1,aes(size=sens),shape=1)+
  facet_wrap( ~ habitat,ncol=2,scales="fixed", labeller = label_parsed)+ 
  scale_fill_gradientn(limits=c(-11.2,-3.8),colours=tim.colors(128),na.value="white",breaks=breaks)+
  guides(color=F,size=F,fill=F)+
  # ggtitle("f = 0.65; CV = 0.5")+
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

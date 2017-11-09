# Script for Paniw et al. XXXXXX - Appendix S4 - Comparison of simulation and empirical studies 

rm(list=ls(all=TRUE))
# load species with long-term data
library(dplyr)
setwd("/Users/mariapaniw/Dropbox/TempAutoProject/SuppMat")

data=read.csv("long_term_studies.csv")

data=data[!is.na(data$vr),]

sp.names=data$species[!duplicated(data$species)]

data.dev=NULL
for(i in 1:length(sp.names)){
  
  sub=droplevels(data[data$species==sp.names[i],])
  
  vr.all=unique(sub$vr.name)
  
  for(y in 1:length(vr.all)){
    
    response=sub$vr[sub$vr.name==vr.all[y]]
    
    
    if(var(response,na.rm=T)>0.0001){
      
      n=length(response)
      predictor=scale(c(NA,response[1:(n-1)]))
      response=response[!is.na(predictor)]
      predictor=predictor[!is.na(predictor)]
      
      # get deviance for survival/transitions 
      if(length(grep("f",vr.all[y]))<1){
        
        if(any(response>1)){response[response>1]=1}
        tryCatch({
          
          mod=glm(response~predictor,family=quasibinomial)

          data.dev=rbind(data.dev,data.frame(dev=(mod$null.deviance-mod$deviance)/mod$deviance,
                                               species=sp.names[i],
                                               vr.name=vr.all[y],
                                             rel.elast=data$rel.elast[data$species==sp.names[i]&data$vr.name==as.character(vr.all[y])][1]))
          
          
        },error=function(e) NA)
        
        
      }else{
        response[response==0]<-0.000001
        
        tryCatch({
          
          # get deviance for reproduction
          mod=glm(response~predictor,family=Gamma)
          mod0=glm(response~1,family=Gamma)
          
          data.def=rbind(data.def,data.frame(dev=(deviance(mod0)-deviance(mod))/deviance(mod0),
                                               species=sp.names[i],
                                               vr.name=vr.all[y],
                                             rel.elast=data$rel.elast[data$species==sp.names[i]&data$vr.name==as.character(vr.all[y])][1]))
          
        },error=function(e) NA)
        
        
      }
    }else{
      
      data.dev=rbind(data.dev,data.frame(dev=0,
                                           species=sp.names[i],
                                           vr.name=vr.all[y],
                                         rel.elast=data$rel.elast[data$species==sp.names[i]&data$vr.name==as.character(vr.all[y])][1]))
      
    }
    
    
  }
  
  
}

data.dev$dev[!is.na(data.dev$dev)&(data.dev$dev>1|data.dev$dev<0)]=NA

# Take average deviance 
dev.mu=aggregate(dev~species,mean,data=data.dev[data.dev$rel.elast>0&!is.na(data.dev$rel.elast),],na.rm=T)

data2=read.csv("Sv1_PCA_varCov.csv")
data2=data2[data2$f==0.65,]
data2$sens=abs(data2$sens)
dev.mu$sens=left_join(dev.mu,data2,by="species")$sens


# Plot 

library(ggplot2)

ggplot(data=dev.mu,aes(log(dev),log(sens)))+
  geom_point(size=3,col="grey20")+
  stat_smooth(method = "lm",fill="grey70",color="black")+
  theme_bw()+
  xlim(-9,-2)+
  xlab(expression(log(bar(D^2)))) +
  ylab(expression(log(S^v1)))+
  theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size=20,colour = "black"))+
  theme(axis.title = element_text(size=26))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.9),
        plot.title = element_text(size=22,lineheight=.8, face="bold"))+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))


####### Part B of plot 4: box plots

data=read.csv("phyloSpecies.csv")

data=data[!duplicated(data$SpeciesAuthor),]
colnames(data)[1]="species"
# get sensitivities:

data2=read.csv("Sv1_PCA.csv")

# subset to cv=0.5 and f=0.65

sens=data2[data2$cv==0.5&data2$f==0.65,]

sens=sens[!is.na(sens$sens),]
sens$sens=abs(sens$sens)
sens=sens[sens$sens>0,]

# get mean sensitivities across vital rates

data2.mu=aggregate(sens~species,sum,data=data2)
data2.mu$sens=log(data2.mu$sens)
data$sens=left_join(data,data2.mu,by=c("species"))$sens

# calculate percentiles of Sv1

data$sens.q=ecdf(data$sens)(data$sens)

# Box plot (subset to studies in which environmental variation was modeled)

sub=data[!is.na(data$envStochModeled)&(data$envStochModeled=="yes"|data$envStochModeled=="mega"),]

ggplot(data=sub[!is.na(sub$TACmodeled),],aes(TACmodeled,sens))+
  geom_boxplot()+
  geom_jitter(shape=1, size=1.5,position=position_jitter(0.15),col="grey")+
  theme_bw()+
  xlab(expression(paste(v[1]," assessed empirically")))+
  
  ylab(expression(paste(log(S^v1))))+
  theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size=20,colour="black"))+
  theme(axis.title = element_text(size=24))+
  theme(axis.title.x = element_text(vjust=2),
        axis.title.y = element_text(vjust=2),
        plot.title = element_text(size=22,lineheight=.8, face="bold"))+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))

# check for significance

t.test(sub$sens[sub$TACmodeled=="yes"&!is.na(sub$TACmodeled)],sub$sens[sub$TACmodeled=="no"&!is.na(sub$TACmodeled)],alternative = "greater")

# Some subsets - Stochastic population dynamics modeled:

### Temporal autocorrelation is important:

sub2=sub[sub$TACmodeled=="yes"&!is.na(sub$TACmodeled),]

nrow(sub2)/nrow(sub) # temporal autocorrelation was modeld in 20 % of studies where stochastic population projections were used


# how many with high Sv1?

nrow(sub2[sub2$sens.q>0.25,])/nrow(sub2) # Sv1 were > 25th percentile for 82 % of studies with stochastic population dynamics considered

nrow(sub2[sub2$sens.q>0.5,])/nrow(sub2) # and 64 were above the 50th percentile


### Temporal autocorrelation is not considered important:

sub3=sub[sub$TACmodeled=="no"&!is.na(sub$TACmodeled),]

# how many with high Sv1?

nrow(sub3[sub3$sens.q>0.25,])/nrow(sub3) # Sv1 were > 25th percentile for 73 % of studies with stochastic population dynamics considered

nrow(sub3[sub3$sens.q>0.5,])/nrow(sub3) # and 47 were above the 50th percentile



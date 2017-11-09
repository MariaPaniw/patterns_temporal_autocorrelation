# Script for Paniw et al. XXXXXX - Appendix S1
#This script peruses through COMPADRE and COMADRE and outputs vital rates and vital rate groups from each average matrix population model of the 455 study populations
#Author: Rob Salguero-Gomez & Maria Paniw
#Created: 19 Aug 2011


#Clean memory
rm(list=ls(all=TRUE))

library(stringr)
#Functions to obtain vital-rate groups from matrices if matrix dimensions are > 4


#Function to collapse a matrix to a given number of sizes given a single range of contiguous classes
#Code developed by Rob Salguero-Gomez (r.salguero@uq.edu.au) based on #Salguero-Gomez & Plotkin Am.
#Nat 2010. matA is a square matrix, collapse is a vector indicating which stages must be collapsed
# and retained as are. E.g. collapse = c("1","2","3-4","5") will collapse stages 3 and 4, whereas
# collapse = c("1","2","3-5") will collapse stages 3 through 5, inclusive. The argument "collapse"
# must contain all stages, and these must be specified as character ("")


collapseMatrix <- function(matU,matF,collapse){

  matA=matU+matF
  collapseUnique = collapse
	originalDim = dim(matA)[1]
	collapseDim = length(collapseUnique)
	P = matrix(0,nrow=collapseDim , ncol= originalDim)

	splitCollapseUnique=strsplit(collapse,"-")
	for (i in 1:collapseDim){
	  columns=as.numeric(splitCollapseUnique[[i]])
	  # P[i,columns]=1
	  
	  if(!is.na(columns[1])){
	    P[i,(columns[1]:columns[length(columns)])]=1
	  }
	  
	}
	
	
	Q=t(P)
	w=Re(eigen(matA)$vectors[,which(Re(eigen(matA)$values)==max(Re(eigen(matA)$values)))])
	w=w/sum(w)
	
	columns=which(colSums(Q)>1)
	  for (j in columns){
	    rows=which(Q[,j]==1)
	    for (i in rows){
	      Q[i,j]=w[i]/sum(w[rows])
	    }
	  }
	collapseA = P %*% matA %*% Q
	collapseU = P %*% matU %*% Q
	collapseF = P %*% matF %*% Q
	collapse=NULL
  collapseList=list("matA"=collapseA,
                "matU"=collapseU,
	              "matF"=collapseF)
	
	return(collapseList)
}
	
#Function to re-arrange stages to propagule (prop), juvenile (preRep), reproductive (rep), and non-reproductive (nonRep)

rearrangeMatrix <- function(matU,matF){
  reArrange=NULL
  matDim=dim(matF)[1]
  Rep=which(colSums(matF)>0)
  #These are stages that are inter-reproductive but are truly non-reproductive
  nonRepInterRep=Rep[1]-1+as.numeric(c(which(colSums(as.matrix(matF[,Rep[1]:Rep[length(Rep)]]))==0)))
  if(length(nonRepInterRep)>0){
    allElseStages=1:matDim
    allElseStages=allElseStages[-which(allElseStages%in%nonRepInterRep)]
    reArrangeStages=c(allElseStages,nonRepInterRep)
    reArrangeMatU=matU[reArrangeStages,reArrangeStages]
    reArrangeMatF=matF[reArrangeStages,reArrangeStages]
    reArrange$matU=reArrangeMatU
    reArrange$matF=reArrangeMatF
  }
  if(length(nonRepInterRep)==0){
    reArrange[[1]]$matU=matU
    reArrange[[1]]$matF=matF
  }
  return(reArrange)
}

#Function to determine pre-reproductive, reproductive and non-reproductive stages from a matrix model
#which pre-formats them to use them as the "collapse" argument of the function above to collapse matrices

reprodStages <- function(matF,matrixStages=NULL){
  propStage=NULL
  if ("prop"%in%matrixStages) {propStage=which(matrixStages=="prop")} else {propStage=NA}
  
  matDim=dim(matF)[1]
  Rep=which(colSums(matF)>0)
  if(min(Rep)==1){preRep=NA}else if(!is.na(propStage[1])&(min(Rep)-max(propStage)==1)){preRep=NA}else{preRep=min(which(matrixStages=="active")):(min(Rep)-1)}
  if(max(Rep)==matDim){nonRep=NA}else{nonRep=(max(Rep)+1):matDim}
  
  if(length(propStage)>1){propStages=paste(propStage[1],"-",propStage[length(propStage)],sep="")}else{propStages=as.character(propStage)}
  if(length(preRep)>1){preRepStages=paste(preRep[1],"-",preRep[length(preRep)],sep="")}else{preRepStages=as.character(preRep)}
  if(length(Rep)>1){repStages=paste(Rep[1],"-",Rep[length(Rep)],sep="")}else{repStages=as.character(Rep)}
  if(length(nonRep)>1){nonRepStages=paste(nonRep[1],"-",nonRep[length(nonRep)],sep="")}else{nonRepStages=as.character(nonRep)}
  
  stages=c(propStages,preRepStages,repStages,nonRepStages)
  
  return(stages)
}


# Set the working directory, then load the COMPADRE data:
dir <- setwd("/Users/mariapaniw/Dropbox/TempAutoProject/SuppMat") # CHANGE THIS TO YOUR DIRECTORY
load(paste(dir,"/COMPADRE_v.4.0.0.RData",sep=""))
load(paste(dir, "/COMADRE_v.2.0.0.RData", sep=""))

# get IDs of species:
sp.data=read.csv("phyloSpecies.csv")

ID=as.character(unique(sp.data$SpeciesAuthor))

#Subsetting all matrices available in COMPADRE (plants/algae)

indexCOMPADRE=which(duplicated(compadre$metadata$SpeciesAuthor)==FALSE&
                      compadre$metadata$SpeciesAuthor%in%ID)   # mean matrices (across sites and times)

compadre2=NULL
compadre2$metadata=compadre$metadata[indexCOMPADRE,]
compadre2$mat=compadre$mat[indexCOMPADRE]
compadre2$matrixClass=compadre$matrixClass[indexCOMPADRE]

#Subsetting all matrices available in COMADRE (animals) - same criteria as for COMPADRE

indexCOMADRE=which(duplicated(comadre$metadata$SpeciesAuthor)==FALSE&
                     comadre$metadata$SpeciesAuthor%in%ID)

comadre2=NULL
comadre2$metadata=comadre$metadata[indexCOMADRE,]
comadre2$mat=comadre$mat[indexCOMADRE]
comadre2$matrixClass=comadre$matrixClass[indexCOMADRE]

#Loop through the subset of COM(P)ADRE matrices and metadata information 
# to obtain relevant vital rates and information on study details and habitat types:

dims=dim(comadre2$metadata)[1]+dim(compadre2$metadata)[1]

mats=vector("list", dims)

count=0
for (loop in 1:2){
  if (loop==1){d=compadre2}else{d=comadre2}
  
  for (i in 1:dim(d$metadata)[1]){
    tryCatch({
      
      count=count+1
      
      ## periodicity
      
      per=sp.data$AnnualPeriodicity[sp.data$SpeciesAuthor==d$metadata$SpeciesAuthor[i]][1]
      
      ### survival
      surv=colSums(d$mat[[i]]$matU)^per
      if(surv[length(surv)]>0.99999) surv[length(surv)]<-0.995
      
      names(surv)=paste("s",1:length(surv),sep="")
      
      U.mat=d$mat[[i]]$matU
      U.mat.g=U.mat
      for(xx in 1:ncol(U.mat.g)){
        U.mat.g[,xx]=U.mat.g[,xx]/colSums(U.mat)[xx]
        U.mat[,xx]=U.mat.g[,xx]*surv[xx]
      }
      U.mat[!is.finite(U.mat)]=0
      U.mat.g[!is.finite(U.mat.g)]=0
      # progression 
      
      gr=U.mat.g[lower.tri(U.mat.g)]
      
      names=NULL
      for(x in 1:length(surv[-1])){
        
        x1=str_pad(x, 2, pad = "0")
        x2=str_pad((x+1):length(surv), 2, pad = "0")
        
        temp=paste("g",paste(x2,x1,sep=""),sep="")
        names=c(names,temp) 
        
        
      }
      
      names(gr)=names
      
      # retrogression 
      
      ret=U.mat.g[upper.tri(U.mat.g)]
      
      names=NULL
      for(x in 2:length(surv)){
        
        x1=str_pad(x, 2, pad = "0")
        x2=str_pad(1:(x-1), 2, pad = "0")
        
        temp=paste("r",paste(x2,x1,sep=""),sep="")
        names=c(names,temp)
      }
      
      names(ret)=names
      
      # reproduction
      placeholder=matrix(1:length(as.numeric(d$mat[[i]]$matF)),dim(d$mat[[i]]$matF)[1],dim(d$mat[[i]]$matF)[1])
      colnames(placeholder)=rownames(placeholder)=1:dim(d$mat[[i]]$matF)[1]
      fec.names=which(d$mat[[i]]$matF>0)
      names=expand.grid(rownames(placeholder),colnames(placeholder))[placeholder%in%fec.names,]
      names=interaction(str_pad(as.numeric(names$Var1),2,pad="0"),str_pad(as.numeric(names$Var2),2,pad="0"),sep="")
      fec=paste("f",names ,sep="")
      fec.value=d$mat[[i]]$matF[d$mat[[i]]$matF>0]*per #account for periodicity
      names(fec.value)=fec
      
      mats[[count]]$vr=c(surv,gr,ret,fec.value)
      
      
      mats[[count]]$matF=d$mat[[i]]$matF*per
      mats[[count]]$matU=U.mat
      mats[[count]]$species=as.character(d$metadata$SpeciesAuthor[i])[1]
      mats[[count]]$speciesAccepted=as.character(d$metadata$SpeciesAccepted[i])[1]
      
      # Adjudt matrix class to represent P, PR, R, NR
      
      if(dim(d$mat[[i]]$matU)[1]>3){
        
        if(length(rearrangeMatrix(d$mat[[i]]$matU,d$mat[[i]]$matF))==1){
          matU.rea=rearrangeMatrix(d$mat[[i]]$matU,d$mat[[i]]$matF)[[1]]$matU
          matF.rea=rearrangeMatrix(matU.rea,d$mat[[i]]$matF)[[1]]$matF
        }else{
          
          matU.rea=rearrangeMatrix(d$mat[[i]]$matU,d$mat[[i]]$matF)$matU
          matF.rea=rearrangeMatrix(matU.rea,d$mat[[i]]$matF)$matF
        }
     
        matrixStages=d$matrixClass[[i]][,1]
        collapse=cbind(reprodStages(matF.rea,matrixStages),c("P","PR","R","NR"))
        
        collapse=matrix(collapse[!is.na(collapse[,1]),],ncol=2)
        collapse[,1]=gsub('-', ':', collapse[,1])
        
        for(cc in 1:nrow(collapse)){
          if(length(grep(":",collapse[cc,1]))>0){
            matrixStages[as.numeric(strsplit(collapse[cc,1], "\\:+")[[1]][1]):as.numeric(strsplit(collapse[cc,1], "\\:+")[[1]][2]),1]=collapse[cc,2]
            
          }else{
            
            matrixStages[as.numeric(collapse[cc,1]),1]=collapse[cc,2]
            
          }
        }
        rm(matU.rea,matF.rea,collapse)
      }else{
        
        matrixStages=d$matrixClass[[i]][,1]
        mtst=matrixStages
        matrixStages[1:(min(which(colSums(mats[[count]]$matF)>0))-1),1]="PR"
        matrixStages[mtst=="prop",1]="P"
        matrixStages[which(colSums(mats[[count]]$matF)>0),1]="R"
        nr=which(colSums(mats[[count]]$matF)==0)
        nr=nr[nr>min(which(colSums(mats[[count]]$matF)>0))]
        if(length(nr)>0) matrixStages[nr,1]="NR"
      }
      
      mats[[count]]$class=matrixStages%>% pull(1)
      
      rm(matU.rea,matF.rea,collapse)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n");print(i)})
  }
}

save(mats,file="matsMean")

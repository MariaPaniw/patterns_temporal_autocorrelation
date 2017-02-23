# Script for Paniw et al. XXXXXX - Appendix S1
#This script peruses through COMPADRE and COMADRE, it collapses matrices into a given number of stages, and outputs vital rates from each matrix
#Author: Rob Salguero-Gomez & Maria Paniw
#Created: 12 May 2016

#Clean memory
rm(list=ls(all=TRUE))

#Functions


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
	
	
#Example:
#collapse1=c("1-2","3-4","5")
#collapse2=c("1-2","3-4-5")
#collapse3=c("1-2-3-4-5")
#matA=matrix(c(0.1,0.0,0.0,5.4,8.3,
#			  0.2,0.2,0.1,0.0,0.0,
#			  0.1,0.4,0.2,0.2,0.0,
#			  0.0,0.1,0.5,0.6,0.2,
#			  0.0,0.1,0.5,0.1,0.7),nrow=5,byrow=T)
#collapseMatrix(matA,collapse)
#eigen(collapseMatrix(matA,collapse1))
#eigen(collapseMatrix(matA,collapse2))
#eigen(collapseMatrix(matA,collapse3))


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



#Function to extract the vital rates of a 3x3 matrix that has been collapsed into pre-, reprod, and post-reprod

extractVitalRates <- function(matU,matF){
  if(dim(matU)[1]!=4) {print("This matrix is not 4x4!")
    stop}
  vitalRates=rep(NA,10)
  surv=colSums(matU,na.rm=T)
  surv[which(is.na(collapse))]=NA
    s1=surv[1]
    s2=surv[2]
    s3=surv[3]
    s4=surv[4]
  matUIndep=matU
  for (i in 1:dim(matU)[1]){matUIndep[,i]=matU[,i]/surv[i]}
    g21=matUIndep[2,1]
    g31=matUIndep[3,1]
    g32=matUIndep[3,2]
    g43=matUIndep[4,3]
    g34=matUIndep[3,4]
    
  f13=matF[1,3]
  f23=matF[2,3]
  f33=matF[3,3]

  vitalRates=data.frame(s1,s2,s3,s4,g21,g31,g32,g43,g34,f13,f23,f33)
  colnames(vitalRates)=c("s1","s2","s3","s4","g21","g31","g32","g43","g34","f13","f23","f33")

  return(vitalRates)    
}


#Code to obtain the vital rates of pre-reproducing, reproducing and non-reproducing
#individuals after matrices have been collapsed

# Set the working directory, then load the COMPADRE data:
dir <- setwd("C:/Users/Maria/Dropbox/TempAutoProject") # CHANGE THIS TO YOUR DIRECTORY
load(paste(dir,"/COMPADRE_v.4.0.0.RData",sep=""))
load(paste(dir, "/COMADRE_v.2.0.0.RData", sep=""))


#Subsetting all matrices available in COMPADRE (plants/algae)

indexCOMPADRE=which(duplicated(compadre$metadata$SpeciesAuthor)==FALSE &   # mean matrices (across sites and times)
                      compadre$metadata$SurvivalIssue<=1 &                 # survival does not exceed 1 (may be due to rounding errors)
                      compadre$metadata$MatrixFec=="Yes" &                 # non-zero fecundities
                      compadre$metadata$MatrixSplit=="Divided" &           # fecundity and transition matrices available separately
                      compadre$metadata$MatrixDimension>3 &
                      compadre$metadata$MatrixTreatment=="Unmanipulated")  # only natural (not experimental or laboratory) sites 

#Discard clonality
compadreClonIndex=NULL
for (i in 1:dim(compadre$metadata)[1]){
  if(sum(compadre$mat[[i]]$matC,na.rm=T)>0){compadreClonIndex[i]=1}else{compadreClonIndex[i]=0}
}
compadreClonIndex=which(compadreClonIndex==1)

indexCOMPADRE=indexCOMPADRE[which(!indexCOMPADRE%in%compadreClonIndex)]

compadre2=NULL
compadre2$metadata=compadre$metadata[indexCOMPADRE,]
compadre2$mat=compadre$mat[indexCOMPADRE]
compadre2$matrixClass=compadre$matrixClass[indexCOMPADRE]

#Subsetting all matrices available in COMADRE (animals) - same criteria as for COMPADRE

indexCOMADRE=which(duplicated(comadre$metadata$SpeciesAuthor)==FALSE &
                      comadre$metadata$SurvivalIssue<=1 &
                      comadre$metadata$MatrixFec=="Yes" &
                      comadre$metadata$MatrixSplit=="Divided" &
                      comadre$metadata$MatrixDimension>3 &
                      comadre$metadata$MatrixTreatment=="Unmanipulated")

#Discard clonality
comadreClonIndex=NULL
for (i in 1:dim(comadre$metadata)[1]){
  if(sum(comadre$mat[[i]]$matC,na.rm=T)>0){comadreClonIndex[i]=1}else{comadreClonIndex[i]=0}
}
comadreClonIndex=which(comadreClonIndex==1)

indexCOMADRE=indexCOMADRE[which(!indexCOMADRE%in%comadreClonIndex)]

comadre2=NULL
comadre2$metadata=comadre$metadata[indexCOMADRE,]
comadre2$mat=comadre$mat[indexCOMADRE]
comadre2$matrixClass=comadre$matrixClass[indexCOMADRE]


#Loop through the subset of COM(P)ADRE matrices and metadata information 
# to obtain relevant vital rates and information on study details and habitat types:

dims=dim(comadre2$metadata)[1]+dim(compadre2$metadata)[1]

# create an empty data frame to hold results 
output=data.frame("SpeciesAuthor"=rep(NA,dims),
                  "SpeciesAccepted"=rep(NA,dims),
                  "Family"=rep(NA,dims),
                  "Kingdom"=rep(NA,dims),
                  "Authors"=rep(NA,dims),
                  "Journal"=rep(NA,dims),
                  "YearPublication"=rep(NA,dims),
                  "StudyDuration"=rep(NA,dims),
                  "habitat"=rep(NA,dims),
                  "s1"=rep(NA,dims),
                  "s2"=rep(NA,dims),
                  "s3"=rep(NA,dims),
                  "s4"=rep(NA,dims),
                  "g21"=rep(NA,dims),
                  "g31"=rep(NA,dims),
                  "g32"=rep(NA,dims),
                  "g43"=rep(NA,dims),
                  "g34"=rep(NA,dims),
                  "f13"=rep(NA,dims),
                  "f23"=rep(NA,dims),
                  "f33"=rep(NA,dims))

count=0
for (loop in 1:2){
  if (loop==1){d=compadre2}else{d=comadre2}
  
for (i in 1:dim(d$metadata)[1]){
  tryCatch({
    count=count+1
    
    output$SpeciesAuthor[count]=d$metadata$SpeciesAuthor[i]
    output$SpeciesAccepted[count]=d$metadata$SpeciesAccepted[i]
    output$Family[count]=as.character(d$metadata$Family[i])
    output$Kingdom[count]=as.character(d$metadata$Kingdom[i])
    output$Authors[count]=d$metadata$Authors[i]
    output$Journal[count]=d$metadata$Journal[i]
    output$YearPublication[count]=d$metadata$YearPublication[i]
    output$StudyDuration[count]=d$metadata$StudyDuration[i]
    output$habitat[count]=d$metadata$Ecoregion[i]
    
    rm(matUcollapse,matFcollapse)
    matU=d$mat[[i]]$matU
    matF=d$mat[[i]]$matF
      matU=rearrangeMatrix(matU,matF)[[1]]$matU
      matF=rearrangeMatrix(matU,matF)[[1]]$matF
    matrixStages=d$matrixClass[[i]]$MatrixClassOrganized
    collapse=reprodStages(matF,matrixStages)
    
    matUcollapse=collapseMatrix(matU,matF,collapse=collapse)$matU
    matUcollapse[is.na(collapse),is.na(collapse)]=NA
    
    matFcollapse=collapseMatrix(matU,matF,collapse=collapse)$matF
    matFcollapse[is.na(collapse),is.na(collapse)]=NA
    
    output[count,c("s1","s2","s3","s4","g21","g31","g32","g43","g34","f13","f23","f33")]=extractVitalRates(matU=matUcollapse,matF=matFcollapse)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
}

output[is.na(output)]=NA
vitalRates=c("s1","s2","s3","s4","g21","g31","g32","g43","g34","f13","f23","f33")

# elimate cases where all vital rates are NAs
eliminate=which(rowSums(is.na(output[,vitalRates]))==12)
output=output[-eliminate,]
dim(output)

#Capping survival to max 1
output$s1[which(output$s1>1)]=1
output$s2[which(output$s2>1)]=1
output$s3[which(output$s3>1)]=1
output$s4[which(output$s4>1)]=1

# Remove studies that analyze the same species from the same or similar sites 

eliminate2=c("Tillandsia_multicaulis","Astragalus_scaphoides","Calochortus_lyallii_2",
             "Cirsium_pitcheri_3","Cirsium_pitcheri_6","Dactylorhiza_lapponica",
             "Erythronium_japonicum","Erythronium_japonicum_3","Haplopappus_radiatus",
             "Lepanthes_rupestris","Lomatium_bradshawii_2","Panax_quinquefolius","Panax_quinquefolius_4",
             "Primula_veris_4","Sarracenia_purpurea_2","Succisa_pratensis_2",
             "Trollius_laxus","Pseudophoenix_sargentii_2", "Lantana_camara_2","Astrophytum_ornatum_2",
             "Mammillaria_huitzilopochtli_2","Mammillaria_pectinifera_4","Neobuxbaumia_mezcalaensis_2",
             "Neobuxbaumia_tetetzo_2","Neobuxbaumia_tetetzo_3","Neobuxbaumia_tetetzo_4",
             "Pinus_lambertiana","Shorea_leprosula_2","Maccullochella_peelii",
             "Macquaria_ambigua","Oncorhynchus_tshawytscha_2","Oncorhynchus_tshawytscha_4",
             "Pimephales_promelas_4","Pylodictis_olivaris_2","Paramuricea_clavata_2",
             "Coragyps_atratus","Sterna_hirundo","Daphnia_pulex","Scolytus_ventralis",
             "Palaemonetes_pugio","Acinonyx_jubatus","Alces_alces_6","Alces_alces_10","Canis_lupus",
             "Cervus_elaphus_4","Eumetopias_jubatus_3","Halichoerus_grypus","Macaca_mulatta",
             "Macaca_mulatta_3","Marmota_flaviventris_3","Phocarctos_hookeri","Sigmodon_hispidus",
             "Ursus_americanus_3","Ursus_americanus_5","Caretta_caretta_3","Chrysemys_picta_2",
             "Chrysemys_picta_3","Crocodylus_johnstoni_3","Xenosaurus_sp.")

output=output[!output$SpeciesAuthor%in%eliminate2,]

# collapse habitat types to 5 categories:

habitat.original=c("TMB","TMB; TDB","TSC; TDB","TDB","TSC","TGV","MAN","DES","DES; MED",
                   "DES; TMB","MED","TBM","TCF","TBM; TSC","TMB; TDB","TMB; TGV",
                   "TBM; TCF","TCF; DES","MED; TBM","TBM; BOR","TGS","FGS","MON",
                   "BOR","TUN","BOR; TUN","TUN; BOR","TEU","POE","LRD","LLE","LRE",
                   "SRE","TRU","TSS","TRC","TMB; DES; MON","TBM; TCF; BOR","TMB; TGV; TDB","SLE")

habitat.collapsed=c("Tropical & Subtropical","Tropical & Subtropical","Tropical & Subtropical",
                    "Tropical & Subtropical","Tropical & Subtropical","Tropical & Subtropical",
                    "Tropical & Subtropical","Arid","Arid","Arid","Temperate","Temperate",
                    "Temperate","Temperate","Temperate","Temperate","Temperate","Temperate",
                    "Temperate","Alpine","Temperate","Temperate","Alpine & Arctic","Alpine & Arctic",
                    "Alpine & Arctic","Alpine & Arctic","Alpine & Arctic","Aquatic","Aquatic",
                    "Aquatic","Aquatic","Aquatic","Aquatic","Aquatic","Aquatic","Aquatic",
                    "Temperate","Temperate","Tropical & Subtropical","Temperate","Arid",
                    "Arid","Aquatic","Temperate","Aquatic","Aquatic","Aquatic",
                    "Temperate","Temperate","Tropical & Subtropical","Aquatic")

for(h in 1:length(habitat.original)){
 
  output$habitat[output$habitat==habitat.original[h]] = habitat.collapsed[h] 
}


# some habitats had to be input manually based on a literature seach:

extraHabitat=read.csv("extraHabitat.csv")

output$habitat[output$SpeciesAuthor%in%extraHabitat$SpeciesAuthor]=as.character(extraHabitat$habitat)
  
output$habitat[output$habitat=="LAB"]=c("Temperate","Aquatic","Tropical & Subtropical","Temperate")

#### 

# Subset data to retain only populations for which we could obtain phylogenetic information

phylo=read.csv("phyloSpecies.csv") 

output=output[output$SpeciesAuthor%in%phylo$SpeciesAuthor,]

nrow(output)# should be 428

# save file

write.csv(output,"vitalRates.csv") # This file is already exists in case of issues with this code


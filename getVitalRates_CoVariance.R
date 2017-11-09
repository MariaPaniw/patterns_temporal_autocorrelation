# Script for Paniw et al. XXXXXX - Appendix S1
#This script peruses through COMPADRE and COMADRE and outputs vital rates, vital rate classes, and vital rate correlation matrix for 109 sub species
#Author: Maria Paniw
#Created: 19 Aug 2011


#Clean memory
rm(list=ls(all=TRUE))

library(stringr)
library(plyr)

# Set the working directory, then load the COMPADRE data:
dir <- setwd("/Users/mariapaniw/Dropbox/TempAutoProject/SuppMat") # CHANGE THIS TO YOUR DIRECTORY
load(paste(dir,"/COMPADRE_v.4.0.0.RData",sep=""))
load(paste(dir, "/COMADRE_v.2.0.0.RData", sep=""))

# load average vital rates
load("matsMean")
# get IDs of species:
sp.data=read.csv("phyloSpecies.csv")

# 109 species with at least 3 annual matrices

ID=c("Primula_elatior","Eryngium_cuneifolium" ,"Agrimonia_eupatoria","Coryphantha_robbinsorum","Petrocoptis_pseudoviscosa_2","Papio_cynocephalus",                 
     "Oenothera_deltoides" ,"Ardisia_elliptica","Lotus_arinagensis", "Mammillaria_napina","Taxus_floridana","Eryngium_alpinum","Propithecus_verreauxi","Cleistes_divaricata_var._bifaria",
     "Cleistes_divaricata_var._divaricata","Phyllanthus_indofischeri","Mimulus_cardinalis","Xenosaurus_grandis" ,"Geum_rivale","Limonium_geronense","Erodium_paularense",                 
      "Mimulus_lewisii","Arabis_fecunda","Atriplex_acanthocarpa","Atriplex_canescens" ,"Astragalus_peckii" , "Sapium_sebiferum",                   
      "Helianthemum_polygonoides","Antirrhinum_lopesianum","Rumex_rupestris","Castanea_dentata","Cytisus_scoparius","Purshia_subintegra", "Calochortus_lyallii","Limonium_malacitanum",   "Xenosaurus_platyceps",               
      "Cimicifuga_elata","Silene_spaldingii","Dicerandra_frutescens","Asplenium_adulterinum","Asplenium_cuneifolium", "Polemonium_van-bruntiae",            
      "Lathyrus_vernus","Pyrrocoma_radiata","Cirsium_vulgare_3", "Cimicifuga_rubifolia","Silene_acaulis","Umbonium_costatum" ,                 
      "Astroblepus_ubidiai","Ramonda_myconi","Cercopithecus_mitis", "Primula_farinosa" ,"Gorilla_beringei","Dioon_caputoi","Anser_caerulescens","Orchis_purpurea","Liatris_scariosa",                   
      "Abies_concolor","Abies_magnifica","Mammillaria_hernandezii_2","Lomatium_bradshawii","Horkelia_congesta","Cecropia_obtusifolia",               
      "Oxytropis_jabalambrensis","Astragalus_scaphoides_2","Primula_veris_2" ,"Armeria_merinoi","Lomatium_cookii","Pinus_strobus",                      
      "Succisa_pratensis_3","Scolytus_ventralis_2","Euterpe_edulis","Anthropoides_paradiseus","Cirsium_palustre", "Lupinus_tidestromii", "Orcinus_orca_2", "Cypripedium_calceolus","Shorea_leprosula",                   
      "Phyllanthus_emblica_3", "Molinia_caerulea", "Calathea_ovandensis","Paramuricea_clavata","Cryptantha_flava","Cypripedium_fasciculatum",           
      "Callospermophilus_lateralis","Colias_alexandra","Brachyteles_hypoxanthus","Mammillaria_huitzilopochtli","Catopsis_compacta","Tillandsia_violacea",                
      "Cebus_capucinus","Santolina_melidensis" ,"Astragalus_tremolsianus" ,"Zea_diploperennis","Astragalus_tyghensis","Actaea_spicata",                     
      "Plantago_media","Ovis_aries_2","Calocedrus_decurrens","Neobuxbaumia_macrocephala","Neobuxbaumia_mezcalaensis","Neobuxbaumia_tetetzo",               
      "Helianthemum_juliae","Vella_pseudocytisus_subsp._paui","Cirsium_pitcheri_8","Ambloplites_rupestris_2","Cottus_bairdi","Etheostoma_flabellare_2",            
      "Tillandsia_macdougallii")

sp.data=droplevels(sp.data[sp.data$SpeciesAuthor%in%ID&sp.data$MatrixComposite=="Individual",])


## take out problematic populations

sp.data=droplevels(sp.data[-which(sp.data$MatrixPopulation%in%c("Transitional Forest","Morningside Nature Center","Bottomland hardwood forest","Young mixed pine-hardwood forest",
                                  "Bull Flat","Campion Crest","Pass","Ridge","Wawona","Abbotts Langdon","Haynes Creek","Sheep Corral Gulch",
                                  "La Pedrera","Plot 1","Plot 2","Plot 4","Site 4","Schinus thicket","S2")),])


### FOR REAL VITAL RATES
matsVarCov=vector("list", length(unique(sp.data$SpeciesAuthor)))

for(i in 1:length(unique(sp.data$SpeciesAuthor))){
  
  # subset species and sites
  sp=as.character(unique(sp.data$SpeciesAuthor)[i])
  site=as.character(unique(sp.data$MatrixPopulation[sp.data$SpeciesAuthor==sp]))
  
  ## periodicity
  
  per=sp.data$AnnualPeriodicity[sp.data$SpeciesAuthor==sp][1]
  
  # empty vector to hold variance (varvar) and covariance (varcov)
  varvar=vector("list", length(site))
  varcov=vector("list", length(site))
  
  for(j in 1:length(site)){
    
    # FOR PLANTS/ALGAE
    if(sp.data$Kingdom[sp.data$SpeciesAuthor==sp][1]=="Plantae"|sp.data$Kingdom[sp.data$SpeciesAuthor==sp][1]=="Chromalveolata"){
      
      index=which(compadre$metadata$SpeciesAuthor==sp &
                    compadre$metadata$MatrixPopulation==site[j] &
                    compadre$metadata$MatrixComposite == "Individual") 
      
      matsVar=compadre$mat[index]
    
      # mean MPM for species
      
      indexMU=which(duplicated(compadre$metadata$SpeciesAuthor)==FALSE&
                      compadre$metadata$SpeciesAuthor==sp) 
      matU=compadre$mat[indexMU][[1]]$matU
      matF=compadre$mat[indexMU][[1]]$matF
    
      # FOR ANIMALS  
    }else{
      index=which(comadre$metadata$SpeciesAuthor==sp &
                    comadre$metadata$MatrixPopulation==site[j] &
                    comadre$metadata$MatrixComposite == "Individual") 
      
      matsVar=comadre$mat[index]
      # mean MPM for species
      
      indexMU=which(duplicated(comadre$metadata$SpeciesAuthor)==FALSE&
                      comadre$metadata$SpeciesAuthor==sp) 
      matU=comadre$mat[indexMU][[1]]$matU
      matF=comadre$mat[indexMU][[1]]$matF
    }
    
    # fix matrices for individual species
    if(sp=="Orcinus_orca_2") matsVar<-matsVar[-25]
    if(sp=="Cirsium_vulgare_3") matsVar<-matsVar[-c(1,4,8)]
    if(sp=="Colias_alexandra"){
      
      matsVar<-matsVar[-c(5,7)]
      matsVar[[1]]$matF[1,7]=0
    } 
    
    
    # Average reproduction across years (to deal with 0 fecundities for some years)
    Fec.mu=matrix(0,dim(matsVar[[1]]$matF)[1],dim(matsVar[[1]]$matF)[1])
    for(bb in 1:length(matsVar)){
      
      Fec.mu=Fec.mu+matsVar[[bb]]$matF
    }
    
    # get vital rates per MPM per site per species
    vr.all=NULL
    for(b in 1:length(matsVar)){
      
      ### survival
      surv=colSums(matsVar[[b]]$matU)^per # account for periodicity
      if(surv[length(surv)]>0.99999) surv[length(surv)]<-0.995 # prevent survival of last stage/age to be 1 - it will make simulations unstable
      
      names(surv)=paste("s",1:length(surv),sep="")
      
      U.mat=matsVar[[b]]$matU
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
      
      # reproduction (if MPM has 0 reproduction)
      if(length(which(matsVar[[b]]$matF>0))==0){
        placeholder=matrix(1:length(as.numeric(matsVar[[b]]$matF)),dim(matsVar[[b]]$matF)[1],dim(matsVar[[1]]$matF)[1])
        colnames(placeholder)=rownames(placeholder)=1:dim(matsVar[[b]]$matF)[1]
        fec.names=which(Fec.mu>0)
        names=expand.grid(rownames(placeholder),colnames(placeholder))[placeholder%in%fec.names,]
        names=interaction(str_pad(as.numeric(names$Var1),2,pad="0"),str_pad(as.numeric(names$Var2),2,pad="0"),sep="")
        fec=paste("f",names ,sep="")
        fec.value=rep(0,length(fec)) 
        names(fec.value)=fec
        
        
      }else{
        placeholder=matrix(1:length(as.numeric(matsVar[[b]]$matF)),dim(matsVar[[b]]$matF)[1],dim(matsVar[[1]]$matF)[1])
        colnames(placeholder)=rownames(placeholder)=1:dim(matsVar[[b]]$matF)[1]
        fec.names=which(Fec.mu>0)
        names=expand.grid(rownames(placeholder),colnames(placeholder))[placeholder%in%fec.names,]
        names2=as.numeric(rownames((names)))
        names=interaction(str_pad(as.numeric(names$Var1),2,pad="0"),str_pad(as.numeric(names$Var2),2,pad="0"),sep="")
        fec=paste("f",names ,sep="")
        fec.value=as.numeric(matsVar[[b]]$matF)[names2]*per #account for periodicity
        names(fec.value)=fec
        
      }
      
      vr=matrix(c(surv,gr,ret,fec.value),ncol=length(c(surv,gr,ret,fec.value)))
      colnames(vr)=names(c(surv,gr,ret,fec.value))
      if(b==1){
        vr.all=rbind(vr.all,vr)
      }else{vr.all=rbind.fill.matrix(vr.all,vr) }
      
    }
    
    # for each site: 
    varcov.sub=cor(vr.all,method="spearman") # correlation
    varcov.sub[is.na(varcov.sub)]=0
    var.sub=diag(var(vr.all))# variance
    var.sub[is.na(var.sub)]=0
    
    varcov[[j]]=varcov.sub
    varvar[[j]]=var.sub
    
  }
  # take mean correlation across sites
  varcov.a=array(unlist(varcov), dim = c(nrow(varcov[[1]]), ncol(varcov[[1]]), length(varcov)))
  varcov.mu=apply(varcov.a,c(1,2),mean,na.rm=T)
  colnames(varcov.mu)=rownames(varcov.mu)=colnames(varcov.sub)
  
  # take mean variance across sites
  varvar.a=array(unlist(varvar), dim = c(1, length(varvar[[1]]), length(varvar)))
  varvar.mu=as.numeric(apply(varvar.a,c(1,2),mean,na.rm=T))
  names(varvar.mu)=names(var.sub)
  
  # remove vital rates with 0 variance of correlation 
  if(any(varvar.mu==0)){
    sub=varcov.mu[-which(varvar.mu==0),-which(varvar.mu==0)]
    sub2=varvar.mu[-which(varvar.mu==0)]
  }else{
    sub=varcov.mu
    sub2=varvar.mu
  }
  
  matsVarCov[[i]]$var=sub2
  
  matsVarCov[[i]]$corr=sub
  matsVarCov[[i]]$matU=matU
  matsVarCov[[i]]$matF=matF
  
  matsVarCov[[i]]$vr.mu=mats[[which(sapply(lapply(mats, function(ch) grep(sp, ch)), function(x) length(x) > 0))]]$vr
  matsVarCov[[i]]$species=sp
}

# save results

save(matsVarCov,file="matsVarCov")

########

corpearson<-list()
Ccorpearson<-list()
resultp<-list()
for(m in 1:t){
  corpearson[[m]]<-matrix(0,dim(sigma)[1],dim(sigma)[2])
  for(i in 1:nrow(sigma)){
    for(j in 1:nrow(sigma)){
      if(sigma[i,j]!=0){
        data.x<-newdata[[m]][,i]
        data.y<-newdata[[m]][,j]
        corpearson[[m]][i,j]<-cor(data.x,data.y,method = "pearson")
      }
    }
  }
  Ccorpearson[[m]]<-abs(corpearson[[m]])
  for(i in 1:nrow(corpearson[[m]])){
    for(j in 1:nrow(corpearson[[m]])){
      Ccorpearson[[m]][i,j]<-log(1/Ccorpearson[[m]][i,j])
    }
  }
  
  for(i in 1:nrow(corpearson[[m]])){
    for(j in 1:nrow(corpearson[[m]])){
      if(Ccorpearson[[m]][i,j]==0){
        Ccorpearson[[m]][i,j]<-Inf
      }
    }
  }
  resultp[[m]]<-kshortestpath( Ccorpearson[[m]],start,destination,paths)
  
}

################
corspearman<-list()
Ccorspearman<-list()
resultS<-list()
for(m in 1:t){
  corspearman[[m]]<-matrix(0,dim(sigma)[1],dim(sigma)[2])
  for(i in 1:dim(sigma)[1]){
    for(j in 1:dim(sigma)[1]){
      if(sigma[i,j]!=0){
        data.x<-newdata[[m]][,i]
        data.y<-newdata[[m]][,j]
        corspearman[[m]][i,j]<-cor(data.x,data.y,method = "spearman")
        if(corspearman[[m]][i,j]==0){corspearman[[m]][i,j]=0.000001}
      }
    }
  }
  Ccorspearman[[m]]<-corspearman[[m]]
  
  for(i in 1:nrow(corspearman[[m]])){
    for(j in 1:nrow(corspearman[[m]])){
      if(corspearman[[m]][i,j]<0){
        Ccorspearman[[m]][i,j]<-abs(Ccorspearman[[m]][i,j])
      }
      Ccorspearman[[m]][i,j]<-log(1/Ccorspearman[[m]][i,j])
    }
  }
  for(i in 1:nrow(corspearman[[m]])){
    for(j in 1:nrow(corspearman[[m]])){
      if(Ccorspearman[[m]][i,j]==0){
        Ccorspearman[[m]][i,j]<-Inf
      }
    }
  }
  
  resultS[[m]]<-kshortestpath( Ccorspearman[[m]],start,destination,paths)
}


#############################3
#install.packages("energy")
library(energy)
cordistance<-list()
Ccordistance<-list()
resultd<-list()
for(m in 1:t){
  cordistance[[m]]<-matrix(0,dim(sigma)[1],dim(sigma)[2])
  for(i in 1:dim(sigma)[1]){
    for(j in 1:dim(sigma)[1]){
      if(sigma[i,j]!=0){
        data.x<-newdata[[m]][,i]
        data.y<-newdata[[m]][,j]
        cordistance[[m]][i,j]<-dcor(data.x,data.y)
        
      }
    }
  }
  Ccordistance[[m]]<-cordistance[[m]]
  
  for(i in 1:nrow(cordistance[[m]])){
    for(j in 1:nrow(cordistance[[m]])){
      Ccordistance[[m]][i,j]<-log(1/Ccordistance[[m]][i,j])
    }
  }
  for(i in 1:nrow(cordistance[[m]])){
    for(j in 1:nrow(cordistance[[m]])){
      if(Ccordistance[[m]][i,j]==0){
        Ccordistance[[m]][i,j]<-Inf
      }
    }
  }
  
  resultd[[m]]<-kshortestpath( Ccordistance[[m]],start,destination,paths)
  cat(m, "\n")
}
Sys.time()
###################MIC5
#install.packages("minerva")
library(minerva)
corMIC<-list()
CcorMIC<-list()
resultMIC<-list()
for(m in 1:t){
  corMIC[[m]]<-matrix(0,dim(sigma)[1],dim(sigma)[2])
  for(i in 1:dim(sigma)[1]){
    for(j in 1:dim(sigma)[1]){
      if(sigma[i,j]!=0){
        data.x<-10*newdata[[m]][,i]
        data.y<-10*newdata[[m]][,j]
        corMIC[[m]][i,j]<-mine(data.x,data.y,alpha=0.3)$MIC
        
      }
    }
  }
  CcorMIC[[m]]<-corMIC[[m]]
  
  for(i in 1:nrow(corMIC[[m]])){
    for(j in 1:nrow(corMIC[[m]])){
      CcorMIC[[m]][i,j]<-log(1/CcorMIC[[m]][i,j])
    }
  }
  for(i in 1:nrow(corMIC[[m]])){
    for(j in 1:nrow(corMIC[[m]])){
      if(CcorMIC[[m]][i,j]==0){
        CcorMIC[[m]][i,j]<-Inf
      }
    }
  }
  
  resultMIC[[m]]<-kshortestpath( CcorMIC[[m]],start,destination,paths)
}

###################################################6

source("KDE.r")
jieguo<-list()
cormi<-list()
Ccormi<-list()
resultmi<-list()
sigmaI<-sigma
sigmaI[lower.tri(sigmaI,diag=T)]=0
edge <- which(sigmaI!=0,arr.ind=T)

for(m in 1:t){
  
  cormi[[m]]<-diag(rep(1,nrow(sigma)))
  jieguo[[m]]<-apply(edge,DataFit=newdata[[m]],1,denPre2D)
  for(i in 1:nrow(edge)){
    cormi[[m]][edge[i,1],edge[i,2]]<-mean(jieguo[[m]][,i])
  }
  Ccormi[[m]]<-cormi[[m]]
  for(i in 1:nrow(cormi[[m]])){
    for(j in 1:nrow(cormi[[m]])){
      if(Ccormi[[m]][i,j]<0){
        Ccormi[[m]][i,j]<-abs(Ccormi[[m]][i,j])
      }
      Ccormi[[m]][i,j]<-log(1/Ccormi[[m]][i,j])
    }
  }
  for(i in 1:nrow(cormi[[m]])){
    for(j in 1:nrow(cormi[[m]])){
      if(Ccormi[[m]][i,j]==0){
        Ccormi[[m]][i,j]<-Inf
      }
      if(Ccormi[[m]][i,j]<0){
        Ccormi[[m]][i,j]<-abs(Ccormi[[m]][i,j])
      }
    }
  }
  
  resultmi[[m]]<-kshortestpath( Ccormi[[m]],start,destination,paths)
  
  cat(m, "\n")
}
Sys.time()
#mcc
library(acepack)
corMCC<-list()
CcorMCC<-list()
resultMCC<-list()
for(m in 1:t){
  corMCC[[m]]<-matrix(0,dim(sigma)[1],dim(sigma)[2])
  for(i in 1:nrow(sigma)){
    for(j in 1:nrow(sigma)){
      if(sigma[i,j]!=0){
        data.x<-newdata[[m]][,i]
        data.y<-newdata[[m]][,j]
        argmax = ace(data.x, data.y)
        corMCC[[m]][i,j]<-cor(argmax$tx, argmax$ty)
      }
    }
  }
  CcorMCC[[m]]<-abs(corMCC[[m]])
  for(i in 1:nrow(corMCC[[m]])){
    for(j in 1:nrow(corMCC[[m]])){
      CcorMCC[[m]][i,j]<-log(1/CcorMCC[[m]][i,j])
    }
  }
  
  for(i in 1:nrow(corMCC[[m]])){
    for(j in 1:nrow(corMCC[[m]])){
      if(CcorMCC[[m]][i,j]==0){
        CcorMCC[[m]][i,j]<-Inf
      }
    }
  }
  resultMCC[[m]]<-kshortestpath( CcorMCC[[m]],start,destination,paths)
  
}

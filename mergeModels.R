mergeModels=function(MV1,MV2) {
  if(any(class(MV1)=='Multilevel')|any(class(MV1)=='Classification')) {
    cat('\nNot yet supported')
    break()
  }
  in1=MV1$inData
  in2=MV2$inData
  nRep1=MV1$inData$nRep
  nRep2=MV2$inData$nRep
  nRep=nRep1+nRep2
  in1$nRep=NULL
  in2$nRep=NULL
  if(!identical(in1,in2)) {
    cat('\nIndata not identical between models')
    stop()
  }
  DA=MV1$inData$DA
  PLS=MV1$inData$method=='PLS'
  yP=MV1$yPred
  yPPR=MV1$yPredPerRep
  VIP=MV1$VIP
  VIPrep=MV1$VIPPerRep
  nV=MV1$nVar
  nVPR=MV1$nVarPerRep
  if (PLS) {
    nC=MV1$nComp
    nCPR=MV1$nCompPerRep
  }
  for(i in 1:3) {
    if (DA) {
      cat('\nNot yet implemented')
    } else {
      yPPR[[i]]=cbind(yPPR[[i]],MV2$yPredPerRep[[i]])
      yP[,i]=apply(yPPR[[i]],1,mean)
      VIPrep[[i]]=cbind(VIP[[i]],MV2$VIPPerRep[[i]])
      VIP[,i]=apply(VIPrep[[i]],1,mean)
      nVPR[[i]]=c(nVPR[[i]],MV2$nVarPerRep[[i]])
      nV[i]=mean(nVPR[[i]])
      if(PLS) {
        nCPR[[i]]=c(nCPR[[i]],MV2$nCompPerRep[[i]])
        nC[i]=mean(nCPR[[i]])
      }
    }
  }
  newMod=list()
  newMod$inData=in1
  newMod$inData$nRep=nRep
  newMod$yPred=yP
  newMod$yPredPerRep=yPPR
  newMod$VIP=VIP
  newMod$VIPPerRep=VIPrep
  newMod$nVar=nV
  newMod$nVarPerRep=nVPR
  if (PLS) {
    newMod$nComp=nC
    newMod$nCompPerRep=nCPR
  }
  class(newMod)=c(class(MV1),'Merged')
  return(newMod)
}
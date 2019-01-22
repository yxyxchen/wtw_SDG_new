simulation = function(modelName, nBlock){
  # library 
  library('ggplot2')
  library('dplyr')
  library('tidyr')
  source('subFxs/simulationFxs.R') # 
  load("wtwSettings.RData")
  source("subFxs/taskFxs.R")
  fileName = sprintf('genData/simulation/%s/realParas.RData', modelName)
  load(fileName)
  
  # generate outfile
  outFile = sprintf('genData/simulation/%s', modelName)
  dir.create(outFile)
  
  # choose modelFun
  modelFun = getSimModelFun(modelName)
  ################ simulation ################
  for(condIdx in 1 : 2){
    cond = conditions[condIdx];
    tMax = tMaxs[condIdx];
    # set seed
    set.seed(123)
    # initialize outputs
    trialData = vector(length = nComb * nRep, mode ='list')
    
    # loop over repetions 
    for(h in 1 : nrow(realParas)){
      para = realParas[h,];
      # calculate wIni
      for(j in 1 : nRep ){
        tempt = modelFun(para, cond, nBlock)
        trialData[[simIdx[h, j]]] = tempt
      }  
    }
    
    if(cond == "HP"){
      trialHPData = trialData
      fileName = sprintf('genData/simulation/%s/trialHPData.RData', modelName)
      save(trialHPData,file = fileName)
    }else{
      trialLPData = trialData
      fileName =  sprintf('genData/simulation/%s/trialLPData.RData', modelName)
      save(trialLPData,file = fileName)
    }
  }
}




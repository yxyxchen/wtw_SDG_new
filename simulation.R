simulation = function(modelName){
  # library 
  library('ggplot2')
  library('dplyr')
  library('tidyr')
  source('subFxs/simulationFxs.R') # 
  load("wtwSettings.RData")
  fileName = sprintf('genData/simulation/%s/realParas.RData', modelName)
  load(fileName)
  
  # generate outfile
  outFile = sprintf('genData/simulation/%s', modelName)
  dir.create(outFile)
  
  # choose modelFun
  modelFun = getModelFun(modelName)
  ################ simulation ################
  nRep = 5
  count = t(matrix(1 : (nComb * nRep), nRep, nComb))
  
  for(condIdx in 1 : 2){
    cond = conditions[condIdx];
    tMax = tMaxs[condIdx];
    # set seed
    set.seed(123)
    # initialize outputs
    TrialEarnings = array(dim = c(nValue^nPara, nRep, blockSecs / iti + 1))
    RewardDelays = array(dim = c(nValue^nPara, nRep, blockSecs / iti + 1))
    TimeWaited = array(dim = c(nValue^nPara, nRep, blockSecs / iti + 1))
    vaQuits = array(dim = c(nValue^nPara, nRep, blockSecs / iti + 1))
    vaWaits = array(dim = c(nValue^nPara, nRep, tMax / stepDuration, blockSecs / iti + 1))
    thisPackData = vector(length = nComb * nRep, mode ='list')
    
    # loop over repetions 
    for(h in 1 : nrow(realParas)){
      para = realParas[h,];
      # calculate wIni
      for(j in 1 : nRep ){
        tempt = modelFun(para, cond)
        TrialEarnings[h, j,] = tempt[['trialEarnings']]
        RewardDelays[h, j,] = tempt[['rewardDelays']]
        TimeWaited[h, j, ] = tempt[['timeWaited']]
        vaQuits[h, j, ] = tempt[['vaQuits']]
        vaWaits[h, j, ,  ] = tempt[['vaWaits']]
        thisPackData[[count[h, j]]] = tempt
      }  
    }
    
    # organize and save outputs 
    outputData = list("timeWaited" = TimeWaited,
                      "rewardDelays" = RewardDelays, "trialEarnings" = TrialEarnings,
                      "vaWaits" = vaWaits, "vaQuits" = vaQuits
    )
    
    outFile = 'simData'
    if(cond == "HP"){
      rawHPData = outputData
      packHPData = thisPackData
      fileName = sprintf('genData/simulation/%s/rawHPData.RData', modelName)
      save(rawHPData,file = fileName) 
      fileName = sprintf('genData/simulation/%s/packHPData.RData', modelName)
      save(packHPData,file = fileName)
    }else{
      rawLPData = outputData
      packLPData = thisPackData
      fileName = sprintf('genData/simulation/%s/rawLPData.RData', modelName)
      save(rawLPData,file = fileName)
      fileName =  sprintf('genData/simulation/%s/packLPData.RData', modelName)
      save(packLPData,file = fileName)
    }
  }
}




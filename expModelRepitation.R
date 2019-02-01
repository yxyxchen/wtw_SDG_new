library("ggplot2")
source("subFxs/plotThemes.R")
load("genData/expDataAnalysis/blockData.RData")
blockData = blockData[blockData$blockNum == 1,]
source("subFxs/simulationFxs.R") # used in simulation
source("subFxs/taskFxs.R") # used in simulation
source("subFxs/analysisFxs.R") # for analysis
load("wtwSettings.RData") # used in simulation 


# input 
modelName = "monteSteep"
pars = c("phi", "tau", "gamma", "steep")
expPara = loadExpPara("monte", c("phi", "tau", "gamma", "steep"))
modelFun = getSimModelFun(modelName)

# 
loadExpPara = function(modelName, pars){
  n = 120
  expPara = matrix(NA, n, (length(pars) + 1))
  for(sIdx in 1 : n){
    fileName = sprintf("genData/expModelFitting/%s/s%d.txt", modelName, sIdx)
    junk = read.csv(fileName)
    expPara[sIdx, ] = apply(junk , MARGIN = 2, mean)
  }
  expPara = data.frame(expPara)
  colnames(expPara) = c(pars, "LL_all")
  expPara$id = unique(blockData$id)
  return(expPara)
}

# simluation 
nRep = 10 # number of repetitions
n = length(unique(blockData$id)) # number of subjects 
trialData = vector(length = n * nRep, mode ='list')
repNo = matrix(1 : (n * nRep), nrow = n, ncol = nRep)
for(sIdx in 1 : n){
  id = expPara$id[sIdx]
  para = as.double(expPara[sIdx, 1 : length(pars)])
  cond = unique(blockData$condition[blockData$id == id])
  for(rIdx in 1 : nRep){
    tempt = modelFun(para, cond, 1)
    trialData[[repNo[sIdx, rIdx]]] = tempt
  }
}

# analysis repetion 
plotTrialwiseData = F
plotKMSC = F
plotWTW = F
# initialize 
totalEarnings_ = matrix(0, n, nRep)
AUC_ = matrix(0, n, nRep)
for(sIdx in 1 : n){
  tMax = ifelse(unique(blockData$condition[sIdx]) == "HP", tMaxs[1], tMaxs[2])
  kmGrid = seq(0, tMax, by=0.1) 
  label = "asda"
  for(rIdx in 1 : nRep){
    if (plotTrialwiseData) {
      trialPlots(thisTrialData,label)
    }
    readline("continue")
    thisTrialData = trialData[[repNo[sIdx, rIdx]]]
    totalEarnings_[sIdx, rIdx] =  sum(thisTrialData$trialEarnings)
    kmscResults = kmsc(thisTrialData,tMax,label,plotKMSC,kmGrid)
    AUC_[sIdx, rIdx] = kmscResults[['auc']]
  }
}


blockData$AUCRep = apply(AUC_, MARGIN = 1, FUN = mean)
blockData$totalEarningsRep = apply(totalEarnings_, MARGIN = 1, FUN = mean)
ggplot(blockData, aes(AUC, AUCRep)) + geom_point() + facet_grid(~condition)


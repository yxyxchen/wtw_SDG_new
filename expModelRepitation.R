library("ggplot2")
source("subFxs/plotThemes.R")
load("genData/expDataAnalysis/blockData.RData")
blockData = blockData[blockData$blockNum == 1,]
source("subFxs/taskFxs.R") # used in repetition
source("subFxs/repetitionFxs.R")
source("subFxs/analysisFxs.R") # for analysis
load("wtwSettings.RData") # used in repetition
source("subFxs/loadFxs.R") # 

# input 
modelName = "monte"
pars = c("phi", "tau", "gamma")
expPara = loadExpPara("monte", c("phi", "tau", "gamma"))

# load raw data 
allData = loadAllData()
hdrData = allData$hdrData  
allIDs = hdrData$ID     
expTrialData = allData$trialData       
  

# simluation 
repModelFun = getRepModelFun(modelName)
nRep = 10 # number of repetitions
n = length(unique(blockData$id)) # number of subjects 
trialData = vector(length = n * nRep, mode ='list')
repNo = matrix(1 : (n * nRep), nrow = n, ncol = nRep)
for(sIdx in 1 : n){
  id = allIDs[sIdx]
  para = as.double(expPara[sIdx, 1 : length(pars)])
  cond = unique(blockData$condition[blockData$id == id])
  thisExpTrialData = expTrialData[[id]]
  schedualeWait = thisExpTrialData$scheduledWait[thisExpTrialData$blockNum == 1]
  for(rIdx in 1 : nRep){
    tempt = repModelFun(para, cond, schedualeWait)
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
      readline("continue")
    }
    thisTrialData = trialData[[repNo[sIdx, rIdx]]]
    totalEarnings_[sIdx, rIdx] =  sum(thisTrialData$trialEarnings)
    kmscResults = kmsc(thisTrialData,tMax,label,plotKMSC,kmGrid)
    AUC_[sIdx, rIdx] = kmscResults[['auc']]
  }
}


blockData$AUCRep = apply(AUC_, MARGIN = 1, FUN = mean)
blockData$totalEarningsRep = apply(totalEarnings_, MARGIN = 1, FUN = mean)
ggplot(blockData, aes(AUC, AUCRep)) + geom_point() + facet_grid(~condition)


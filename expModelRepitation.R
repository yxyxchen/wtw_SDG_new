library("ggplot2")
library(stringi)
source("subFxs/plotThemes.R")
load("genData/expDataAnalysis/blockData.RData")
blockData = blockData[blockData$blockNum == 1,]
source("subFxs/taskFxs.R") # used in repetition
source("subFxs/repetitionFxs.R")
source("subFxs/analysisFxs.R") # for analysis
load("wtwSettings.RData") # used in repetition
source("subFxs/loadFxs.R") # 

loadExpPara = function(modelName, pars){
  load("genData/expDataAnalysis/blockData.RData")
  noStressIDList = unique(blockData$id[blockData$stress == "no stress"]) 
  nNoStress = length(noStressIDList)
  nE = (length(pars) + 2) 
  expPara = matrix(NA, nNoStress, nE * 4)
  for(i in 1 : nNoStress){
    ID = noStressIDList[i]
    fileName = sprintf("genData/expModelFitting/%s/s%d_summary.txt", modelName, ID)
    junk = read.csv(fileName, header = F)
    if(ncol(junk) == 11){
      junk = read.csv(fileName, header = T, row.names = 1)
    }
    
    expPara[i, 1:nE] = junk[,1]
    expPara[i, (nE + 1) : (2 * nE)] = junk[,2]
    expPara[i, (2*nE + 1) : (3 * nE)] = junk[,9]
    expPara[i, (3 * nE + 1) : (4 * nE)] = junk[,10]
  }
  expPara = data.frame(expPara)
  
  junk = c(pars, "LL_all", "lp__")
  colnames(expPara) = c(junk, paste0(junk, "SD"), paste0(junk, "Effe"), paste0(junk, "Rhat"))
  expPara$id = noStressIDList  # needed for left_join
  return(expPara)
}
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
noStressIDList = unique(blockData$id[blockData$stress == "no stress"]) 
nNoStress = length(noStressIDList)
trialData = vector(length = nNoStress * nRep, mode ='list')
repNo = matrix(1 : (nNoStress * nRep), nrow = nNoStress, ncol = nRep)
for(sIdx in 1 : nNoStress){
  id = noStressIDList[[sIdx]]
  para = as.double(expPara[sIdx, 1 : length(pars)])
  cond = unique(blockData$condition[blockData$id == id])
  thisExpTrialData = expTrialData[[id]]
  schedualeWait = thisExpTrialData$scheduledWait
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
totalEarnings_ = matrix(0, nNoStress, nRep)
AUC_ = matrix(0, nNoStress, nRep)
timeWaited_ = vector(mode = "list", length = nNoStress)
for(sIdx in 1 : nNoStress){
  id = noStressIDList[[sIdx]]
  tMax = ifelse(unique(blockData$condition[blockData$id == id]) == "HP", tMaxs[1], tMaxs[2])
  kmGrid = seq(0, tMax, by=0.1) 
  label = "asda"
  nTrial = nrow(expTrialData[[id]])
  timeWaitedMatrix = matrix(0, nTrial, nRep)
  for(rIdx in 1 : nRep){
    thisTrialData = trialData[[repNo[sIdx, rIdx]]]
    if (plotTrialwiseData) {
      trialPlots(thisTrialData,label)
      readline("continue")
    }
    junk = thisTrialData$timeWaited
    timeWaitedMatrix[,rIdx] = junk
    totalEarnings_[sIdx, rIdx] =  sum(thisTrialData$trialEarnings)
    kmscResults = kmsc(thisTrialData,tMax,label,plotKMSC,kmGrid)
    AUC_[sIdx, rIdx] = kmscResults[['auc']]
  }
  timeWaited_[[sIdx]] = apply(timeWaitedMatrix, MARGIN = 1, mean)
}

# AUC prediction
plotData = blockData[blockData$stress == "no stress",]
plotData$AUCRep = apply(AUC_, MARGIN = 1, FUN = mean)
useID = noStressIDList[expPara$gammaRhat < 2]
ggplot(plotData[plotData$id %in% useID, ],
       aes(AUC, AUCRep)) + geom_point() + facet_grid(~condition) + geom_abline(slope = 1, xintercept = 0)


# trialData prediction
for(sIdx in 1 : nNoStress){
  thisID = noStressIDList[sIdx]
  nTrial = nrow(expTrialData[[thisID]])
  thisExpTrialData = expTrialData[[thisID]]
  if(thisID %in% useID){
    para = as.double(expPara[sIdx, 1 : length(pars)])
    label = sprintf('Subject %s, %s, %s',thisID, hdrData$condition[sIdx], hdrData$stress[sIdx])
    label = paste(label, paste(round(para, 3), collapse = "", seq = " "))
    
    # prepara data 
    timeWaited = thisExpTrialData$timeWaited
    trialEarnings = thisExpTrialData$trialEarnings
    scheduledWait = thisExpTrialData$scheduledWait
    timeWaited[trialEarnings >0] = scheduledWait[trialEarnings >0]
    plotData = data.frame(trialNum = rep(1 : nTrial, 2), timeWaited = c(timeWaited,
                                                                     timeWaited_[[sIdx]]),
                          quitIdx = rep(trialEarnings == 0, 2), source = rep(c("exp", "rep"), each = nTrial))
    p = ggplot(plotData, aes(trialNum, timeWaited)) + geom_line(alpha = 0.7, aes(color = source))  +
      geom_point(data = plotData[plotData$quitIdx == 1 & plotData$source == "exp", ], aes(trialNum, timeWaited)) + ggtitle(label)+
      displayTheme
      
    print(p)
    readline("continue")
  }
}


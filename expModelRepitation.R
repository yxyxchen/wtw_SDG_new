library("ggplot2")
library(stringi)
source("subFxs/plotThemes.R")
load("genData/expDataAnalysis/blockData.RData")
load("genData/expDataAnalysis/subData.RData")
blockData = blockData[blockData$blockNum == 1,]
source("subFxs/taskFxs.R") # used in repetition
source("subFxs/repetitionFxs.R")
source("subFxs/analysisFxs.R") # for analysis
load("wtwSettings.RData") # used in repetition
source("subFxs/loadFxs.R") # 

# input 
modelName = "monteRP"
pars = c("phi", "tau", "gamma")
pars = c("phiR","phiP", "tau", "gamma")
expPara = loadExpPara(modelName, pars)
useID = noStressIDList[!(expPara$gammaRhat > 1.1 | expPara$phiRhat >1.1 | expPara$tauRhat > 1.1 | 
                           expPara$phiEffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100)]
useID = noStressIDList[!(expPara$gammaRhat > 1.1 | expPara$phiRhat >1.1 | expPara$tauRhat > 1.1 | expPara$steepRhat > 1.1 |
                   expPara$phiEffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100 | expPara$steepEffe < 100)]
useID = noStressIDList[!(expPara$gammaRhat > 1.1 | expPara$phiRRhat >1.1 | expPara$tauRhat > 1.1 | expPara$phiPRhat>1.1 | 
                         expPara$phiREffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100 | expPara$phiPEffe < 100)]
# load raw data 
allData = loadAllData()
hdrData = allData$hdrData  
allIDs = hdrData$ID     
expTrialData = allData$trialData       
  
# simluation 
repModelFun = getRepModelFun(modelName)
nRep = 10# number of repetitions
noStressIDList = unique(blockData$id[blockData$stress == "no stress"]) 
nNoStress = length(noStressIDList)
trialData = vector(length = nNoStress * nRep, mode ='list')
repNo = matrix(1 : (nNoStress * nRep), nrow = nNoStress, ncol = nRep)
for(sIdx in 1 : nNoStress){
  id = noStressIDList[[sIdx]]
  # para = as.double(expPara[sIdx, 1 : length(pars)])
  paraList = read.table(sprintf("genData/expModelFitting/%s/s%d.txt", modelName, id),
                        sep = ",", row.names = NULL)
  cond = unique(blockData$condition[blockData$id == id])
  thisExpTrialData = expTrialData[[id]]
  schedualeWait = thisExpTrialData$scheduledWait
  for(rIdx in 1 : nRep){
    para = as.double(paraList[sample(1 : nrow(paraList), 1), 1 : length(pars)])
    tempt = repModelFun(para, cond, schedualeWait)
    trialData[[repNo[sIdx, rIdx]]] = tempt
  }
}

# calculate AUC 
plotKMSC = F
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
    junk = thisTrialData$timeWaited
    timeWaitedMatrix[,rIdx] = junk
    totalEarnings_[sIdx, rIdx] =  sum(thisTrialData$trialEarnings)
    kmscResults = kmsc(thisTrialData,tMax,label,plotKMSC,kmGrid)
    AUC_[sIdx, rIdx] = kmscResults[['auc']]
  }
  timeWaited_[[sIdx]] = apply(timeWaitedMatrix, MARGIN = 1, mean)
}

# AUC prediction
plotData = subData[subData$stress == "no stress",]
plotData$AUCRep = apply(AUC_, MARGIN = 1, FUN = mean)
plotData$AUCRepSd = apply(AUC_, MARGIN = 1, FUN = sd)
plotData$AUCRepMin = plotData$AUCRep - plotData$AUCRepSd / sqrt(length(useID))
plotData$AUCRepMax = plotData$AUCRep + plotData$AUCRepSd / sqrt(length(useID))
ggplot(plotData[plotData$id %in% useID, ],
       aes(AUC, AUCRep)) + geom_point() + facet_grid(~condition) + 
  geom_abline(slope = 1, intercept = 0) + geom_errorbar(aes(ymin = AUCRepMin, ymax = AUCRepMax)) + 
  displayTheme


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


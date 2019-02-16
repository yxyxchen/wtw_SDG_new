library("ggplot2")
library(stringr)
source("subFxs/plotThemes.R")
load("genData/expDataAnalysis/blockData.RData")
source("subFxs/taskFxs.R") # used in repetition
source("subFxs/repetitionFxs.R")
source("subFxs/analysisFxs.R") # for analysis
load("wtwSettings.RData") # used in repetition
source("subFxs/loadFxs.R") # 
load("wtwSettings.RData")
load("genData/expDataAnalysis/kmOnGrid.RData")

expPara = loadExpPara(modelName, pars)
idList = unique(blockData$id)
n = length(idList)
RhatCols = which(str_detect(colnames(expPara), "hat"))[1 : length(pars)]
EffeCols = which(str_detect(colnames(expPara), "Effe"))[1 : length(pars)]
useID = idList[apply(expPara[,RhatCols] < 1.1, MARGIN = 1, sum) == length(pars) &
                         apply(expPara[,EffeCols] >100, MARGIN = 1, sum) == length(pars)]
# load raw data 
allData = loadAllData()
hdrData = allData$hdrData  
allIDs = hdrData$ID     
expTrialData = allData$trialData       
  
# simluation 
set.seed(123)
repModelFun = getRepModelFun(modelName)
nRep = 10 # number of repetitions
n = length(idList)
trialData = vector(length = n * nRep, mode ='list')
repNo = matrix(1 : (n * nRep), nrow = n, ncol = nRep)
for(sIdx in 1 : n){
  id = idList[[sIdx]]
  # para = as.double(expPara[sIdx, 1 : length(pars)])
  paraList = read.table(sprintf("genData/expModelFitting/%s/s%d.txt", modelName, id),
                        sep = ",", row.names = NULL)
  cond = unique(blockData$condition[blockData$id == id])
  thisExpTrialData = expTrialData[[id]]
  thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
  scheduledWait = thisExpTrialData$scheduledWait
  for(rIdx in 1 : nRep){
    para = as.double(paraList[sample(1 : nrow(paraList), 1), 1 : length(pars)])
    tempt = repModelFun(para, cond, scheduledWait)
    trialData[[repNo[sIdx, rIdx]]] = tempt
  }
}

# calculate AUC and timeWaited
plotKMSC = F
# initialize 
totalEarningsRep_ = matrix(0, n, nRep)
AUCRep_ = matrix(0, n, nRep)
timeWaitedRep_ = vector(mode = "list", length = n)
timeWaitedRepSd_ = vector(mode = "list", length = n)
kmOnGridRep_ = vector(mode = "list", length = n)
kmOnGridRepSd_ = vector(mode = "list", length = n)
for(sIdx in 1 : n){
  id = idList[[sIdx]]
  thisExpTrialData = expTrialData[[id]]
  thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ] 
  tMax = ifelse( unique(thisExpTrialData$condition) == "HP", tMaxs[1], tMaxs[2])
  kmGrid = seq(0, tMax, by=0.1) 
  label = "asda"
  nTrial = nrow(thisExpTrialData)
  # initialize
  timeWaitedMatrix = matrix(0, nTrial, nRep)
  kmOnGridMatrix = matrix(0, length(kmGrid), nRep)
  for(rIdx in 1 : nRep){
    thisTrialData = trialData[[repNo[sIdx, rIdx]]]
    junk = thisTrialData$timeWaited
    timeWaitedMatrix[,rIdx] = junk
    
    totalEarningsRep_[sIdx, rIdx] =  sum(thisTrialData$trialEarnings)
    kmscResults = kmsc(thisTrialData,tMax,label,plotKMSC,kmGrid)
    AUCRep_[sIdx, rIdx] = kmscResults[['auc']]
    kmOnGridMatrix[,rIdx] = kmscResults$kmOnGrid
  }
  timeWaitedRep_[[sIdx]] = apply(timeWaitedMatrix, MARGIN = 1, mean)
  timeWaitedRepSd_[[sIdx]] = apply(timeWaitedMatrix, MARGIN = 1, sd)
  kmOnGridRep_[[sIdx]] = apply(kmOnGridMatrix, MARGIN = 1, mean)
  kmOnGridRepSd_[[sIdx]] = apply(kmOnGridMatrix, MARGIN = 1, sd)
}

# AUC prediction
plotData = blockData[blockData$blockNum == 1,]
plotData$AUCRep = apply(AUCRep_, MARGIN = 1, FUN = mean)
plotData$AUCRepSd = apply(AUCRep_, MARGIN = 1, FUN = sd)
plotData$AUCRepMin = plotData$AUCRep - plotData$AUCRepSd / sqrt(length(useID))
plotData$AUCRepMax = plotData$AUCRep + plotData$AUCRepSd / sqrt(length(useID))
ggplot(plotData[plotData$id %in% useID, ],
       aes(AUC, AUCRep)) + geom_point() + facet_grid(~condition) + 
  geom_abline(slope = 1, intercept = 0) + geom_errorbar(aes(ymin = AUCRepMin, ymax = AUCRepMax)) + 
  saveTheme
fileName = sprintf("figures/expModelRepitation/AUC_AUCRep_%s.pdf", modelName)
ggsave(filename = fileName,  width = 6, height = 4)


# survival curve prediction
for(sIdx in 1 : n){
  thisID = idList[sIdx]
  if(thisID %in% useID){
    para = as.double(expPara[sIdx, 1 : length(pars)])
    label = sprintf('Subject %s, %s, %s, LL = %.1f',thisID, hdrData$condition[sIdx], hdrData$stress[sIdx], expPara$LL_all[sIdx])
    label = paste(label, paste(round(para, 3), collapse = "", seq = " "))
    # prepara data 
    cond = subData$condition[[which(subData$id == thisID)]]
    kmOnGrid = kmOnGrid_[[which(subData$id == thisID)]]
    tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
    kmGrid = seq(0, tMax, by = 0.1)
    kmOnGrid = kmOnGrid_[[which(subData$id == thisID)]]
    kmOnGridRep = kmOnGridRep_[[which(idList== thisID)]]
    junk = data.frame(time = kmGrid, exp = kmOnGrid, rep= kmOnGridRep)
    plotData = gather(junk, source, survival_rate, -time)
    p = ggplot(plotData, aes(time, survival_rate, color = source)) + geom_line() + ggtitle(label) + displayTheme
    print(p)
    readline("continue")
  }
}


# trialData prediction
for(sIdx in 1 : n){
  thisID = idList[sIdx]
  thisExpTrialData = expTrialData[[thisID]]
  thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum == 1,]
  nTrial = nrow(thisExpTrialData)
  if(thisID %in% useID){
    para = as.double(expPara[expPara$id == thisID, 1 : length(pars)])
   
    # prepara data 
    timeWaited = thisExpTrialData$timeWaited
    trialEarnings = thisExpTrialData$trialEarnings
    scheduledWait = thisExpTrialData$scheduledWait
    timeWaited[trialEarnings >0] = scheduledWait[trialEarnings >0]
    nAction = sum(round(ifelse(trialEarnings >0, ceiling(timeWaited / stepDuration), floor(timeWaited / stepDuration) + 1)))
    label = sprintf('Subject %s, %s, %s, -LL = %.2f',thisID, unique(blockData$condition[blockData$id == thisID]),
                    unique(blockData$stress[blockData$id == thisID]), -expPara$LL_all[expPara$id == thisID] / nAction)
    label = paste(label, paste(round(para, 3), collapse = "", seq = " "))
    
    # prepare plotData
    plotData = data.frame(trialNum = rep(1 : nTrial, 2), timeWaited = c(timeWaited,
                                                                     timeWaitedRep_[[which(expPara$id == thisID)]]),
                          quitIdx = rep(trialEarnings == 0, 2), source = rep(c("exp", "rep"), each = nTrial))
    p = ggplot(plotData, aes(trialNum, timeWaited)) + geom_line(aes(color = source, alpha = source))  +
      geom_point(data = plotData[plotData$quitIdx == 1 & plotData$source == "exp", ], aes(trialNum, timeWaited)) + ggtitle(label)+
      displayTheme + scale_alpha_manual(values = c(0.8, 0.5))
    print(p)
    readline("continue")
  }
}

# case study
thisID = 21
cond = unique(blockData$condition[blockData$id ==thisID])
# extract expData
thisExpTrialData = expTrialData[[thisID]]
thisExpTrialData = thisExpTrialData[thisExpTrialData$blockNum ==1, ]
timeWaited = thisExpTrialData$timeWaited
trialEarnings = thisExpTrialData$trialEarnings
scheduledWait = thisExpTrialData$scheduledWait
timeWaited[trialEarnings >0] = scheduledWait[trialEarnings >0]
nTrial = nrow(thisExpTrialData)

## simulate 
para = c(0.1, 8, 0.95, 2, 1)
tempt = monteWini(para, cond, scheduledWait)

## transfer temptdata
# plot
label = sprintf('Subject %s, %s, %s',thisID, unique(blockData$condition[blockData$id == thisID]),
                unique(blockData$stress[blockData$id == thisID]))
label = paste(label, paste(round(para, 3), collapse = "", seq = " "))
plotData = data.frame(trialNum = rep(1 : nTrial, 2), timeWaited = c(timeWaited,
                                                                    tempt$timeWaited),
                      quitIdx = rep(trialEarnings == 0, 2), source = rep(c("exp", "rep"), each = nTrial))

ggplot(plotData, aes(trialNum, timeWaited)) + geom_line(aes(color = source, alpha = source))  +
  geom_point(data = plotData[plotData$quitIdx == 1 & plotData$source == "exp", ], aes(trialNum, timeWaited)) + ggtitle(label)+
  displayTheme +scale_alpha_manual(values = c(0.8, 0.5))


actionValueViewer(tempt$vaWaits, tempt$vaQuits, tempt)

apply(tempt$vaWaits, MARGIN = 2, which.min)
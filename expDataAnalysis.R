# load libraries
source('subFxs/loadFxs.R') # for loading data 
source('subFxs/analysisFxs.R') # for analysis 
source("subFxs/plotThemes.R")
library("ggplot2")
library('dplyr')
# library(Hmisc)

# load setting parameters 
load("wtwSettings.RData")

# load all data
allData = loadAllData()
hdrData = allData$hdrData           
trialData = allData$trialData       
# list with a named element for each subject ID.
allIDs = hdrData$ID                   # column of subject IDs
n = length(allIDs)                    # n
nBlock = 3
cat('Analyzing data for n','=',n,'subjects.\n')

# control which individual-level plots to generate
plotTrialwiseData = T
plotKMSC = F
plotWTW = F
plotTimeEarnings = F   # no good effect 
plotTrialEarnings =  F  # no good effect 

# initialize outputs, organised by block
tGrid = seq(0, blockSecs, by = 0.1)
AUC = numeric(length =n * nBlock)
totalEarnings =  numeric(length =n * nBlock)
wtwEarly = numeric(length =n * nBlock)
# descriptive statistics for individual subjects and blocks
for (sIdx in 1 : n) {
  thisID = allIDs[sIdx]
  for (bkIdx in 1:nBlock){
    # select data 
    thisTrialData = trialData[[thisID]]
    thisBlockIdx = (thisTrialData$blockNum == bkIdx)
    thisTrialData = thisTrialData[thisBlockIdx,]
      
    # generate arguments for later analysis 
    label = sprintf('Subject %s, Cond %s, Stress %s)',thisID, hdrData$condition[sIdx], hdrData$stress[sIdx])
    noIdx = (sIdx - 1) * nBlock + bkIdx # 
    tMax = ifelse(hdrData$condition[sIdx] == conditions[1], tMaxs[1], tMaxs[2])
    kmGrid = seq(0, tMax, by=0.1) # grid on which to average survival curves.
    
    # calcualte totalEarnings
    totalEarnings[noIdx] =  sum(thisTrialData$trialEarnings)
    
    # plot trial-by-trial data
    if (plotTrialwiseData) {
      trialPlots(thisTrialData,label)
    }
    
    # survival analysis
    kmscResults = kmsc(thisTrialData,tMax,label,plotKMSC,kmGrid)
    AUC[noIdx] = kmscResults[['auc']]


    # WTW time series
    wtwCeiling = tMax
    wtwtsResults = wtwTS(thisTrialData, tGrid, wtwCeiling, label, plotWTW)
    wtwEarly[noIdx] = mean(wtwtsResults[1 : (1 * 60 * 10)])
    # accumutative timeEarnings 
    timeEarnings = getTimeEarnings(thisTrialData, tGrid, label, plotTimeEarnings)
    
    # accumutative trialEarnings 
    plotData = data.frame(trialNum = thisTrialData$trialNum,
                          cumEarnings = cumsum(thisTrialData$trialEarnings))
    if(plotTrialEarnings){
      p = ggplot(plotData, aes(trialNum, cumEarnings)) + geom_line()   
      print(p)
    }
    
    # wait for input before continuing, if individual plots were requested
    if (any(plotTrialwiseData, plotKMSC, plotWTW, plotTimeEarnings, plotTrialEarnings)) {
      readline(prompt = paste('subject',thisID, "block", bkIdx, '(hit ENTER to continue)'))
      graphics.off()
    }
  } # loop over blocks
}


# organize and save blockWiseData
blockData = data.frame(id = rep(allIDs, each = nBlock), blockNum = rep( t(1 : nBlock), n),
                       cbal = rep(hdrData$cbal, each = nBlock), condition = factor(rep(hdrData$condition, each = nBlock), levels = c("HP", "LP")),
                       stress = factor(rep(hdrData$stress, each = nBlock), levels = c("no stress", "stress")), AUC = AUC, wtwEarly = wtwEarly,
                       totalEarnings = totalEarnings)
save(blockData, file = 'genData/expDataAnalysis/blockData.RData')


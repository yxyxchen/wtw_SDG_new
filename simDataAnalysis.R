# this script analysized the simulation data on the group level
modelName = "monteCfd"
load(sprintf("genData/simulation/%s/realParas.RData", modelName)) # for nComb, paraNames, nPara
############ load data and functions #########
# generally loading 
library("ggplot2")
library("dplyr")
library("tidyr")
source(file = './subFxs/plotThemes.R')
load("wtwSettings.RData")

# specific loading 
library('scales')
source('subFxs/analysisFxs.R')
load(sprintf("genData/simulation/%s/trialHPData.RData", modelName)) # trialHPData
load(sprintf("genData/simulation/%s/trialLPData.RData", modelName)) # trialLPData
nRep = 5

# trial-level analysis 
plotTrialwiseData = F
plotKMSC = F
plotWTW = F
# loop over conditions 
for(cIdx in 1 : 2){
  cond = conditions[cIdx]
  if(cIdx == 1){
    trialData = trialHPData
  }else{
    trialData = trialLPData
  }
  # generate arguments for later analysis 
  tMax = tMaxs[cIdx]
  kmGrid = seq(0, tMax, by=0.1) 
  
  # initialize outputs
  totalEarnings_ = matrix(0, nComb, nRep)
  AUC_ = matrix(0, nComb, nRep)
  
  # loop over combIdx
  for(combIdx in 1 : nComb){
    # loop over repetitions 
    for(rIdx in 1 : nRep){
      # select data
      thisTrialData = trialData[[simIdx[combIdx, rIdx]]]
      realPara = realParas[combIdx,]
      
      # generate arguments for later analysis 
      label = sprintf('%s %s', conditions[cIdx], paste0(paraNames, realPara, collapse = " "))
      
      # summarise totalEarnings 
      totalEarnings_[combIdx, rIdx] =  sum(thisTrialData$trialEarnings)
      
      # plot trial-by-trial data
      if (plotTrialwiseData) {
        trialPlots(thisTrialData,label)
      }
      
      # survival analysis
      kmscResults = kmsc(thisTrialData,tMax,label,plotKMSC,kmGrid)
      AUC_[combIdx, rIdx] = kmscResults[['auc']]
      
      # WTW time series
      wtwCeiling = tMax
      wtwtsResults = wtwTS(thisTrialData, tGrid, wtwCeiling, label, plotWTW)
    }# end of repetitions
  }# end of para combinations 
  # organize and save blockData
  AUC = rowSums(AUC_) / nRep
  totalEarnings = rowSums(totalEarnings_) / nRep
  blockData = data.frame(id = 1 : combIdx, blockNum = rep(1, nComb),
                         condition = factor(rep(cond, each = nComb), levels = c("HP", "LP")),
                         AUC = AUC, 
                         totalEarnings = totalEarnings)
  if(cIdx == 1){
    blockHPData = blockData
  }else{
    blockLPData = blockData
  }
}

####### plot distribution of totalEarnings
blockData = data.frame(rbind(blockHPData, blockLPData), rbind(realParas))
colnames(blockData)[(ncol(blockData) - nPara + 1) : ncol(blockData) ] = paraNames
ggplot(blockData, aes(totalEarnings)) + geom_histogram(bins = 15) +
  facet_wrap(~condition, nrow = 1) + xlab('Total earnings') + ylab("Num of simulations") + saveTheme + xlim(c(0, 600))
fileName = sprintf("figures/simDataAnalysis/%s/totalEarnings.pdf", modelName)
ggsave(fileName, width = 16, height = 8)

# calculate range
dplyr::summarise(group_by(blockData, condition),
          minEarning = min(totalEarnings),
          maxEarning = max(totalEarnings))

############ summarise para effects on total earnings ###########
paraValues = 1:5
paraData = data.frame(condition = rep(c("HP", "LP"), each = nValue, nPara),
                         paraNames = rep(paraNames, each = nValue * 2),
                         paraValues = rep(paraValues, nPara * 2))
paraData$paraNames = factor(paraData$paraNames, levels = paraNames)


# summarise mu 
muByPhi = summarise_at(group_by(blockData, condition, phi), vars(AUC:totalEarnings), mean)
muByTau = summarise_at(group_by(blockData, condition, c), vars(AUC:totalEarnings), mean)
muByGamma = summarise_at(group_by(blockData, condition, gamma), vars(AUC:totalEarnings), mean)

# summarise sd
stdByPhi = summarise_at(group_by(blockData, condition, phi), vars(AUC:totalEarnings), sd)
stdByTau = summarise_at(group_by(blockData, condition, tau), vars(AUC:totalEarnings), sd)
stdByGamma = summarise_at(group_by(blockData, condition, gamma), vars(AUC:totalEarnings), sd)

# 
mu = rbind(muByPhi, muByTau, muByGamma);
mu = mu[, 3:4]
std = rbind(stdByPhi, stdByTau, stdByGamma)
std = std[,3:4]
max= mu + std
min = mu -std

summaryEarnData = cbind(paraData, mu[,2], std[,2], max[,2], min[,2]);
summaryAUCData = cbind(paraData, mu[,1], std[,1], max[,1], min[,1]);
colnames(summaryEarnData) = c(colnames(paraData), 'mu', 'std', 'max', 'min')
colnames(summaryAUCData) = c(colnames(paraData), 'mu', 'std', 'max', 'min')

# plot 
for(c in 1:2){
  cond = conditions[c]
  ggplot(summaryAUCData[summaryAUCData$condition == cond,], aes(factor(paraValues), mu)) +
    geom_bar(stat = "identity", width=0.5, fill = conditionColors[c]) + geom_errorbar(aes(ymin = min, ymax = max), width=.2)+
    facet_wrap(~paraNames, nrow = 1)+ saveTheme +
    xlab("Parameter value") + ylab("AUC / s") + ggtitle(cond) 
  fileName = sprintf("figures/simDataAnalysis/%s/paraAUCEffect%s.pdf", modelName, cond)
  ggsave(fileName, width = 16, height = 8) 
}

# plot 
for(c in 1:2){
  cond = conditions[c]
  ggplot(summaryEarnData[summaryEarnData$condition == cond,], aes(factor(paraValues), mu)) +
    geom_bar(stat = "identity", width=0.5, fill = conditionColors[c]) + geom_errorbar(aes(ymin = min, ymax = max), width=.2)+
    facet_wrap(~paraNames, nrow = 1)+ saveTheme +
    xlab("Parameter value") + ylab("Total Earnings") + ggtitle(cond) 
  fileName = sprintf("figures/simDataAnalysis/%s/paraEarnEffect%s.pdf",modelName, cond)
  ggsave(fileName, width = 16, height = 8) 
}

######### plot AUC against totalEarnings #######
# prepare data
plotData = blockData %>% arrange(totalEarnings) %>%group_by(condition) %>%
  mutate(earningRank = rank(totalEarnings, ties.method = "first"))

# plot for LP
ggplot(plotData[plotData$condition == 'LP',], aes(AUC, totalEarnings)) + geom_point(size = 1.5) +
  saveTheme + ylab('Total earnings') + xlim(c(0, tMaxs[2])) + ylim(c(0, 500))
fileName = sprintf("figures/simDataAnalysis/%s/AUCLP_earningsLP.pdf", modelName)
ggsave(fileName, width = 6, height = 4)

# plot for HP
ggplot(plotData[plotData$condition == 'HP',], aes(AUC, totalEarnings)) + geom_point(size = 1.5) +
  saveTheme + ylab('Total earnings') + xlim(c(0, tMaxs[1])) + ylim(c(0, 500))
fileName = sprintf("figures/simDataAnalysis/%s/AUCLP_earningsHP.pdf", modelName)
ggsave(fileName, width = 6, height = 4)


######## plot the timeseries of wtw #######
# meanValues = c(apply(rawWTW$HP, MARGIN = 3, FUN = mean), 
#                apply(rawWTW$LP, MARGIN = 3, FUN = mean))
# stdValues = c(apply(rawWTW$HP, MARGIN = 3, FUN = sd), 
#                apply(rawWTW$LP, MARGIN = 3, FUN = sd))
# plotData = data.frame(meanValues, stdValues,
#                       time = rep(tGrid, time = 2),
#                       condition = rep(c('HP', 'LP'), each = length(tGrid)),
#                       minValues = meanValues - stdValues / sqrt(dim(rawWTW$HP)[1]),
#                       maxValues = meanValues + stdValues / sqrt(dim(rawWTW$HP)[1]))
# 
# ggplot(plotData, aes(time, meanValues, color = condition)) + 
#   geom_ribbon(data = plotData[plotData$condition == 'HP',], aes(ymin=minValues, ymax=maxValues),linetype=0, alpha = 0.1, color = "#bababa") +
#   geom_ribbon(data = plotData[plotData$condition == 'LP',], aes(ymin=minValues, ymax=maxValues),linetype=0, alpha = 0.1, color = "#bababa") + 
#   geom_line(size = 1) + xlab('Time in block / s') + ylab('WTW / s') + saveTheme 
# fileName = file.path(outFile, "wtwTimeSeries.pdf")
# ggsave(fileName, width = 12, height = 8)

########### plot HPAUC against LPAUC ##############
plotData = data.frame(HPAUC = blockHPData$AUC, LPAUC = blockLPData$AUC)
ggplot(plotData, aes(HPAUC, LPAUC)) + geom_point(shape = 3 ) + geom_smooth(method = lm) +
  xlab('HP AUC/s') + ylab("LP AUC /s") + saveTheme
fileName = fileName = sprintf("figures/simDataAnalysis/%s/HPAUC_LPAUC.pdf", modelName)
ggsave(fileName, width = 8, height = 8)



















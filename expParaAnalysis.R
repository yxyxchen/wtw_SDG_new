library("ggplot2")
source("subFxs/plotThemes.R")
source("subFxs/loadFxs.R")
# load experimental data
load("genData/expDataAnalysis/blockData.RData")
blockData = blockData[blockData$blockNum == 1,]
idList = unique(blockData$id) 
n = length(idList)
plotParaAUC = function(expPara, paraName, blockData, useID){
  expPara = expPara[(expPara$id %in% useID),]
  blockData = blockData[(blockData$id %in% useID),]
  plotData = cbind(expPara, blockData)
  corHP = cor.test(plotData[plotData$condition == "HP", paraName], plotData[plotData$condition == "HP",]$AUC, method = "spearman")
  corLP = cor.test(plotData[plotData$condition == "LP", paraName], plotData[plotData$condition == "LP",]$AUC, method = "spearman")
  rhoHP = round(corHP$estimate, 3)
  rhoLP= round(corLP$estimate, 3)
  pHP = round(corHP$p.value, 3)
  pLP = round(corLP$p.value, 3)

  textData = data.frame(label = c(paste(rhoHP, "(p =", pHP, ")"), paste(rhoLP, "(p =", pLP, ")")),
                        condition = c("HP", "LP"))
  plotData = left_join(expPara, blockData, by = "id")
  p = ggplot(plotData, aes_string(paraName, "AUC")) + geom_point() + facet_grid(~condition)  + geom_text(data = textData,
    aes(x = -Inf,y = -Inf, 
    label = label),
    hjust   = -0.5,
    vjust   = -3,
    color = "blue",
    size = 5
    ) + displayTheme
  print(p)
} 

####### monte
ggplot(expPara, aes(QwaitIni)) + geom_histogram() + facet_grid(~condition)
expPara = loadExpPara("cons_partFlex_gamma", c("phi", "tau", "gamma", "QwaitIni"))
useID = idList[!(expPara$gammaRhat > 1.1 | expPara$phiRhat >1.1 | expPara$tauRhat > 1.1 | 
                           expPara$phiEffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100)]
plotParaAUC(expPara, "QwaitIni", blockData, useID)
plotParaAUC(expPara, "tau", blockData, useID)
plotParaAUC(expPara, "gamma", blockData, useID)


####### monteRP
expPara = loadExpPara("monteRP", c("phiR", "phiP", "tau", "gamma"))
useID = idList[!(expPara$gammaRhat > 1.1 | expPara$phiRRhat >1.1 | expPara$tauRhat > 1.1 | expPara$phiPRhat > 1.1| 
                           expPara$phiREffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100 | expPara$phiPEffe < 100)]
plotParaAUC(expPara, "phiR", blockData, useID)
plotParaAUC(expPara, "phiP", blockData, useID)
plotParaAUC(expPara, "tau", blockData, useID)
plotParaAUC(expPara, "gamma", blockData, useID)
plotParaAUC(expPara, "optimism", blockData, useID)
ggplot(expPara, aes(phiR, phiP)) + geom_point() + geom_abline(slope = 1, intercept = 0)
expPara$optimism = expPara$phiR - expPara$phiP


#####
#monteSteep
expPara = loadExpPara("monteSteep", c("phi", "tau", "gamma", "steep")) # numPara and nPara
useID = idList[!(expPara$gammaRhat > 1.1 | expPara$phiRhat >1.1 | expPara$tauRhat > 1.1 | expPara$steepRhat > 1.1 |
                           expPara$phiEffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100 | expPara$steepEffe < 100)]
plotParaAUC(expPara, "tau", blockData, useID)
plotParaAUC(expPara, "steep", blockData, useID)
plotParaAUC(expPara, "phi", blockData, useID)
plotParaAUC(expPara, "gamma", blockData, useID)
hist(expPara$steep)


#####
#monteRatio
expPara = loadExpPara("monteRatio", c("phi", "tau", "gamma", "ratio")) # numPara and nPara
RhatCols = which(str_detect(colnames(expPara), "hat"))[1 : length(pars)]
EffeCols = which(str_detect(colnames(expPara), "Effe"))[1 : length(pars)]
useID = idList[apply(expPara[,RhatCols] < 1.1, MARGIN = 1, sum) == length(pars) &
                         apply(expPara[,EffeCols] >100, MARGIN = 1, sum) == length(pars)]
plotParaAUC(expPara, "tau", blockData, useID)
plotParaAUC(expPara, "ratio", blockData, useID)
plotParaAUC(expPara, "phi", blockData, useID)
plotParaAUC(expPara, "gamma", blockData, useID) # log gamma
hist(expPara$ratio)

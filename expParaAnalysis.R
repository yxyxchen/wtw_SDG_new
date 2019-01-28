library("ggplot2")
source("subFxs/plotThemes.R")
load("genData/expDataAnalysis/blockData.RData")

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

plotParaAUC = function(expPara, paraName, blockData){
  rhoHP = vector(length = 3)
  rhoLP = vector(length = 3)
  pHP = vector(length = 3)
  pLP = vector(length = 3)
  for(i in 1 : 3){
    plotData = cbind(expPara, blockData[blockData$blockNum == i, ])
    corHP = cor.test(plotData[plotData$condition == "HP", paraName], plotData[plotData$condition == "HP",]$AUC, method = "spearman")
    corLP = cor.test(plotData[plotData$condition == "LP", paraName], plotData[plotData$condition == "LP",]$AUC, method = "spearman")
    rhoHP[i] = round(corHP$estimate, 3)
    rhoLP[i] = round(corLP$estimate, 3)
    pHP[i] = round(corHP$p.value, 3)
    pLP[i] = round(corLP$p.value, 3)
  }
  textData = data.frame(label = c(paste0(rhoHP, "(p =", pHP, ")"), paste0(rhoLP, "(p =", pLP, ")")),
                        condition = rep(c("HP", "LP"), each = 3),
                        blockNum = rep(1:3, 2))
  plotData = left_join(expPara, blockData, by = "id")
  p = ggplot(plotData, aes_string(paraName, "AUC")) + geom_point() + facet_grid(condition ~ blockNum)  + geom_text(
      data    = textData,
      mapping = aes(x = -Inf, y = -Inf, label = label),
      hjust   = -0.5,
      vjust   = -3,
      color = "blue",
      size = 5
    )
  print(p)
} 

####### monte
expPara = loadExpPara("monte")
plotParaAUC(expPara, "phi", blockData)
plotParaAUC(expPara, "tau", blockData)
plotParaAUC(expPara, "gamma", blockData)
#####

expPara = loadExpPara("monteRP", c("phiR", "phiP", "tau", "gamma")) # numPara and nPara
plotParaAUC(expPara, "tau", blockData)
plotParaAUC(expPara, "gamma", blockData)
expPara$optimism = expPara$phiR - expPara$phiP
plotParaAUC(expPara, "optimism", blockData)
plotParaAUC(expPara, "phiR", blockData)
plotParaAUC(expPara, "phiP", blockData)
##################
# monteBias
expPara = loadExpPara("monteBias", c("phi", "tau", "gamma", "quitBias")) # numPara and nPara
plotParaAUC(expPara, "tau", blockData)
plotParaAUC(expPara, "gamma", blockData)
plotParaAUC(expPara, "phi", blockData)
plotParaAUC(expPara, "quitBias", blockData)

###################
#monteSteep
expPara = loadExpPara("monteSteep", c("phi", "tau", "gamma", "steep")) # numPara and nPara
plotParaAUC(expPara, "tau", blockData)
plotParaAUC(expPara, "gamma", blockData)
plotParaAUC(expPara, "phi", blockData)
plotParaAUC(expPara, "steep", blockData)





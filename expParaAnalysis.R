library("ggplot2")
source("subFxs/plotThemes.R")
load("genData/expDataAnalysis/blockData.RData")

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

ggplot(expPara, aes(phi)) + geom_histogram()
ggplot(expPara, aes(tau)) + geom_histogram()
ggplot(expPara, aes(gamma)) + geom_histogram()
ggplot(expPara, aes(gammaSD)) + geom_histogram()
ggplot(expPara, aes(gammaRhat)) + geom_histogram() +
  geom_vline(xintercept = 1, color = "red")
ggplot(expPara, aes(gammaRhat)) + geom_histogram() +
  geom_vline(xintercept = 1, color = "red")

ggplot(expPara, aes(gammaEffe)) + geom_histogram()


plotParaAUC = function(expPara, paraName, blockData, deleteID){
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
      size = 2
    )
  print(p)
} 

####### monte

expPara = loadExpPara("monte", c("phi", "tau", "gamma"))
useID = noStressIDList[expPara$gammaRhat < 2]
plotParaAUC(expPara[expPara$id %in% useID, c(1:nE, 21)], "phi",
            blockData[blockData$stress == "no stress" & blockData$id %in% useID, ])
plotParaAUC(expPara[expPara$id %in% useID, c(1:nE, 21)], "gamma",
            blockData[blockData$stress == "no stress" & blockData$id %in% useID, ])
plotParaAUC(expPara[expPara$id %in% useID, c(1:nE, 21)], "tau",
            blockData[blockData$stress == "no stress" & blockData$id %in% useID, ])
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





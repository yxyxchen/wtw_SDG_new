library("ggplot2")
source("subFxs/plotThemes.R")
load("genData/expDataAnalysis/subData.RData")


ggplot(expPara, aes(phi)) + geom_histogram()
ggplot(expPara, aes(tau)) + geom_histogram()
ggplot(expPara, aes(gamma)) + geom_histogram()
ggplot(expPara, aes(gammaSD)) + geom_histogram()
ggplot(expPara, aes(gammaRhat)) + geom_histogram() +
  geom_vline(xintercept = 1, color = "red")
ggplot(expPara, aes(gammaRhat)) + geom_histogram() +
  geom_vline(xintercept = 1, color = "red")
ggplot(expPara, aes(gammaEffe)) + geom_histogram()


plotParaAUC = function(expPara, paraName, subData, useID){
  expPara = expPara[(expPara$id %in% useID),]
  subData = subData[(subData$id %in% useID),]
  plotData = cbind(expPara, subData)
  corHP = cor.test(plotData[plotData$condition == "HP", paraName], plotData[plotData$condition == "HP",]$AUC, method = "spearman")
  corLP = cor.test(plotData[plotData$condition == "LP", paraName], plotData[plotData$condition == "LP",]$AUC, method = "spearman")
  rhoHP = round(corHP$estimate, 3)
  rhoLP= round(corLP$estimate, 3)
  pHP = round(corHP$p.value, 3)
  pLP = round(corLP$p.value, 3)

  textData = data.frame(label = c(paste(rhoHP, "(p =", pHP, ")"), paste(rhoLP, "(p =", pLP, ")")),
                        condition = c("HP", "LP"))
  plotData = left_join(expPara, subData, by = "id")
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
expPara = loadExpPara("monte", c("phi", "tau", "gamma"))
useID = noStressIDList[!(expPara$gammaRhat > 1.1 | expPara$phiRhat >1.1 | expPara$tauRhat > 1.1 | 
                           expPara$phiEffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100)]
plotParaAUC(expPara, "phi", subData, useID)
plotParaAUC(expPara, "tau", subData, useID)
plotParaAUC(expPara, "gamma", subData, useID)


####### monteRP
expPara = loadExpPara("monteRP", c("phiR", "phiP", "tau", "gamma"))
useID = noStressIDList[!(expPara$gammaRhat > 1.1 | expPara$phiRRhat >1.1 | expPara$tauRhat > 1.1 | expPara$phiPRhat > 1.1| 
                           expPara$phiREffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100 | expPara$phiPEffe < 100)]
plotParaAUC(expPara, "phiR", subData, useID)
plotParaAUC(expPara, "phiP", subData, useID)
plotParaAUC(expPara, "phiR", subData, useID)
plotParaAUC(expPara, "tau", subData, useID)
plotParaAUC(expPara, "gamma", subData, useID)
plotParaAUC(expPara, "optimism", subData, useID)
ggplot(expPara, aes(phiR, phiP)) + geom_point() + geom_abline(slope = 1, intercept = 0)
expPara$optimism = expPara$phiR - expPara$phiP


#####
#monteSteep
expPara = loadExpPara("monteSteep", c("phi", "tau", "gamma", "steep")) # numPara and nPara
useID = noStressIDList[!(expPara$gammaRhat > 1.1 | expPara$phiRhat >1.1 | expPara$tauRhat > 1.1 | expPara$steepRhat > 1.1 |
                           expPara$phiEffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100 | expPara$steepEffe < 100)]
plotParaAUC(expPara, "tau", subData, useID)
plotParaAUC(expPara, "steep", subData, useID)
plotParaAUC(expPara, "phi", subData, useID)
plotParaAUC(expPara, "gamma", subData, useID)
hist(expPara$steep)

n = 120
modelName = "monteRP"
samplePara = vector("list", length = n)
<<<<<<< HEAD
summaryPara = data.frame("phi" = vector(length = n), "phiP" = vector(length = n),
=======
summaryPara = data.frame("phiR" = vector(length = n), "phiP" = vector(length = n),
>>>>>>> 76dded834aafefccba54ad4be637489b7f59cbba
                         "tau" = vector(length = n), "gamma" = vector(length = n),
                         "LL_all" = vector(length = n))
for(sIdx in 1 : n){
  fileName = sprintf("genData/expModelFitting/%s/s%d.txt", modelName, sIdx)
  samplePara[[sIdx]] = read.csv(fileName)
  summaryPara[sIdx, ] = apply(samplePara[[sIdx]] , MARGIN = 2, mean)
}
summaryPara$optimism = summaryPara$phiR - summaryPara$phiP
  
library(ggplot2)
ggplot(summaryPara, aes(phiR, phiP)) + geom_point() + xlim(c(0,0.2)) +
  ylim(c(0, 0.2)) + geom_abline(slope = 1, yintercept = 0)

hist(summaryPara$optimism)

### relation ship between AUC and 
plotData = cbind(summaryPara, blockData[blockData$blockNum == 3, ])
ggplot(plotData, aes(optimism, AUC)) + geom_point() + facet_grid(~condition)

cor.test(plotData$optimism, plotData$AUC)

<<<<<<< HEAD

##################


n = 120
modelName = "monteBias"
samplePara = vector("list", length = n)
summaryPara = data.frame("phi" = vector(length = n), 
                         "tau" = vector(length = n), "gamma" = vector(length = n), "bias" = vector(length = n),
                         "LL_all" = vector(length = n))
for(sIdx in 1 : n){
  fileName = sprintf("genData/expModelFitting/%s/s%d.txt", modelName, sIdx)
  samplePara[[sIdx]] = read.csv(fileName)
  summaryPara[sIdx, ] = apply(samplePara[[sIdx]] , MARGIN = 2, mean)
}

hist(summaryPara$bias)

### relation ship between AUC and 
library("ggplot2")
load("genData/expDataAnalysis/blockData.RData")
plotData = cbind(summaryPara, blockData[blockData$blockNum == 1, ])
ggplot(plotData, aes(bias, AUC)) + geom_point() + facet_grid(~condition)

cor.test(plotData$optimism, plotData$AUC)
=======
>>>>>>> 76dded834aafefccba54ad4be637489b7f59cbba
# this script extract waic for all models. 
library("R.matlab")
library("stringr")
library("ggplot2")
n = 120
modelNames = unlist(lapply(list.files(path="stanModels", pattern="*.stan"),
                    FUN = function(x) str_sub(x, 1, -1 - 5)))
modelNames = factor(modelNames,
                   levels = c("baseline", "monteNoDiscount", "monte",
                              "monteBias", "monteRP", "monteRatio",
                              "monteSteep", "monteRPBias", "monteRPBiasSteep"),
                   ordered = T)
nModel = length(modelNames)


logEvidenceList = matrix(NA, n, nModel)
pWaicList = matrix(NA, n, nModel)
for(m in 1 : nModel){
  modelName = modelNames[m]
  logEvidence = vector(length = n)
  pWaic = vector(length = n)
  for(sIdx in 1 : n){
    fileName = sprintf("genData/expModelFitting/%s/s%d.RData", modelName, sIdx)
    load(fileName)
    logEvidence[sIdx] = WAIC$elpd_waic
    pWaic[sIdx] = WAIC$p_waic
  }
  logEvidenceList[,m] = logEvidence
  pWaicList[,m] = pWaic
}
waicList = -2 * logEvidenceList
f= "genData/expModelFitting/logEvidenceList.mat"
writeMat(f, logEvidenceList = logEvidenceList, modelNames = modelNames)

### analysis 
load("genData/expDataAnalysis/blockData.RData")
waicList[,which(modelNames == "baseline")] - waicList[,which(modelNames == "monte")] 
hist(waicList[,which(modelNames == "baseline")] - waicList[,which(modelNames == "monte")] )
plotData = data.frame(blockData[blockData$blockNum == 2,],
                      monteAdvantage = waicList[,which(modelNames == "baseline")] - waicList[,which(modelNames == "monte")])
ggplot(plotData, aes(monteAdvantage, AUC)) + geom_point() + facet_grid(~condition)

### analysis here 
paraInfo = data.frame(muWaic = apply(waicList, MARGIN = 2, FUN = median),
                      seWaic = apply(waicList, MARGIN = 2, FUN = function(x) sd(x) / sqrt(n)),
                      muP = apply(pWaicList, MARGIN = 2, FUN = median),
                      seP = apply(pWaicList, MARGIN = 2, FUN = function(x) sd(x) / sqrt(n)),
                      modelName = modelNames)

ggplot(paraInfo, aes(modelNames, muP)) + geom_bar(stat="identity")


ggplot(paraInfo, aes(modelNames, muWaic)) + geom_bar(stat="identity")



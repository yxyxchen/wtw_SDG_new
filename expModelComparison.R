# this script extract waic for all models. 
library("stringr")
library("ggplot2")
modelNames = unlist(list.files(path="genData/expModelFitting/"))
modelNames = factor(modelNames,
                   levels = c("baseline", "monte", "monteRP", "monteSteep",
                              "monteRPSteep"),
                   ordered = T)

load("genData/expDataAnalysis/blockData.RData")
noStressIDList = unique(blockData$id[blockData$stress == "no stress"]) 
nNoStress = length(noStressIDList)

# select ID
source("subFxs/loadFxs.R")
expPara = loadExpPara("monte", c("phi", "tau", "gamma"))
a = noStressIDList[expPara$gammaRhat > 1.1 | expPara$phiRhat >1.1 | expPara$tauRhat > 1.1 | 
                     expPara$phiEffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100 ]
expPara = loadExpPara("monteSteep", c("phi", "tau", "gamma", "steep")) # numPara and nPara
b = noStressIDList[expPara$gammaRhat > 1.1 | expPara$phiRhat >1.1 | expPara$tauRhat > 1.1 | expPara$steepRhat>1.1 | 
                     expPara$phiEffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100 | expPara$steepEffe < 100]
expPara = loadExpPara("monteRP", c("phiR", "phiP", "tau", "gamma")) # numPara and nPara
c = noStressIDList[expPara$gammaRhat > 1.1 | expPara$phiRRhat >1.1 | expPara$tauRhat > 1.1 | expPara$phiPRhat>1.1 | 
                     expPara$phiREffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100 | expPara$phiPEffe < 100]

useID = noStressIDList[!noStressIDList  %in% unique(c(a, b ,c))]
nUse = length(useID)
nModel = length(modelNames)
logEvidenceList = matrix(NA, nUse, nModel)
pWaicList = matrix(NA, nUse, nModel)
for(m in 1 : nModel){
  modelName = modelNames[m]
  logEvidence = vector(length = nUse)
  pWaic = vector(length = nUse)
  for(sIdx in 1 : nUse ){
    id = useID[sIdx]
    fileName = sprintf("genData/expModelFitting/%s/s%d_waic.RData", modelName, id)
    load(fileName)
    logEvidence[sIdx] = WAIC$elpd_waic
    pWaic[sIdx] = WAIC$p_waic
  }
  logEvidenceList[,m] = logEvidence
  pWaicList[,m] = pWaic
}
waicList = -2 * logEvidenceList
f= "genData/expModelFitting/logEvidenceList.csv"
write.table(file = f, waicList, sep = ",", col.names = F, row.names = F)



### compara waic 
paraInfo = data.frame(muWaic = apply(waicList, MARGIN = 2, FUN = mean),
                      seWaic = apply(waicList, MARGIN = 2, FUN = function(x) sd(x) / sqrt(nUse)),
                      muP = apply(pWaicList, MARGIN = 2, FUN = mean),
                      seP = apply(pWaicList, MARGIN = 2, FUN = function(x) sd(x) / sqrt(nUse)),
                      modelName = modelNames)
paraInfo$minWaic = paraInfo$muWaic - paraInfo$seWaic
paraInfo$maxWaic = paraInfo$muWaic + paraInfo$seWaic

ggplot(paraInfo, aes(modelNames, muP)) + geom_bar(stat="identity") 
ggplot(paraInfo, aes(modelNames, muWaic)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin = minWaic,
                                                                                          ymax = maxWaic))


# compare waic in hist 
hist(waicList[,1] - waicList[,3])


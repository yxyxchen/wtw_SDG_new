# this script extract waic for all models. 

# libraries and scripts
library("stringr")
library("ggplot2")

# load model names
modelNames = unlist(list.files(path="genData/expModelFitting/"))
modelNames = factor(modelNames,
                   levels = c("baseline", "monte", "monteRP", "monteSteep",
                              "monteSteepExp"),
                   ordered = T)
nModel = length(modelNames)
# load experimental data
load("genData/expDataAnalysis/blockData.RData")
idList = unique(blockData$id) 
n = length(idList)

# define a function for convinence
# select useID
useID_[[i]] = vector(mode = "list", length = nModel)
useID = idList
source("subFxs/loadFxs.R")
for(i in 1 : nModel){
  modelName = modelNames[i]
  expPara = loadExpPara(modelName, getPars(modelName))
  RhatCols = which(str_detect(colnames(expPara), "hat"))[1 : length(pars)]
  EffeCols = which(str_detect(colnames(expPara), "Effe"))[1 : length(pars)]
  useID_[[i]] = idList[apply(expPara[,RhatCols] < 1.1, MARGIN = 1, sum) == length(pars) &
                   apply(expPara[,EffeCols] >100, MARGIN = 1, sum) == length(pars)]
  useID = idList[idList %in% useID_[[i]] & idList %in% useID]
}
nUse = length(useID)

# extract logEvidence per trial
# here logEvidence is on the loglikelyhood scale, so it is negtive 
# also, it accounts for model complexity 
logEvidence_ = matrix(NA, nUse, nModel)
logEvidenceSe_ = matrix(NA, nUse, nModel)
logEvidencePerAct_ = matrix(NA, nUse, nModel) # save later for model comparison 
pWaic_ = matrix(NA, nUse, nModel)
for(m in 1 : nModel){
  modelName = modelNames[m]
  for(sIdx in 1 : nUse ){
    id = useID[sIdx]
    nAction = blockData$nAction[blockData$id == id & blockData$blockNum == 1]
    fileName = sprintf("genData/expModelFitting/%s/s%d_waic.RData", modelName, id)
    load(fileName)
    logEvidence_[sIdx, m] = WAIC$elpd_waic
    logEvidenceSe_[sIdx, m] = WAIC$se_elpd_waic
    logEvidencePerAct_[sIdx, m] = WAIC$elpd_waic / nAction
    pWaic_[sIdx, m] = WAIC$p_waic
  }
}
waic_ = -2 * logEvidence_
f= "genData/expModelFitting/logEvidenceList.csv"
write.table(file = f, logEvidence_, sep = ",", col.names = F, row.names = F)



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


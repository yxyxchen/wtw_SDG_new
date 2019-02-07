source("expModelFitting.R")
modelName = "monte"
pars = c("phi", "tau", "gamma")
expModelFitting(modelName, pars)

modelName = "monteSteep"
pars = c("phi", "tau", "gamma", "steep")
expModelFitting(modelName, pars)


modelName = "monteRP"
pars = c("phiR", "phiP", "tau", "gamma")
expModelFitting(modelName, pars)

source("subFxs/loadFxs.R")
expPara = loadExpPara("monte", c("phi", "tau", "gamma"))
# 9 16 42 51 69 100
# 3 7 9 24 45 69 100 now
a = noStressIDList[expPara$gammaRhat > 1.1 | expPara$phiRhat >1.1 | expPara$tauRhat > 1.1 | 
                 expPara$phiEffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100 ]

expPara = loadExpPara("monteSteep", c("phi", "tau", "gamma", "steep")) # numPara and nPara
b = noStressIDList[expPara$gammaRhat > 1.1 | expPara$phiRhat >1.1 | expPara$tauRhat > 1.1 | expPara$steepRhat>1.1 | 
                 expPara$phiEffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100 | expPara$steepEffe < 100]
# 3no   7   9  24  42no  44no 45  52no  56no
# 65no 71no  84  92 100 105 109 110

expPara = loadExpPara("monteRP", c("phiR", "phiP", "tau", "gamma")) # numPara and nPara
c = noStressIDList[expPara$gammaRhat > 1.1 | expPara$phiRRhat >1.1 | expPara$tauRhat > 1.1 | expPara$phiPRhat>1.1 | 
                               expPara$phiREffe < 100 | expPara$tauEffe < 100 | expPara$tauEffe < 100 | expPara$phiPEffe < 100]
# 11  16  20  25  26  38  42  53  79  85  94 112

conds = blockData$condition[blockData$blockNum == 1 & blockData$id %in% unique(c(a, b ,c))]
sort(conds)

conds = blockData$condition[blockData$blockNum == 1 & blockData$id %in% a]
sort(conds)


conds = blockData$condition[blockData$blockNum == 1 & blockData$id %in% b]
sort(conds)

conds = blockData$condition[blockData$blockNum == 1 & blockData$id %in% c]
sort(conds)

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

# 
source("expModelFitting.R")
modelName = "monteRatio"
pars = c("phi", "tau", "gamma", "ratio")
source("expModelFitting.R")

modelName = "monteSteepExp"
pars = c("phi", "tau", "gamma", "steep")
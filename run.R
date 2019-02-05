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
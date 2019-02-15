getPars = function(modelName){
  if(modelName == "monte") pars = c("phi", "tau", "gamma")
  else if(modelName == "monteRP") pars = c("phiR", "phiP", "tau", "gamma")
  else if(modelName == "monteSteep") pars = c("phi", "tau", "gamma", "steep")
  else if(modelName == "monteSteepExp") pars = c("phi", "tau", "gamma", "steep")
  else return("wrong model name")
  return(pars)
}

getPars = function(modelName){
  if(modelName == "cons_arbitrary") pars = c("phi", "tau", "gamma")
  else if(modelName == "cons_theoretic") pars = c("phi", "tau", "gamma")
  else if(modelName == "cons_partFlex_gamma") pars = c("phi", "tau", "gamma", "QwaitIni")
  else if(modelName == "monteSteep") pars = c("phi", "tau", "gamma", "steep")
  else if(modelName == "monteSteepExp") pars = c("phi", "tau", "gamma", "steep")
  else return("wrong model name")
  return(pars)
}

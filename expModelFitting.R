expModelFitting = function(modelName, pars){
  dir.create("genData")
  dir.create("genData/expModelFitting")
  dir.create(sprintf("genData/expModelFitting/%s", modelName))
 #  load libraries and set environments
  options(warn=-1, message =-1) # default settings borrowed somewhere
  library('plyr'); library(dplyr); library(ggplot2);library('tidyr');library('rstan') #load libraries
  Sys.setenv(USE_CXX14=1) # making rstan working on this device 
  rstan_options(auto_write = TRUE) # default settings borrowed somewhere
  options(mc.cores = parallel::detectCores())# enable multi-core precessors 
  library("loo")
  # source scripts
  source('subFxs/modelFittingFxs.R') # for fitting single case 
  source('subFxs/loadFxs.R') # for load data
  load("wtwSettings.RData")
  library("coda") # calculate psr in modelFittingFxs
  # compile the stan model 
  dir.create(sprintf("genData/expModelFitting/%s", modelName))
  model = stan_model(file = sprintf("stanModels/%s.stan", modelName))
  
  # load expData
  allData = loadAllData()
  hdrData = allData$hdrData           
  trialData = allData$trialData       
  # list with a named element for each subject ID.
  allIDs = hdrData$ID                   # column of subject IDs
  n = length(allIDs)                    # n
  
  # extract input arguments from no stree subjects 
  load("genData/expDataAnalysis/blockData.RData")
  idList = unique(blockData$id) 
  n = length(idList)

  # loop over suvject
  for(i in 70 : 89){
      thisID = idList[[i]]
      thisTrialData = trialData[[thisID]]
      thisTrialData = thisTrialData[thisTrialData$blockNum == 1,]
      timeWaited = thisTrialData$timeWaited
      scheduledWait = thisTrialData$scheduledWait
      trialEarnings = thisTrialData$trialEarnings
      timeWaited[trialEarnings > 0] = scheduledWait[trialEarnings > 0]
      cond = unique(thisTrialData$condition)
      wIni = ifelse(cond == "HP", wInisTheory[1], wInisTheory[2]) # wIni is the theoratical initial values 
      fileName = sprintf("genData/expModelFitting/%s/s%d", modelName, thisID)
      modelFitting(cond, wIni, timeWaited, trialEarnings, scheduledWait, fileName, pars, model)
  } 
}

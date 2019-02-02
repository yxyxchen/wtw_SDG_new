expModelFitting = function(modelName, pars){
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
  load("wtwSettings.R")
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
  
  trialData = trialData[(1 : length(trialData)) %in% allIDs]
  timeWaitedList = sapply(1 : n, function(sIdx) {
    tempt = trialData[[sIdx]]
    junk = tempt$timeWaited[1 : sum(tempt$blockNum == 1)]
  })
  trialEarningsList = sapply(1 :n, function(sIdx) {
    tempt = trialData[[sIdx]]
    junk = tempt$trialEarnings[1 : sum(tempt$blockNum == 1)]
  })
  scheduledWaitList = sapply(1 :n, function(sIdx) {
    tempt = trialData[[sIdx]]
    junk = tempt$scheduledWait[1 : sum(tempt$blockNum == 1)]
  })
  condList = sapply(1 :n, function(sIdx) {
    tempt = trialData[[sIdx]]
    junk = unique(tempt$condition) 
  })
  wIniList = ifelse(condList == "HP", wInis[1], wInis[2])
  timeWaitedList = sapply(1:n, function(sIdx){
    ifelse(trialEarningsList[[sIdx]] > 0, scheduledWaitList[[sIdx]], timeWaitedList[[sIdx]])
  })
  
  # loop over suvject
  for(sIdx in 1 : 20){
    wIni = wIniList[[sIdx]]
    cond = condList[[sIdx]]
    trialEarnings= trialEarningsList[[sIdx]]
    timeWaited = timeWaitedList[[sIdx]]
    fileName = sprintf("genData/expModelFitting/monteNew/s%d", sIdx)
    modelFitting(cond, wIni, timeWaited, trialEarnings, fileName, pars, model)
  }
}

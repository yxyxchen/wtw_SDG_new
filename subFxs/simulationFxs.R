# select modelFun by modelName
getSimModelFun = function(modelName){
  if(modelName == "monte"){
    modelFun = monte
  }else if(modelName == "monteCfd"){
    modelFun = monteCfd
  }else if(modelName == "monteRP"){
    modelFun = monteRP
  }else if(modelName == "monteSteep"){
    modelFun = monteSteep
  }
  return(modelFun)
}

################ monte ######################
monte = function(para, cond, nBlock){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]

  # determine simBlockSecs
  simBlockSecs = blockSecs * nBlock
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(wInis[1], nTimeStep)
  Qquit = wInis[2]
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, simBlockSecs / iti + 1);
  vaWaits[,1] = Qwait
  vaQuits = vector(length = simBlockSecs / iti + 1);
  vaQuits[1] = Qquit
  rewardDelays = rep(0, simBlockSecs / iti + 1)
  
  # initialize totalSecs
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, simBlockSecs / iti + 1)
  timeWaited = rep(0, simBlockSecs / iti + 1)
  sellTime = rep(0, simBlockSecs / iti + 1)
  
  # loop over trials
  tIdx = 1
  while(totalSecs < simBlockSecs) {
    # sample rewardDelay
    rewardDelay = drawSample(cond)
    # calculaye available time steps
    # since we use floor there maybe 0.5 sec error (less than 90 s)
    nAvaStep = min(floor((simBlockSecs - totalSecs) / stepDuration), nTimeStep)
    
    # loop over steps
    t = 1
    while(t <= nAvaStep){
      # determine action
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # next reward 
      rewardOccur = rewardDelay <= timeTicks[t + 1] && rewardDelay > timeTicks[t] 
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
    
      # dertime whether to continue
      trialGoOn= (action == 'wait' && !rewardOccur && t < nAvaStep)
      if(!trialGoOn){
        # if the trial stops, track trialEarnings, timeWaited and rewardDelays
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, rewardDelay, ifelse(action == "quit", timeTicks[t], timeTicks[t+1]))
        rewardDelays[tIdx] = rewardDelay
        sellTime[tIdx] = totalSecs + ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
        break
      }else{
        # otherwise, continue
        t = t + 1
      }
    }# end of the trial 
    
    # update totalSecs
    totalSecs = totalSecs + iti+ ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
    
    # update Qwait and Qquit and go to the next trail if totalSecs < simBlockSecs
    if(totalSecs < simBlockSecs){
      if(nextReward == 0){
        nextWaitRateHat =  1 / sum(1  + exp((Qquit - Qwait[1])* tau))
        trialReward = nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) +
          (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration)
      }else{
        trialReward = nextReward
      }
      
      if(action == 'wait'){
        Qwait[1 : t] = (1 - phi) * Qwait[1 : t] + phi * trialReward * gamma ^ rev((0 : (t - 1 )))
      }else{
        Qquit =  (1 - phi) * Qquit + phi *  trialReward
        if(t > 1){
          Qwait[1 : (t - 1)] = (1 - phi) * Qwait[1 : (t - 1)] +
            phi * trialReward * gamma ^ rev((1 : (t - 1 )))
        }
      }
      # track vaWaits and vaQuits 
      vaWaits[,tIdx + 1] = Qwait
      vaQuits[tIdx + 1] = Qquit
      
      # update tIdx
      tIdx = tIdx + 1
    }# end of the update
    
    if(Qquit < 0 || sum(Qwait < 0) > 0){
      browser()
    }
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : tIdx,
     "trialEarnings" = trialEarnings[1 : tIdx],
     "timeWaited" = timeWaited[1 : tIdx],
     "sellTime" = sellTime[1 : tIdx], # used in wtw analysis
     "scheduledWait" = rewardDelays[1 : tIdx],
     "vaWaits" = vaWaits[, 1 : tIdx],
     "vaQuits" = vaQuits[1 : tIdx]
    )
  return(outputs)
} #end of the function

################ monteSteep ###################### 
monteSteep = function(para, cond, nBlock){
  # parse para
  phi = para[1]
  steep = para[2]
  tau = 10
  gamma = 0.95
  
  # determine simBlockSecs
  simBlockSecs = blockSecs * nBlock
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  wIni = ifelse(cond == "HP", wInis[1], wInis[2])
  
  # initialize action values
  Qwait = wIni * exp(-steep * (1 : nTimeStep - 1))
  Qquit = wIni * gamma ^(iti / stepDuration)
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, simBlockSecs / iti + 1);
  vaWaits[,1] = Qwait
  vaQuits = vector(length = simBlockSecs / iti + 1);
  vaQuits[1] = Qquit
  rewardDelays = rep(0, simBlockSecs / iti + 1)
  
  # initialize totalSecs
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, simBlockSecs / iti + 1)
  timeWaited = rep(0, simBlockSecs / iti + 1)
  sellTime = rep(0, simBlockSecs / iti + 1)
  
  # loop over trials
  tIdx = 1
  while(totalSecs < simBlockSecs) {
    # sample rewardDelay
    rewardDelay = drawSample(cond)
    # calculaye available time steps
    # since we use floor there maybe 0.5 sec error (less than 90 s)
    nAvaStep = min(floor((simBlockSecs - totalSecs) / stepDuration), nTimeStep)
    
    # loop over steps
    t = 1
    while(t <= nAvaStep){
      # determine action
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # next reward 
      rewardOccur = rewardDelay <= timeTicks[t + 1] && rewardDelay > timeTicks[t] 
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether to continue
      trialGoOn= (action == 'wait' && !rewardOccur && t < nAvaStep)
      if(!trialGoOn){
        # if the trial stops, track trialEarnings, timeWaited and rewardDelays
        trialEarnings[tIdx] = ifelse(nextReward == tokenValue, tokenValue, 0);
        timeWaited[tIdx] = ifelse(getReward, rewardDelay, ifelse(action == "quit", timeTicks[t], timeTicks[t+1]))
        rewardDelays[tIdx] = rewardDelay
        sellTime[tIdx] = totalSecs + ifelse(getReward, rewardDelay, timeWaited[tIdx])
        break
      }else{
        # otherwise, continue
        t = t + 1
      }
    }# end of the trial 
    
    # update totalSecs
    totalSecs = totalSecs + iti+ ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
  
    # update Qwait and Qquit and go to the next trail if totalSecs < simBlockSecs
    if(totalSecs < simBlockSecs){
      if(nextReward == 0){
        nextWaitRateHat =  1 / sum(1  + exp((Qquit - Qwait[1])* tau))
        trialReward = nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) +
          (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration)
      }else{
        trialReward = nextReward
      }
      
      if(action == 'wait'){
        Qwait[1 : t] = (1 - phi) * Qwait[1 : t] + phi * trialReward * gamma ^ rev((0 : (t - 1 )))
      }else{
        Qquit =  (1 - phi) * Qquit + phi *  trialReward
        if(t > 1){
          Qwait[1 : (t - 1)] = (1 - phi) * Qwait[1 : (t - 1)] +
            phi * trialReward * gamma ^ rev((1 : (t - 1 )))
        }
      }
      # track vaWaits and vaQuits 
      vaWaits[,tIdx + 1] = Qwait
      vaQuits[tIdx + 1] = Qquit
      
      # update tIdx
      tIdx = tIdx + 1
    }# end of the update
    
    if(Qquit < 0 || sum(Qwait < 0) > 0){
      browser()
    }
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : tIdx,
    "trialEarnings" = trialEarnings[1 : tIdx],
    "timeWaited" = timeWaited[1 : tIdx],
    "sellTime" = sellTime[1 : tIdx], # used in wtw analysis, cumSum timeWaited and iti
    "scheduledWait" = rewardDelays[1 : tIdx],
    "vaWaits" = vaWaits[, 1 : tIdx],
    "vaQuits" = vaQuits[1 : tIdx]
  )
  return(outputs)
} #end of the function
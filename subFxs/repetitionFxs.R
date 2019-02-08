# select modelFun by modelName
getRepModelFun = function(modelName){
  if(modelName == "monte"){
    repModelFun = monte
  }else if(modelName == "monteCfd"){
    repModelFun = monteCfd
  }else if(modelName == "monteRP"){
    repModelFun = monteRP
  }else if(modelName == "monteSteep"){
    repModelFun = monteSteep
  }
  return(repModelFun)
}


################ monte ######################
monte = function(para, cond, scheduledWait){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]
  
  # determine number of trials 
  nTrial = length(schedualeWait)
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  wIni = ifelse(cond == "HP", wInis[1], wInis[2])
  
  # initialize action values
  Qwait = rep(wIni, nTimeStep) 
  Qquit = wIni * gamma ^(iti / stepDuration)
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, nTrial);
  vaWaits[,1] = Qwait
  vaQuits = vector(length = nTrial);
  vaQuits[1] = Qquit
  rewardDelays = rep(0, nTrial)
  
  # initialize totalSecs 
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # sample rewardDelay
    rewardDelay = scheduledWait[tIdx]
    
    # loop over steps
    t = 1
    while(t <= nTimeStep){
      # determine action
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # next reward 
      rewardOccur = rewardDelay <= timeTicks[t + 1] && rewardDelay > timeTicks[t] 
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether to continue
      trialGoOn= (action == 'wait' && !rewardOccur && t < nTimeStep)
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
    
    # update Qwait and Qquit and go to the next trail if t < nTimeStep
    if(tIdx < nTrial){
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
      
    }# end of the update
    
    if(Qquit < 0 || sum(Qwait < 0) > 0){
      browser()
    }
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = rewardDelays,
    "vaWaits" = vaWaits,
    "vaQuits" = vaQuits
  )
  return(outputs)
} #end of the function


################ monte ######################
monteSteep = function(para, cond, scheduledWait){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]
  steep = para[4]
  
  # determine number of trials 
  nTrial = length(schedualeWait)
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  wIni = ifelse(cond == "HP", wInis[1], wInis[2])
  
  # initialize action values
  Qwait = wIni *(( 1 : nTimeStep - 1) * steep + 1)
  Qquit = wIni * gamma ^(iti / stepDuration)
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, nTrial);
  vaWaits[,1] = Qwait
  vaQuits = vector(length = nTrial);
  vaQuits[1] = Qquit
  rewardDelays = rep(0, nTrial)
  
  # initialize totalSecs 
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # sample rewardDelay
    rewardDelay = scheduledWait[tIdx]
    
    # loop over steps
    t = 1
    while(t <= nTimeStep){
      # determine action
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # next reward 
      rewardOccur = rewardDelay <= timeTicks[t + 1] && rewardDelay > timeTicks[t] 
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether to continue
      trialGoOn= (action == 'wait' && !rewardOccur && t < nTimeStep)
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
    
    # update Qwait and Qquit and go to the next trail if t < nTimeStep
    if(tIdx < nTrial){
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
      
    }# end of the update
    
    # if(Qquit < 0 || sum(Qwait < 0) > 0){
    #   browser()
    # }
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = rewardDelays,
    "vaWaits" = vaWaits,
    "vaQuits" = vaQuits
  )
  return(outputs)
} #end of the function


################ monte ######################
monteRP = function(para, cond, scheduledWait){
  # parse para
  phiR = para[1]
  phiP = para[2]
  tau = para[3]
  gamma = para[4]
  
  
  # determine number of trials 
  nTrial = length(schedualeWait)
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  wIni = ifelse(cond == "HP", wInis[1], wInis[2])
  
  # initialize action values
  Qwait = rep(wIni, nTimeStep) 
  Qquit = wIni * gamma ^(iti / stepDuration)
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, nTrial);
  vaWaits[,1] = Qwait
  vaQuits = vector(length = nTrial);
  vaQuits[1] = Qquit
  rewardDelays = rep(0, nTrial)
  
  # initialize totalSecs 
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, nTrial)
  timeWaited = rep(0, nTrial)
  sellTime = rep(0, nTrial)
  
  # loop over trials
  for(tIdx in 1 : nTrial) {
    # sample rewardDelay
    rewardDelay = scheduledWait[tIdx]
    
    # loop over steps
    t = 1
    while(t <= nTimeStep){
      # determine action
      waitRate =  1 / sum(1  + exp((Qquit - Qwait[t])* tau))
      action = ifelse(runif(1) < waitRate, 'wait', 'quit')
      # next reward 
      rewardOccur = rewardDelay <= timeTicks[t + 1] && rewardDelay > timeTicks[t] 
      getReward = (action == 'wait' && rewardOccur);
      nextReward = ifelse(getReward, tokenValue, 0) 
      
      # dertime whether to continue
      trialGoOn= (action == 'wait' && !rewardOccur && t < nTimeStep)
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
    
    # update Qwait and Qquit and go to the next trail if t < nTimeStep
    if(tIdx < nTrial){
      if(nextReward == 0){
        nextWaitRateHat =  1 / sum(1  + exp((Qquit - Qwait[1])* tau))
        trialReward = nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) +
          (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration)
      }else{
        trialReward = nextReward
      }
      
      if(action == 'wait'){
        delta = trialReward * gamma ^ rev((0 : (t - 1 ))) - Qwait[1 : t]
        Qwait[1 : t] =  Qwait[1 : t] + ifelse(delta > 0,  phiR, phiP) * delta
      }else{
        delta = trialReward - Qquit
        Qquit =   Qquit +  ifelse(delta > 0, phiR, phiP) * delta
        if(t > 1){
          delta = trialReward * gamma ^ rev((1 : (t - 1 ))) - Qwait[1 : (t-1)]
          Qwait[1 : (t - 1)] =  Qwait[1 : (t - 1)] + ifelse(delta > 0, phiR, phiP) * delta
        }
      }
      # track vaWaits and vaQuits 
      vaWaits[,tIdx + 1] = Qwait
      vaQuits[tIdx + 1] = Qquit
      
    }# end of the update
    
    if(Qquit < 0 || sum(Qwait < 0) > 0){
      browser()
    }
  } # end of all trials 
  
  outputs = list( 
    "trialNum" = 1 : nTrial,
    "trialEarnings" = trialEarnings,
    "timeWaited" = timeWaited,
    "sellTime" = sellTime, # used in wtw analysis
    "scheduledWait" = rewardDelays,
    "vaWaits" = vaWaits,
    "vaQuits" = vaQuits
  )
  return(outputs)
} #end of the function
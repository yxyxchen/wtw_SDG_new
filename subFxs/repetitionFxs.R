# select modelFun by modelName
getRepModelFun = function(modelName){
  if(modelName == "cons_arbitrary"){
    repModelFun = cons_arbitrary
  }else if(modelName == "cons_theoretic"){
    repModelFun = cons_theoretic
  }else if(modelName == "cons_partFlex_gamma"){
    repModelFun = cons_partFlex_gamma
  }else if(modelName == "monteSteep"){
    repModelFun = monteSteep
  }else if(modelName == "monteRatio"){
    repModelFun = monteRatio
  }else if(modelName == "monteSteepExp"){
    repModelFun = monteSteepExp
  }else if(modelName == "monteIncre"){
    repModelFun = monteIncre
  }
  return(repModelFun)
}

################ monte ######################
cons_partFlex_gamma = function(para, cond, scheduledWait){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]
  QwaitIni = para[4]
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = ifelse(cond == "HP", wInisTheory[1], wInisTheory[2])
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(QwaitIni, nTimeStep) 
  Qquit = QwaitIni * gamma ^(iti / stepDuration)
  
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
cons_theoretic = function(para, cond, scheduledWait){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  wIni = ifelse(cond == "HP", wInisTheory[1], wInisTheory[2])
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
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
cons_arbitrary = function(para, cond, scheduledWait){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(wInis[1], nTimeStep) 
  Qquit = wInis[2] 
  
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
  nTrial = length(scheduledWait)
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = wInis[1] *(( 1 : nTimeStep - 1) * steep + 1)
  Qquit = wInis[2]
  
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
  nTrial = length(scheduledWait)
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(wInis[1], nTimeStep)
  Qquit = wInis[2]
  
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


################ monte ######################
monteRatio = function(para, cond, scheduledWait, ratio){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]
  ratio = para[4]
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(wInis[1], nTimeStep)
  Qquit = wInis[1] * ratio
  
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
monteSteepExp = function(para, cond, scheduledWait){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]
  steep = para[4]
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = wInis[1] *(exp(-steep * (0 : (nTimeStep - 1))))
  Qquit = wInis[2] 
  
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

################ monteIncre ######################
monteIncre = function(para, cond, scheduledWait){
  # parse para
  tau = para[1]
  gamma = para[2]
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(wInis[1], nTimeStep) 
  Await = rep(1, nTimeStep)
  Qquit = wInis[2]
  Aquit = 1
  
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
        Qwait[1 : t] = (1 - 1 / Await[1:t]) * Qwait[1 : t] + (1 / Await[1:t]) * trialReward * gamma ^ rev((0 : (t - 1 )))
        Await[1:t] = Await[1:t] + 1    
      }else{
        Qquit =  (1 - 1 / Aquit) * Qquit + 1 / Aquit *  trialReward
        if(t > 1){
          Qwait[1 : (t - 1)] = (1 - 1 / Await[1 : (t -1)]) * Qwait[1 : (t - 1)] +
            (1 / Await[1 : (t-1)]) * trialReward * gamma ^ rev((1 : (t - 1 )))
          Await[1:(t-1)] = Await[1:(t-1)] + 1  
          Aquit = Aquit + 1
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


################ monteWini ######################
monteWini = function(para, cond, scheduledWait){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]
  wIni = para[4]
  wIni2 = para[5]
  
  # determine number of trials 
  nTrial = length(scheduledWait)
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  
  # initialize action values
  Qwait = rep(wIni, nTimeStep) 
  Qquit = wIni2
  
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

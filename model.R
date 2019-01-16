# add number of repeation
monte = function(para, cond){
  # parse para
  phi = para[1]
  tau = para[2]
  gamma = para[3]

  # load setting parameters and drawSample()
  source('subFxs/taskFxs.R')
  source("subFxs/wtwSettings.R")
  
  # determine parameters for this condition
  tMax= ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  timeTicks = seq(0, tMax, by = stepDuration)
  nTimeStep = tMax / stepDuration
  wIni = ifelse(cond == "HP", wInis[1], wInis[2])
  
  # initialize action values
  Qwait = rep(wIni, nTimeStep) 
  Qquit = wIni * gamma ^(iti / stepDuration)
  
  # initialize varibles for recording
  vaWaits = matrix(NA, nTimeStep, blockSecs / iti + 1);
  vaWaits[,1] = Qwait
  vaQuits = vector(length = blockSecs / iti + 1);
  vaQuits[1] = Qquit
  rewardDelays = rep(0, blockSecs / iti + 1)
  
  # initialize totalSecs
  totalSecs = 0
  
  # initialize outputs 
  trialEarnings = rep(0, blockSecs / iti + 1)
  timeWaited = rep(0, blockSecs / iti + 1)
  
  # loop over trials
  tIdx = 1
  while(totalSecs < blockSecs) {
    # sample rewardDelay
    rewardDelay = drawSample(cond)
    # calculaye available time steps
    # since we use floor there maybe 0.5 sec error (less than 90 s)
    nAvaStep = min(floor((blockSecs - totalSecs) / stepDuration), nTimeStep)
    
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
        break
      }else{
        # otherwise, continue
        t = t + 1
      }
    }# end of the trial 
    
    # update totalSecs
    totalSecs = totalSecs + iti+ ifelse(getReward, rewardDelay, timeWaited[tIdx]) 
    
    # update Qwait and Qquit if it is not the last trial
    if(totalSecs < blockSecs){
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
    
    # update tIdx
    tIdx = tIdx + 1
    
    if(Qquit < 0 || sum(Qwait < 0) > 0){
      browser()
    }
  } # end of all trials 
  outputs = list(
                 "trialEarnings" = trialEarnings,
                 "timeWaited" = timeWaited,
                 "rewardDelays" = rewardDelays,
                 "vaWaits" = vaWaits,
                 "vaQuits" = vaQuits)
  return(outputs)
} #end of the function


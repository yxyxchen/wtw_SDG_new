
data {
  // depending on the condition
  real wIni;
  int tMax;
  int nTimeStep; // since round returns real here, so nTimeStep != tMax / stepDuration
  
  // depending on each subject
  int N; // number of trials
  vector[N] timeWaited;
  vector[N] trialEarnings;
  int nTimePoints[N]; // list of available time points 
}
transformed data {
  // constant
  real stepDuration = 0.5;
  real iti = 2;
  real tokenValue = 10;
  }
  parameters {
  real<lower = 0, upper = 1> phiR;
  real<lower = 0, upper = 1> phiP;
  real<lower = 1, upper = 30> tau;
  real<lower = 0, upper = 1> gamma;
}
transformed parameters{
  // initialize action values 
  vector[nTimeStep] Qwait = rep_vector(wIni, nTimeStep);
  real Qquit = wIni * gamma ^(iti / stepDuration);
  
  // initialize recordings of action values 
  matrix[nTimeStep, N] Qwaits = to_matrix(rep_vector(0, nTimeStep * N), nTimeStep, N);
  vector[N] Qquits = rep_vector(0, N);
  
  
  // initialize trialReward and nextWaitRateHat
  real trialReward;
  real nextWaitRateHat = 1.0;

  // define gamma List
  vector[nTimeStep] gammaList;
  for(i in 1 : nTimeStep){
    gammaList[i] = gamma ^ (nTimeStep - i);//larger and larger
  }
  
  // fill the first trial of Qwaits and Quits
  Qwaits[,1] = Qwait;
  Qquits[1] = Qquit;

  //loop over trial
  for(tIdx in 1 : (N -1)){
    // determine nTimePoint
    int nTimePoint = nTimePoints[tIdx]; 
    // update and track action values
    if(trialEarnings[tIdx] > 0){
      trialReward = tokenValue;
      for(t in 1 : nTimePoint){
	if(trialReward * gammaList[nTimeStep - nTimePoint + t] > Qwait[t]){
		Qwait[t] = (1 - phiR) * Qwait[t] + phiR * trialReward * gammaList[nTimeStep - nTimePoint + t];
	}else{
		Qwait[t] = (1 - phiP) * Qwait[t] + phiP * trialReward * gammaList[nTimeStep - nTimePoint + t];
	}
      }
    }else{
      nextWaitRateHat =  1 / (1  + exp((Qquit - Qwait[1])* tau));
      trialReward = nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) + (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration);
      if(trialReward > Qquit){
      	Qquit =  (1 - phiR) * Qquit + phiR *  trialReward;
      }else{
      	Qquit =  (1 - phiP) * Qquit + phiP *  trialReward;
      }
      // update Qwait if nTimePoint > 1
      if(nTimePoint > 1){
        for(t in 1 : (nTimePoint - 1)){
          if(trialReward * gammaList[nTimeStep - nTimePoint + t] > Qwait[t]){
            Qwait[t] = (1 - phiR) * Qwait[t] + phiR * trialReward * gammaList[nTimeStep - nTimePoint + t];
	   }else{
            Qwait[t] = (1 - phiP) * Qwait[t] + phiP * trialReward * gammaList[nTimeStep - nTimePoint + t];
	   }
        }
      }
    }
    Qwaits[,tIdx+1] = Qwait;
    Qquits[tIdx+1] = Qquit;
  }// end of the loop
}
model {
  phiR ~ uniform(0, 1);
  phiP ~ uniform(0, 1);
  tau ~ uniform(1, 30);
  gamma ~ uniform(0, 1);
  
  // calculate the likelihood 
  for(tIdx in 1 : N){
    int action;
    vector[2] values;
    for(i in 1 : nTimePoints[tIdx]){
    if(trialEarnings[tIdx] == 0 && i == nTimePoints[tIdx]){
      action = 2; // quit
    }else{
      action = 1; // wait
    }
      values[1] = Qwaits[i, tIdx] * tau;
      values[2] = Qquits[tIdx] * tau;
      action ~ categorical_logit(values);
    } 
  }
}
generated quantities {
// initialize log_lik
  matrix[nTimeStep, N] log_lik = to_matrix(rep_vector(0, nTimeStep * N), nTimeStep, N);
  vector[2] values;
  real LL_all;
  // loop over trials
  for(tIdx in 1 : N){
    int action;
    for(i in 1 : nTimePoints[tIdx]){
      if(trialEarnings[tIdx] == 0 && i == nTimePoints[tIdx]){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      values[1] = Qwaits[i, tIdx] * tau;
      values[2] = Qquits[tIdx] * tau;
      log_lik[i, tIdx] =categorical_logit_lpmf(action | values);
    }
  }// end of the loop
  LL_all =sum(log_lik);
}




data {
  // depending on the condition
  real wInis[2];
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
  real<lower = 0, upper = 0.3> phi;
  real<lower = 2, upper = 22> tau;
  real<lower = 0.7, upper = 1> gamma;
  real<lower = 0.2, upper = 1> ratio; 
}
transformed parameters{
  // initialize action values 
  vector[nTimeStep] Qwait = rep_vector(wIni, nTimeStep);
  real Qquit = wIni * ratio;
  
  // initialize recordings of action values 
  matrix[nTimeStep, N] Qwaits = to_matrix(rep_vector(0, nTimeStep * N), nTimeStep, N);
  vector[N] Qquits = rep_vector(0, N);
  
  
  // initialize trialReward and nextWaitRateHat
  real trialReward;
  real nextWaitRateHat = 1.0;

  // define gamma List
  vector[nTimeStep] gammaList;
  for(i in 1 : nTimeStep){
    gammaList[i] = gamma ^ (nTimeStep - i);
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
      Qwait[1 : nTimePoint] = (1 - phi) * Qwait[1 : nTimePoint] + phi * trialReward * gammaList[(nTimeStep - nTimePoint + 1):nTimeStep];
    }else{
      nextWaitRateHat =  1 / (1  + exp((Qquit - Qwait[1])* tau));
      trialReward = nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) + (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration);
      Qquit =  (1 - phi) * Qquit + phi *  trialReward;
      if(nTimePoint > 1){
        Qwait[1 : (nTimePoint - 1)] = (1 - phi) * Qwait[1 : (nTimePoint - 1)] + phi * trialReward * gammaList[(nTimeStep - nTimePoint + 1):(nTimeStep - 1)];
      }
    }
    Qwaits[,tIdx+1] = Qwait;
    Qquits[tIdx+1] = Qquit;
  }// end of the loop
}
model {
  phi ~ uniform(0, 0.3);
  tau ~ uniform(2, 22);
  gamma ~ uniform(0.7, 1);
  ratio ~ uniform(0.2, 1);
  
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
      target += categorical_logit_lpmf(action | values);
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




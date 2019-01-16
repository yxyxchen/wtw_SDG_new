
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
  real tau = 100;
  }
  parameters {
  real<lower = 0, upper = 1> phi;
  real<lower = 1, upper = 30> c;
  real<lower = 0, upper = 1> gamma;
}
transformed parameters{
  // initialize action values 
  vector[nTimeStep] Qwait = rep_vector(wIni, nTimeStep);
  vector[nTimeStep] QwaitCount = rep_vector(0, nTimeStep);
  real Qquit = wIni * gamma ^(iti / stepDuration);
  real QquitCount = 0;
  
  // initialize recordings of action values 
  matrix[nTimeStep, N] Qwaits = to_matrix(rep_vector(0, nTimeStep * N), nTimeStep, N);
  matrix[nTimeStep, N] QwaitCounts = to_matrix(rep_vector(0, nTimeStep * N), nTimeStep, N);
  vector[N] Qquits = rep_vector(0, N);
  vector[N] QquitCounts = rep_vector(0, N);
  
  // initialize trialReward and nextWaitRateHat
  real trialReward;
  real nextWaitRateHat = 1.0;
  // generate gammaList
  vector[nTimeStep] gammaList;
  for(i in 1 : nTimeStep){
    gammaList[i] = gamma ^ (nTimeStep - i);
  }

  // fill the first trial of Qwaits and Quits
  Qwaits[,1] = Qwait;
  Qquits[1] = Qquit;
  QwaitCounts[,1] = QwaitCount;
  QquitCounts[1] = QquitCount;
  
  //loop over trial
  for(tIdx in 1 : (N-1)){
    // determine nTimePoint
    int nTimePoint = nTimePoints[tIdx]; 
    // update and track action values
    if(trialEarnings[tIdx] > 0){
      trialReward = tokenValue;
      Qwait[1 : nTimePoint] = (1 - phi) * Qwait[1 : nTimePoint] + phi * trialReward * gammaList[(nTimeStep - nTimePoint + 1):nTimeStep];
      QwaitCount[1:nTimePoint] = QwaitCount[1:nTimePoint] + 1; // update count
    }else{
      nextWaitRateHat =  1 / (1  + exp((Qquit - Qwait[1])* tau));
      trialReward = nextWaitRateHat * Qwait[1] * gamma ^(iti / stepDuration) + (1 - nextWaitRateHat) * Qquit * gamma ^(iti / stepDuration);
      Qquit =  (1 - phi) * Qquit + phi *  trialReward;
      QquitCount = QquitCount + 1;
      if(nTimePoint > 1){
        Qwait[1 : (nTimePoint - 1)] = (1 - phi) * Qwait[1 : (nTimePoint - 1)] + phi * trialReward * gammaList[(nTimeStep - nTimePoint + 1):(nTimeStep - 1)];
        QwaitCount[1:(nTimePoint-1)] = QwaitCount[1:(nTimePoint - 1)] + 1;
      }
    }
    Qwaits[,tIdx+1] = Qwait;
    Qquits[tIdx+1] = Qquit;
    QwaitCounts[,tIdx+1] = QwaitCount;
    QquitCounts[tIdx + 1] = QquitCount;
  }// end of the loop
}
model {
  phi ~ uniform(0, 1);
  c ~ uniform(1, 30);
  gamma ~ uniform(0, 1);
  
  // count the number of actions
  int t = 0;
  // calculate the likelihood 
  for(tIdx in 1 : N){
    int action;
    vector[2] values;
    for(i in 1 : nTimePoints[tIdx]){
      t = t + 1;
      if(trialEarnings[tIdx] == 0 && i == nTimePoints[tIdx]){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
        values[1] = (Qwaits[i, tIdx] + c * sqrt(log(t) / QwaitCounts[i, tIdx])) * tau;
        values[2] = (Qquits[tIdx]  + c * sqrt(log(t) / QquitCounts[tIdx])) * tau;
        action ~ categorical_logit(values);
    } 
  }
}
generated quantities {
// initialize log_lik
  matrix[nTimeStep, N] log_lik = to_matrix(rep_vector(0, nTimeStep * N), nTimeStep, N);
  vector[2] values;
  real LL_all;
  // count the number of actions
  int t = 0;
  // loop over trials
  for(tIdx in 1 : N){
    int action;
    for(i in 1 : nTimePoints[tIdx]){
      t = t + 1;
      if(trialEarnings[tIdx] == 0 && i == nTimePoints[tIdx]){
        action = 2; // quit
      }else{
        action = 1; // wait
      }
      values[1] = (Qwaits[i, tIdx] + c * sqrt(log(t) / QwaitCounts[i, tIdx])) * tau;
      values[2] = (Qquits[tIdx]  + c * sqrt(log(t) / QquitCounts[tIdx])) * tau;
      log_lik[i, tIdx] =categorical_logit_lpmf(action | values);
    }
  }// end of the loop
  LL_all =sum(log_lik);
}





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
  real<lower = 0, upper = 1> waitRate;
}
model {
  waitRate ~ uniform(0, 1);
  // calculate the likelihood 
  for(tIdx in 1 : N){
    int action;
    for(i in 1 : nTimePoints[tIdx]){
    if(trialEarnings[tIdx] == 0 && i == nTimePoints[tIdx]){
      action = 0; // quit
    }else{
      action = 1; // wait
    }
      action ~ bernoulli(waitRate);
    } 
  }
}
generated quantities {
// initialize log_lik
  matrix[nTimeStep, N] log_lik = to_matrix(rep_vector(0, nTimeStep * N), nTimeStep, N);
  real LL_all;
  // loop over trials
  for(tIdx in 1 : N){
    int action;
    for(i in 1 : nTimePoints[tIdx]){
      if(trialEarnings[tIdx] == 0 && i == nTimePoints[tIdx]){
        action = 0; // quit
      }else{
        action = 1; // wait
      }
      log_lik[i, tIdx] = bernoulli_lpmf(action | waitRate);
    }
  }// end of the loop
  LL_all =sum(log_lik);
}




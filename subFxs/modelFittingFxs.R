# set up 
modelFitting = function(cond, wIni, timeWaited, trialEarning, fileName, pars, model){
  tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  condIdx = ifelse(cond =="HP", 1, 2)
  nChain = 3
  nIter = 3000
  nTimeStep = tMax / stepDuration
  nTimePoints = round(ifelse(trialEarnings >0, ceiling(timeWaited / stepDuration), floor(timeWaited / stepDuration) + 1))
  data_list <- list(tMax = tMax,
                    wIni = wIni,
                    nTimeStep = nTimeStep,
                    N = length(timeWaited),
                    timeWaited = timeWaited,
                    trialEarnings = trialEarnings,
                    nTimePoints = nTimePoints)
  # init = list(list('phi' = 0.5, 'tau' = 15, 'gamma' = 0.5),
  #             list('phi' = 0.1, 'tau' = 5, 'gamma' = 0.1))
  fit = sampling(object = model, data = data_list, cores = nChain, chains = nChain,
               iter = nIter, 
               show_messages = F) 
  # extract parameters
  extractedPara = fit %>%
    rstan::extract(permuted = F, pars = c(pars, "LL_all"))
  # save sampling sequences
  tempt = extractedPara %>%
    adply(2, function(x) x) %>%  # change arrays into 2-d dataframe 
    select(-chains) 
  write.csv(tempt, file = sprintf("%s.txt", fileName), row.names=FALSE)
  # calculate and save WAIC
  log_lik = extract_log_lik(fit)
  WAIC = waic(log_lik)
  save("WAIC", file = sprintf("%s.RData", fileName))
  # calculate potential scale reduction 
  n = nInter / 2 # number of samples used in each seq
  m = nChain
grandMean = extractedPara %>% apply(FUN = mean, MARGIN = c(3))
seqMean = extractedPara %>% apply(FUN = mean, MARGIN = c(2,3))
seqVar = extractedPara %>% apply(FUN = var, MARGIN = c(2,3)) # sample variance with n-1 degree
B = seqMean %>% apply(FUN = var, MARGIN = 2) * n
W =  seqVar %>% apply(FUN = function(x) sum(x) / m, MARGIN = 2)
sigmaSquare = (n-1) / n * W + B / n
V = sigmaSquare + B/(nChain * n)
covS2X2 = mapply(FUN = cov, lapply(1 : length(pars), function(i) seqVar[,i]),
                 lapply(1 : length(pars), function(i) seqMean[,i]^2))
covS2X = mapply(FUN = cov, lapply(1 : length(pars), function(i) seqVar[,i]),
                 lapply(1 : length(pars), function(i) seqMean[,i]))
VVar = ((n-1) / n) ^ 2 /  m * apply(seqVar, FUN = var, MARGIN = 2) +
  ((m+1) / m / n)^2 *2 /(m - 1) * B ^ 2 +
  2 * (m+1) * (n-1) / (m * n^2) * n / m * (covS2X2 - 2 * grandMean * covS2X)
df =  2 * V^2 / var(V)
R = V / W * df / (df - 2)
}







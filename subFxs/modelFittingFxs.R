# set up 
modelFitting = function(cond, wIni, timeWaited, trialEarning, fileName, pars, model){
  tMax = ifelse(cond == "HP", tMaxs[1], tMaxs[2])
  condIdx = ifelse(cond =="HP", 1, 2)
  nChain = 4
  nIter = 5000
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
  fit = sampling(object = model, data = data_list, cores = min(nChain, 3), chains = nChain,
               iter = nIter) 
  # extract parameters
  extractedPara = fit %>%
    rstan::extract(permuted = F, pars = c(pars, "LL_all"))
  # save sampling sequences
  tempt = extractedPara %>%
    adply(2, function(x) x) %>%  # change arrays into 2-d dataframe 
    select(-chains) 
  write.csv(tempt, file = sprintf("%s.txt", fileName), row.names=FALSE)
  # calculate and save WAIC
  log_lik = extract_log_lik(fit) # quit time consuming 
  WAIC = waic(log_lik)
  save("WAIC", file = sprintf("%s_waic.RData", fileName))
  # calculate potential scale reduction 
  chain1 = mcmc(data= extractedPara[,1,])
  chain2 = mcmc(data= extractedPara[,2,])
  chain3 = mcmc(data= extractedPara[,3,])
  chain4 = mcmc(data= extractedPara[,4,])
  mcOb = mcmc.list(chain1, chain2, chain3, chain4)
  psr = gelman.diag(mcOb, confidence = 0.95, transform=FALSE, autoburnin=F,
                  multivariate=TRUE)$psrf
  save("psr", file = sprintf("%s_psr.RData", fileName))
}



# set up 
modelFitting = function(cond, wIni, timeWaited, trialEarnings, fileName, pars, model){
  load("wtwSettings.RData")
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
  write.table(matrix(unlist(tempt), ncol = length(pars) + 1), file = sprintf("%s.txt", fileName), sep = ",",
              col.names = F, row.names=FALSE)
  # calculate and save WAIC
  log_lik = extract_log_lik(fit) # quit time consuming 
  WAIC = waic(log_lik)
  save("WAIC", file = sprintf("%s_waic.RData", fileName))
  # save summarized fit 
  fitSumary <- summary(fit,pars = c(pars, "lp__", "LL_all"), use_cache = F)$summary
  write.table(matrix(fitSumary, nrow = length(pars) + 2), file = sprintf("%s_summary.txt", fileName),  sep = ",",
            col.names = F, row.names=FALSE)
}



# monte
nPara = 3
paraNames = c('phi', 'tau', 'gamma')
nValue = 5
nComb = nValue ^ nPara
realParas = matrix(NA, nValue^nPara, nPara)
realParas[,1] = rep(seq(0.01, 0.05, length.out = 5), nValue^(nPara - 1)) # phi
realParas[,2] = rep(rep(seq(10, 22, length.out = 5), each = nValue), nValue^(nPara - 2)) # tau
realParas[,3] = rep(seq(0.8, 0.98, length.out = 5), each = nValue^2)
save("realParas", "nComb", "nValue", "nPara", "paraNames", file = "genData/simulation/monte/realParas.RData")

# monteCfd
nPara = 3
paraNames = c('phi', 'c', 'gamma')
nValue = 5
nComb = nValue ^ nPara
realParas = matrix(NA, nValue^nPara, nPara)
realParas[,1] = rep(seq(0.01, 0.05, length.out = 5), nValue^(nPara - 1)) # phi
realParas[,2] = rep(rep(seq(0, 10, length.out = 5), each = nValue), nValue^(nPara - 2)) # tau
realParas[,3] = rep(seq(0.8, 0.98, length.out = 5), each = nValue^2)
dir.create("genData/simulation/monteCfd")
save("realParas", "nComb", "nValue", "nPara", "paraNames", file = "genData/simulation/monteCfd/realParas.RData")
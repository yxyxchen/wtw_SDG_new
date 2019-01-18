# this script extract waic for all models. 
library("R.matlab")
n = 120
nModel = 2
<<<<<<< HEAD
modelNames = c("monte", "monteBias")
=======
modelNames = c("monte", "monteRP")
>>>>>>> 76dded834aafefccba54ad4be637489b7f59cbba

waicList = matrix(NA, n, nModel)

for(m in 1 : nModel){
  modelName = modelNames[m]
  waic = vector(length = n)
  for(sIdx in 1 : n){
    fileName = sprintf("genData/expModelFitting/%s/s%d.RData", modelName, sIdx)
    load(fileName)
    waic[sIdx] = WAIC$waic
  }
  waicList[,m] = waic
}
f= "genData/expModelFitting/waicList.mat"
writeMat(f, waicList = waicList)


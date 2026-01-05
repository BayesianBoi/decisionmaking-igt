IGTbehavAdv_v2 = function( inputData, numTrials = 100 )  {
# calculates adv. deck choice probs in each subjec
# inputData: input data across all subjects
# Caution!! colname for Deck = "Deck", colname for subjID = "subjID"
#   Assumes each subject has the same number of trials (e.g., 100)
# numTrials: number of trials
  
  #numSubjs = dim(inputData)[1]/numTrials
  numSubjs = length( unique( inputData$subjID))
  
  numDecks = 4
  #numTrials = 100
  numBlocks = 5
  numTrialsPblock = numTrials / numBlocks
  
  advProb = array( NA, c(numSubjs, numBlocks) )
  allSubjs = unique(inputData[,"subjID"])
  
  for (subjIdx in 1:numSubjs) {
    currSubj = allSubjs[subjIdx]
    currData = subset( inputData, subjID == currSubj)
    lastTrial = dim(currData)[1]  # last trial
    for (blIdx in 1:numBlocks) {
      if (blIdx < numBlocks) {
        currBl = currData[(1+(blIdx-1)*numTrialsPblock):(blIdx*numTrialsPblock), ]
          
        tmpProb = ( sum(currBl[,"deck"]==3) + sum(currBl[,"deck"]==4) ) / numTrialsPblock  
      } else {
        currBl = currData[(1+(blIdx-1)*numTrialsPblock):lastTrial, ]  # use only up to lastTrial (e.g., 91...)
        tmpTrials = dim(currBl)[1]
        tmpProb = ( sum(currBl[,"deck"]==3) + sum(currBl[,"deck"]==4) ) / tmpTrials  
      }
      
      advProb[subjIdx, blIdx] = tmpProb
    }    
  }

return( advProb )
}


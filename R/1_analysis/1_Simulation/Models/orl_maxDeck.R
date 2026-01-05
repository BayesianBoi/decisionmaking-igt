sampleIdx <- sample(samp_vec, 1)

# True estimated parameters used to simulate data
Arew  <- simul_pars$Arew[sampleIdx,subj]
Apun  <- simul_pars$Apun[sampleIdx,subj]
K     <- simul_pars$K[sampleIdx,subj]
betaF <- simul_pars$betaF[sampleIdx,subj]
betaP <- simul_pars$betaP[sampleIdx,subj]

# Initialize values
ev        <- c(0, 0, 0, 0)               # Initial expected values
ef        <- c(0, 0, 0, 0)               # Initial expected frequency values
pers      <- c(0, 0, 0, 0)               # Initial perseverance values
pChoose   <- c(0.25, 0.25, 0.25, 0.25)   # Initial probability of choosing each deck
deckCount <- c(0, 0, 0, 0)               # Initial deck choice count 
use_decks <- c(1, 2, 3, 4)

for (t in 1:n_trials) {
  
  # Probability of selecting deck
  for (d in 1:4) {
    pChoose[d] <- exp( (ev[d] + ef[d] * betaF + pers[d] * betaP) ) / sum( exp( (ev[use_decks] + ef[use_decks] * betaF + pers[use_decks] * betaP) ) )
    if (!(d %in% use_decks)) {
      pChoose[d] <- 0
    }
  }
  
  # Store probability of choosing deck in array
  pDeck_subj[t,, iter] <- pChoose
  
  # Selected deck given probability of selecting deck
  tmpChoice = sample(1:4, size = 1, replace = T, prob = pChoose  )
  # Add choice count to deck
  deckCount[tmpChoice] <- deckCount[tmpChoice] + 1
  
  # Amount won
  gain <- all_decks[deckCount[tmpChoice], 1, tmpChoice]
  # Amount lost
  loss <- all_decks[deckCount[tmpChoice], 2, tmpChoice]
  # Total outcome
  outcome <- gain - abs(loss)
  
  # Binary outcome 
  if ( outcome > 0 ) {
    binOutcome <- 1
  } else if ( outcome < 0) {
    binOutcome <- -1
  } else {
    binOutcome <- 0
  }
  
  # Prediction error
  PEval  <- outcome - ev[tmpChoice]
  PEfreq <- binOutcome - ef[tmpChoice]
  PEfreq_fic  <- -binOutcome/3 - ef[-tmpChoice]
  
  # Value update
  if ( outcome >= 0) {  # x(t) >= 0
    ev[tmpChoice] <- ev[tmpChoice] + Arew * PEval     # expected value update
    ef[tmpChoice] <- ef[tmpChoice] + Arew * PEfreq    # expected frequency update
    ef[-tmpChoice] <- ef[-tmpChoice] + Apun * PEfreq_fic   # fictive frequency update
  } else {              # x(t) < 0
    ev[tmpChoice] <- ev[tmpChoice] + Apun * PEval     # expected value update
    ef[tmpChoice] <- ef[tmpChoice] + Apun * PEfreq    # expected frequency update
    ef[-tmpChoice] <- ef[-tmpChoice] + Arew * PEfreq_fic   # fictive frequency update
  }
  
  # Update perseverance 
  pers[tmpChoice] <- 1   # perseverance term
  pers <- pers / (1 + K); # decay
  
  # useable decks
  if (any(deckCount >= max_trials)) {
    use_decks <- dplyr::setdiff(1:4, c(1:4)[which(deckCount>=max_trials)])
  }
}
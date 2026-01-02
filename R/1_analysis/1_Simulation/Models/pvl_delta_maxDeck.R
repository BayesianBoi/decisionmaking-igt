sampleIdx <- sample(samp_vec, 1)

# True estimated parameters used to simulate data
A      <- simul_pars$A[sampleIdx,subj]
alpha  <- simul_pars$alpha[sampleIdx,subj]
cons   <- simul_pars$cons[sampleIdx,subj]
lambda <- simul_pars$lambda[sampleIdx,subj]

# Initialize values
ev        <- c(0, 0, 0, 0)               # Initial expected values
pChoose   <- c(0.25, 0.25, 0.25, 0.25)   # Initial probability of choosing each deck
deckCount <- c(0, 0, 0, 0)               # Initial deck choice count 
use_decks <- c(1, 2, 3, 4)

for (t in 1:n_trials) {
  
  # Choice sensitivity
  sens <- 3**cons -1
  
  # Probability of selecting deck
  for (d in 1:4) {
    pChoose[d] <- exp( (ev[d] * sens) ) / sum( exp( (ev[use_decks]  * sens) ) )
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
  
  # Utility of outcome
  if ( outcome >= 0) {  # x(t) >= 0
    curUtil = outcome**alpha
  } else {              # x(t) < 0
    curUtil = -1 * lambda *  (-1*outcome)**alpha
  }
  
  # Prediction error and expected value update
  PE <- curUtil - ev[tmpChoice]; 
  ev[tmpChoice] <- ev[tmpChoice] + A * PE 
  
  # useable decks
  if (any(deckCount >= max_trials)) {
    use_decks <- dplyr::setdiff(1:4, c(1:4)[which(deckCount>=max_trials)])
  }
}

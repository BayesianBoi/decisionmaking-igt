stan_rhat <- function (fit = NULL, less = NULL) 
{
  summaryData <- rstan::summary(fit)
  rhatData <- data.frame(Rhat = summaryData[["summary"]][, 
                                                         "Rhat"])
  if (!is.null(less)) {
    if (all(rhatData$Rhat <= less)) {
      cat("TRUE: All Rhat values are less than ", less, 
          sep = "")
    }
    else {
      cat("FALSE: Some Rhat values are greater than ", 
          less, sep = "")
    }
  }
  else {
    return(rhatData)
  }
}
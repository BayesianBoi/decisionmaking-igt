# Plot the differnce between all parameters in models. This will only work if hBayesDM is installed!
plot_diff <- function (..., data = FALSE, len = FALSE, cols = FALSE, title = NULL, xlim = NULL, bins = rep(30,10)) {
  if (data == FALSE) {
    data <- list(...)
  }
  tmp1 <- data[[1]]
  tmp2 <- data[[2]]
  rm(data)
  if (len == F) {
    len = min(dim(tmp1[[1]]), dim(tmp2[[1]]))
  }
  numVars <- length(tmp1)
  if (cols == F) {
    cols <- round(numVars/2)
  }
  plotVec <- vector("list", length = numVars)
  
  if (is.null(title)) {
    title <- paste0(names(tmp2[[i]]))
  } else {
    title <- title
  }
  for (i in 1:numVars) {
    tmpDist <- tmp1[[i]][1:len] - tmp2[[i]][1:len]
    plotVec[[i]] <- plotHDI(tmpDist, Title = title[i], xlim = xlim[[i]], binSize = bins[i]) + 
      ggplot2::geom_vline(xintercept = 0, lty = 2, color = "red")
  }
  rm(tmp1, tmp2)
  cowplot::plot_grid(plotlist = plotVec, ncol = cols)
}
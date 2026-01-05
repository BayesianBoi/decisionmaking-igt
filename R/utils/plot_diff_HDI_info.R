plotHDI_info <- function (sample = NULL, credMass = 0.95, Title = NULL, xLab = "Value", 
          yLab = "Density", fontSize = NULL, binSize = 30, ...) 
{
  ..density.. <- NULL
  HDI <- HDIofMCMC(as.vector(t(sample)), credMass = credMass)
  sample_df <- data.frame(sample)
  h1 <- ggplot(sample_df, aes(x = sample)) +
    geom_histogram(aes(y = ..density..), colour = "black", 
                   fill = "grey", bins = binSize, ...) +
    geom_segment(aes(x = HDI[1], y = 0, xend = HDI[2], yend = 0), size = 1.5, colour = "red") + 
    theme_minimal(base_size = 15) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
    theme(axis.title.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          plot.title = ggplot2::element_blank())
  # print(paste0(credMass * 100, "% Highest Density Interval (HDI):"))
  # print(paste0("Lower bound=", round(HDI[1], 4), ", Upper bound=", 
  #              round(HDI[2], 4)))
  print(paste0(Title, " Greater: ", mean(sample_df>0)))
  print(paste0(Title, " Lesser: ", mean(sample_df<0)))
  return(h1)
}
plot_diff_info <- function (..., data = FALSE, len = FALSE, cols = FALSE, title = NULL) {
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
    plotVec[[i]] <- plotHDI_info(tmpDist, Title = title[i]) + 
      ggplot2::geom_vline(xintercept = 0, lty = 2, size = 1, color = "red")
  }
  rm(tmp1, tmp2)
  cowplot::plot_grid(plotlist = plotVec, ncol = cols)
}
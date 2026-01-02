library(zoo)

# Raw data to smooth
raw <- rnorm(100, 0, 1) + seq(1, 10, length.out = 100)

# Ahn 2008 method
nt <- 100
smooth_1 <- vector(length = nt)
fsize = 3     # (moving) filter size (2*fsize(=3)+1 = 7)
for (i in 1:nt)  {
  if ( i - fsize < 1 )  {
    mstart = 1
  } else {
    mstart = i - fsize
  }
  if ( i + fsize > nt )  {
    mend = nt
  } else {
    mend = i + fsize
  }
  smooth_1[i] = mean( raw[mstart : mend] )
}

# Rollapply method
smooth_2 <- rollapply(raw, FUN = mean, width = 7, partial = T, align = "center")

plot(smooth_1, smooth_2); cor(smooth_1, smooth_2)

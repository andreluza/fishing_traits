

###  kernel function from Cooke et al. 2019 (check it out here : https://github.com/03rcooke/hyper_pca/commit/2b7df79a30242d3d479e75382a8865df3f5a6f7d)

cl <- function(df, prob) {
  dx <- diff(df$x[1:2])
  dy <- diff(df$y[1:2])
  sz <- sort(df$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}

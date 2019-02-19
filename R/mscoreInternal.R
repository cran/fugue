mscoreInternal<-function (ymat, inner , trim)
{
  nrows<-dim(ymat)[1]
  ncols<-dim(ymat)[2]
  out <- matrix(NA, nrows, ncols)
  one <- rep(1, ncols - 1)
  difs <- array(NA, c(nrows, ncols, ncols - 1))
  for (j in 1:ncols) {
    difs[, j, ] <- outer(as.vector(unlist(ymat[, j])), one,
                         "*") - ymat[, -j]
  }
  ms <- as.vector(difs)
  ab <- pmin(trim, pmax(0, trim*(abs(ms) - inner))/(trim - inner))
  ms <- sign(ms) * ab

  ms <- array(ms, c(nrows, ncols, ncols - 1))
  ms <- apply(ms, c(1, 2), sum, na.rm = TRUE)
  ms[is.na(ymat)] <- NA
  ms
}

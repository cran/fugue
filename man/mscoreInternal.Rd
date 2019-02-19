\name{mscoreInternal}
\alias{mscoreInternal}

\title{
Computes M-scores for M-tests and estimates.
}
\description{
Of very limited interest to most users, function mscoreInternal() computes the M-scores used by
functions sen().
}
\usage{
mscoreInternal(ymat, inner, trim)
}

\arguments{
  \item{ymat}{
A matrix of outcomes scaled for use in an M-statistic; see the
discussion of the parameter lambda in the documentation for the sen
function.
If the largest matched set has K controls, and there are I matched sets,
then ymat is an I x (K+1) matrix.  Each row is a matched set.  The first
column contains the treated individual in the matched set.  The remaining
columns contain the controls.  If a set has fewer than K controls, then
its last columns are NAs.
}
  \item{inner}{
inner is the inner[i] parameter described in the documentation for sen().
}
  \item{trim}{
trim is the trim[i] parameter described in the documentation for sen().
}
}

\value{
Generally, a matrix with the same dimensions as ymat containing the M-scores.
}
\references{
Rosenbaum, P. R. (2015).  Two R packages for sensitivity analysis
in observational studies.  Observational Studies, v. 1. (Free on-line.)
}
\author{
Xinran Li and Paul R. Rosenbaum
}



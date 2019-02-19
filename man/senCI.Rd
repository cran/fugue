\name{senCI}
\alias{senCI}
\concept{Sensitivity analysis}
\concept{Confidence interval}
\concept{Observational study}
\concept{Causal inference}
\concept{Matching}
\concept{rbounds}
\concept{Fugue statistic}
\title{
Sensitivity Analysis for Point Estimates and Confidence Intervals in an Observational Study.
}
\description{
Each matched set contains one treated
individual and one or more controls.
Uses Huber's M-statistic as the basis for
the test; see Maritz (1979).  Matched sets of different sizes
use different \eqn{\psi}-functions, creating what is called a fugue statistic.
Performs either a randomization
test (\eqn{\Gamma=1}) or an analysis of sensitivity to departures from random
assignment (\eqn{\Gamma>1}).  For hypothesis tests, use function sen().
The method is described in Li and Rosenbaum (2019); see also Rosenbaum (2007,2013).
}
\usage{
senCI(y, z, mset, gamma = 1, inner = NULL, trim = NULL, lambda = 1/2,
     alpha = 0.05, alternative = "greater")
}
\arguments{
  \item{y}{ A vector of responses with no missing data.}
  \item{z}{ Treatment indicator, z=1 for treated, z=0 for control with
  length(z)==length(y).
}
  \item{mset}{ Matched set indicator, 1, 2, ..., sum(z) with
  length(mset)==length(y).  Matched set indicators should be either integers or a factor.
}
  \item{gamma}{
gamma is the sensitivity parameter \eqn{\Gamma}, where \eqn{\Gamma \ge 1}.  Setting
\eqn{\Gamma = 1} is equivalent to assuming ignorable treatment assignment given
the matched sets, and it performs a within-set randomization test.
}
  \item{inner}{
inner and trim together define the \eqn{\psi}-function for the M-statistic.
If the largest matched set has k controls, then inner is either a scalar or
a vector with k=length(inner).  If inner is a scalar, then
the same value of inner is used, regardless of the number of controls.
Otherwise, inner[1] is used with one control, inner[2] is used with
two controls, etc.  If inner is NULL, default values of
inner=c(.8,.8,.6,.4,0,0,0,...,0) are used.
}
  \item{trim}{
inner and trim together define the \eqn{\psi}-function for the M-statistic.
If the largest matched set has k controls, then trim is either a scalar or
a vector with k=length(trim).  If trim is a scalar, then
the same value of trim is used, regardless of the number of controls.
Otherwise, trim[1] is used with one control, trim[2] is used with
two controls, etc.  If trim is NULL, default values of
trim=c(3,3,...,3) are used.  For each i, 0 <= inner[i] < trim[i] < Inf.
}
  \item{lambda}{
Before applying the \eqn{\psi}-function to treated-minus-control differences,
the differences are scaled by dividing by the lambda quantile of all
within set absolute differences.  Typically, lambda = 1/2 for the median.
The value of lambda has no effect if trim=Inf and inner=0.  See Maritz (1979)
for the paired case and Rosenbaum (2007) for matched sets.

An error will result unless 0 < lambda < 1.
}
  \item{alpha}{
The coverage rate of the confidence interval is 1-\eqn{\alpha}.  If the bias in
treatment assignment is at most \eqn{\Gamma}, then the confidence interval covers
at rate 1-\eqn{\alpha}.
}
  \item{alternative}{
If alternative="greater" or alternative="less", the a one-sided confidence interval
is returned.  If alternative="twosided", a somewhat conservative two-sided confidence
interval is returned.  See the discussion of two-sided tests in the documentation
for sen().
  }
}
\details{
The confidence interval inverts the test provided by sen().  See the
documentation for sen() for more information.

The upper bound on the P-value is based on the separable approximation
described in Gastwirth, Krieger and Rosenbaum (2000); see also
Rosenbaum (2007, 2018).
}
\value{
\item{point.estimates }{For \eqn{\Gamma>1}, an interval of point estimates
is returned.  \eqn{\Gamma=1}, the interval is a point.}
\item{confidence.interval }{The confidence interval.}
}
\references{
Cox, D. R. (1977).  The role of signficance tests (with Discussion).
Scand. J. Statist. 4, 49-70.

Huber, P. (1981) Robust Statistics.  New York: John Wiley.
(M-estimates based on M-statistics.)

Li, Xinran and Rosenbaum, P. R. (2019) Maintaining high constant
design sensitivity in observational studies with matched sets of
varying sizes.  Manuscript.

Maritz, J. S. (1979). A note on exact robust confidence
intervals for location. Biometrika 66 163--166. (Introduces
exact permutation tests based on M-statistics by redefining
the scaling parameter.)

Rosenbaum, P. R. (2007). Sensitivity analysis for m-estimates, tests
and confidence intervals in matched observational studies.
Biometrics 63 456-64. (R package sensitivitymv) <doi:10.1111/j.1541-0420.2006.00717.x>

Rosenbaum, P. R. (2013). Impact of multiple matched controls on
design sensitivity in observational studies. Biometrics 69 118-127.
(Introduces inner trimming.) <doi:10.1111/j.1541-0420.2012.01821.x>

Rosenbaum, P. R. (2015).  Two R packages for sensitivity analysis
in observational studies.  Observational Studies, v. 1. (Free on-line.)

Rosenbaum, P. R. (2016) Using Scheffe projections for multiple outcomes
in an observational study of smoking and periondontal disease.  Annals of
Applied Statistics, 10, 1447-1471. <doi:10.1214/16-AOAS942>

Rosenbaum, P. R. (2018). Sensitivity analysis for stratified comparisons
in an observational study of the effect of smoking on homocysteine levels.
The Annals of Applied Statistics, 12(4), 2312-2334.
<doi:10.1214/18-AOAS1153>

}
\author{
Xinran Li and Paul R. Rosenbaum.
}

\note{The examples reproduce some results from Li and Rosenbaum (2019).}

\examples{
\dontrun{
# Reproduces results from Table 3 of Li and Rosenbaum (2019)
data(nh1and3)
attach(nh1and3)
senCI(homocysteine,z,mset,gamma=1)
senCI(homocysteine,z,mset,gamma=1.9)
senCI(homocysteine,z,mset,inner=0,gamma=1.9)
amplify(1.9,c(3,3.5,4))

# Relationships between confidence intervals and P-value bounds
senCI(homocysteine,z,mset,alternative="twosided",gamma=1.75)
sen(homocysteine,z,mset,alternative="less",tau=2.21721733,gamma=1.75)
senCI(homocysteine,z,mset,alternative="less",gamma=1.75)
sen(homocysteine,z,mset,alternative="less",tau=2.159342,gamma=1.75)
detach(nh1and3)
}
}

\keyword{ htest }
\keyword{ robust }

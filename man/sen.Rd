\name{sen}
\alias{sen}
\concept{Sensitivity analysis}
\concept{Observational study}
\concept{Causal inference}
\concept{Matching}
\concept{rbounds}
\title{
Sensitivity Analysis for a Matched Comparison in an Observational Study.
}
\description{
Each matched set contains one treated
individual and one or more controls.
Uses Huber's M-statistic as the basis for
the test, for instance, a mean.  Matched sets of different sizes
use different \eqn{\psi}-functions, creating what is called a fugue statistic.
Performs either a randomization
test (Gamma=1) or an analysis of sensitivity to departures from random
assignment (Gamma>1).  For confidence intervals, use function senCI().
The method is described in Li and Rosenbaum (2019); see also Rosenbaum (2007,2013).
}
\usage{
sen(y, z, mset, gamma = 1, inner = NULL, trim = NULL, lambda = 1/2,
     tau = 0, alternative = "greater")
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
  \item{tau}{
The null hypothesis asserts that the treatment has an additive effect, tau.
By default, tau=0, so by default the null hypothesis is Fisher's sharp null
hypothesis of no treatment effect.
}
  \item{alternative}{
If alternative="greater", the null hypothesis of a treatment effect of tau
is tested against the alternative of a treatment effect larger than tau.  If
alternative="less", the null hypothesis of a treatment effect of tau
is tested against the alternative of a treatment effect smaller than tau.
In particular, alternative="less" is equivalent to: (i) alternative="greater",
(ii) y replaced by -y, and (iii) tau replaced by -tau.  See the note for
discussion of two-sided sensitivity analyses.
  }
}
\details{
The novel element in the fugue package is the automatic use of different
\eqn{\psi}-functions for matched sets of different sizes.  These
\eqn{\psi}-functions have been selected to approximately equate the
design sensitivities in sets of unequal sizes when the errors are
Normal and the additive effect is half the standard deviation of
a matched pair difference; see Li and Rosenbaum (2019).  If you
disable this automatic feature by manually setting a single value for inner
and trim, then the results will agree with senm() in the R
package sensitivitymult.  For instance, using both sen() in the
fugue package and senm() in the sensitivitymult package will
yield the same deviate and P-value in:
data(nh1and3)
attach(nh1and3)
sen(homocysteine,z,mset,inner=0,gamma=1.9)
senm(homocysteine,z,mset,inner=0,trim=3,gamma=1.9)
Note that the sensitivitymult package is intended to
implement methods from Rosenbaum (2016,2019) that are
not implemented in the fugue package.

For the given \eqn{\Gamma}, sen() computes the upper bound on the 1-sided
P-value testing the null hypothesis
of an additive treatment effect tau against the alternative hypothesis of
a treatment effect larger than tau.  By default, sen() tests the null hypothesis of
no treatment effect against the alternative of a positive treatment effect.
The P-value is an approximate P-value
based on a Normal approximation to the null distribution; see Rosenbaum (2007).

Matched sets of unequal size are weighted using weights that would
be efficient in a randomization test under a simple model with
additive set and treatment effects and errors with constant variance;
see Rosenbaum (2007).

The upper bound on the P-value is based on the separable approximation
described in Gastwirth, Krieger and Rosenbaum (2000); see also
Rosenbaum (2007, 2018).
}
\value{
\item{pval }{The upper bound on the 1-sided P-value.}
\item{deviate }{The deviate that was compared to the Normal
distribution to produce pval.}
\item{statistic }{The value of the M-statistic.}
\item{expectation }{The maximum expectation of the
M-statistic for the given \eqn{\Gamma}.}
\item{variance }{The maximum variance of the
M-statistic among treatment assignments
that achieve the maximum expectation.
Part of the separable approximation.}
}
\references{
Cox, D. R. (1977).  The role of signficance tests (with Discussion).
Scand. J. Statist. 4, 49-70.

Huber, P. (1981) Robust Statistics.  New York: John Wiley.
(M-estimates based on M-statistics.)

Li, X. and Rosenbaum, P. R. (2019) Maintaining high constant
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

Rosenbaum, P. R. (2014). Weighted M-statistics with superior design
sensitivity in matched observational studies with multiple controls.
J. Am. Statist. Assoc. 109 1145-1158. (R package sensitivitymw)
<doi:10.1080/01621459.2013.879261>

Rosenbaum, P. R. (2015).  Two R packages for sensitivity analysis
in observational studies.  Observational Studies, v. 1. (Free on-line.)

Rosenbaum, P. R. (2016) Using Scheffe projections for multiple outcomes
in an observational study of smoking and periondontal disease.  Annals of
Applied Statistics, 10, 1447-1471. <doi:10.1214/16-AOAS942>

Rosenbaum, P. R. (2018). Sensitivity analysis for stratified comparisons
in an observational study of the effect of smoking on homocysteine levels.
The Annals of Applied Statistics, 12(4), 2312-2334.
<doi:10.1214/18-AOAS1153>

Rosenbaum, P. R. (2019). Combining planned and discovered comparisons
in observational studies.
Biostatistics, to appear. <doi.org/10.1093/biostatistics/kxy055>

}
\author{
Xinran Li and Paul R. Rosenbaum.
}
\note{
The function sen() performs 1-sided tests.  One approach
to a 2-sided, \eqn{\alpha}-level test does both 1-sided tests
at level \eqn{\alpha/2}, and rejects the null hypothesis if either
1-sided
test rejects.  Equivalently, a bound on the two sided
P-value is the smaller of 1 and twice the smaller of the two 1-sided
P-values.  This approach views a 2-sided test as two 1-sided tests
with a Bonferroni correction; see Cox (1977, Section 4.2).  In all
cases, this approach is a valid large sample test: a true
null hypothesis is falsely
rejected with probability at most \eqn{\alpha} if the bias in
treatment assignment is at most \eqn{\Gamma}; so, this procedure
is entirely safe to use.  For a randomization test, \eqn{\Gamma=1}, this
Bonferroni procedure is not typically conservative. For large \eqn{\Gamma},
this Bonferroni procedure tends to be somewhat conservative.

The examples reproduce some results from Li and Rosenbaum (2019).
}


\examples{
# Reproduces results from Table 3 of Li and Rosenbaum (2019)
data(nh1and3)
attach(nh1and3)
sen(homocysteine,z,mset,gamma=1)
sen(homocysteine,z,mset,gamma=1.9)
sen(homocysteine,z,mset,inner=0,gamma=1.9)
amplify(1.9,c(3,3.5,4))
detach(nh1and3)
}

\keyword{ htest }
\keyword{ robust }

\name{nh1and3}
\alias{nh1and3}
\concept{Causal inference}
\concept{Observational study}
\docType{data}
\title{
Smoking Matched Sets with 1 or 3 Controls
}
\description{
Data from NHANES 2005-2006 concerning homocysteine levels in daily
smokers (z=1) and never smokers (z=0), aged 20 and older.
}
\usage{data("nh1and3")}
\format{
  A data frame with 1370 observations consisting of 353 matched pairs and 166
  matched sets with 3 controls.
  \describe{
    \item{\code{SEQN}}{NHANES ID number}
    \item{\code{z}}{=1 for a daily smoker, =0 for a never smoker}
    \item{\code{mset}}{Matched set indicator, for 519 sets, 1, 2, ..., 519}
    \item{\code{homocysteine}}{Blood homocysteine level}
    \item{\code{cigsperday30}}{Cigarettes smoked per day}
    \item{\code{cotinine}}{Cotinine is a biomarker for exposure to nicotine}
    \item{\code{female}}{=1 for female, =0 for male}
    \item{\code{age}}{Age in years}
    \item{\code{black}}{=1 for black, =0 for other}
    \item{\code{education}}{NHANES 1-5 score.  3 is a high school degree.}
    \item{\code{povertyr}}{Ratio of family income to the poverty level, capped at 5.}
  }
}
\details{
Data from NHANES 2005-2006 concerning homocysteine levels in daily
smokers (z=1) and never smokers (z=0), aged 20 and older. Daily
smokers smoked every day for the last 30 days, smoking an average
of at least 10 cigarettes per day. Never smokers smoked fewer than
100 cigarettes in their lives, do not smoke now, and had no tobacco
use in the previous 5 days.
}
\source{
NHANES 2005-2006
}
\references{
Bazzano, L. A., He, J., Muntner, P., Vupputuri, S. and Whelton, P. K. (2003)
Relationship between cigarette smoking and novel risk
factors for cardiovascular disease in the United States.
Annals of Internal Medicine, 138, 891-897.

Pimentel, S. D., Small, D. S. and Rosenbaum, P. R. (2016) Constructed
second control groups and attenuation of unmeasured biases.
Journal of the American Statistical Association, 111, 1157-1167.
<doi.org/10.1080/01621459.2015.1076342>
}
\examples{
data(nh1and3)
attach(nh1and3)
table(table(nh1and3$mset))
par(mfrow=c(1,2))
boxplot(homocysteine[1:166]~z[1:166],ylim=c(0,70),main="1-1 match",
  ylab="homocysteine",names=c("Control","Smoker"))
boxplot(homocysteine[167:1370]~z[167:1370],ylim=c(0,70),
  main="1-3 match",ylab="homocysteine",names=c("Control","Smoke"))
detach(nh1and3)
}
\keyword{datasets}

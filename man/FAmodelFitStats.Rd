\name{FAmodelFitStats}
\alias{FAmodelFitStats}

\title{Calculate Summary Statistics with given FA Model Parameters}
\usage{
    FAmodelFitStats(B, Phi, omega, S, N)
}
\arguments{
    \item{B}{loadings.}
    \item{Phi}{cov. matrix of factors.}
    \item{omega}{vector of error variances}
    \item{S}{sample covariance matrix.}
    \item{N}{sample size.}
}

\value{a vector of various fit statistics.}
\description{Calculates various statistics with given Paramaters of an FA Model.}
\details{
  This function is used by FAfitStats and would not normally be called by a user.}
\seealso{
   \code{\link{FAfitStats}}
}
\author{Paul Gilbert and Erik Meijer}
\concept{TSFA}
\keyword{internal}


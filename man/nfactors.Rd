\name{nfactors}
\alias{nfactors}
\alias{nfactors.FAmodel}
%%\alias{nfactors.TSFmodel}
%%\alias{nfactors.TSFestModel}
\alias{nfactors.TSFfactors}
\alias{nfactors.EstEval}

\title{Extract the Number of Factors from an Object}
\usage{
    nfactors(x)
    \method{nfactors}{FAmodel}(x)
%%    \method{nfactors}{TSFmodel}(x)
%%    \method{nfactors}{TSFestModel}(x)
    \method{nfactors}{TSFfactors}(x)
    \method{nfactors}{EstEval}(x)
}
\arguments{
    \item{x}{an object.}
}
\description{
   Extract the number of factors from an object.
}

\value{an integer.}

\seealso{
   \code{\link{factors}},
   \code{\link{factorNames}},
   \code{\link{TSFmodel}},
}
\author{Paul Gilbert}
\concept{TSFA}
\keyword{ts}


\name{factorNames}
\alias{factorNames}
\alias{factorNames.TSFmodel}
\alias{factorNames.TSFestModel}
\alias{factorNames.TSFfactors}
\alias{factorNames.EstEval}
\alias{seriesNames.TSFestModel}

\title{Extract the Factors Names from an Object}
\usage{
    factorNames(x)
    \method{factorNames}{TSFmodel}(x)
    \method{factorNames}{TSFestModel}(x)
    \method{factorNames}{TSFfactors}(x)
    \method{factorNames}{EstEval}(x)
    \method{seriesNames}{TSFestModel}(x)
}
\arguments{
    \item{x}{an object.}
}
\description{
   Extract the factor (or series) names from an object.
}

\value{character vector of names.}

\seealso{
   \code{\link{factors}},
   \code{\link{nfactors}},
   \code{\link[tframe]{seriesNames}},
   \code{\link{TSFmodel}},
}
\author{Paul Gilbert}
\concept{TSFA}
\keyword{ts}


\name{summary.TSFmodel}
\alias{summary.TSFmodel}
\alias{print.summary.TSFmodel}
\alias{summary.FAmodel}
\alias{print.summary.FAmodel}
\alias{summary.TSFmodelEstEval}
\alias{print.summary.TSFmodelEstEval}

\title{summary.TSFmodel Method for Base Generic}
\usage{
    \method{summary}{TSFmodel}(object, ...)
    \method{summary}{FAmodel}(object, ...)
    \method{summary}{TSFmodelEstEval}(object, ...)
    \method{print}{summary.TSFmodel}(x, ...)
    \method{print}{summary.FAmodel}(x, ...)
    \method{print}{summary.TSFmodelEstEval}(x, digits = options()$digits, ...)
}
\arguments{
    \item{object}{an object to summarize.}
    \item{x}{an object to print.}
    \item{digits}{precision of printed numbers.}
    \item{...}{further arguments passed to other methods.}
}

\value{a summary object.}
\description{
   Summary method for object in \pkg{tsfa}, such as the object 
   returned by the estimation method \code{estTSF.ML}.
   See \code{\link{FAfitStats}} for details on the results from
   \code{summary.TSFmodel}.
}
\seealso{
   \code{\link{estTSF.ML}},
   \code{\link{FAfitStats}},
   \code{\link[base]{summary}}
}
\author{Paul Gilbert and Erik Meijer}
\concept{TSFA}
\keyword{ts}

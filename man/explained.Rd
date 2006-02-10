\name{explained.TSFmodel}
\alias{explained}
\alias{explained.TSFmodel}
\alias{explained.TSFestModel}

\title{Calculate Explained Portion of Data}
\usage{
    explained(object, ...)
    \method{explained}{TSFmodel}(object, names=object$names, ...)
    \method{explained}{TSFestModel}(object, ...)
}
\arguments{
    \item{object}{A TSFmodel or TSFestModel.}
    \item{names}{A vector of strings to use for the output series.}
    \item{...}{arguments passed to other methods.}
}
\description{
   Calculate portion of the data (indicators) explained by the factors.
}

\value{A time series matrix.}

\seealso{
   \code{\link{TSFmodel}},
   \code{\link{predict}},
   \code{\link{estTSF.ML}},
   \code{\link[dse1]{simulate}},
   \code{\link{tfplot.TSFmodel}}, 
}
\author{Paul Gilbert}
\concept{TSFA}
\keyword{ts}


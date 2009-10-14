\name{explained}
\alias{explained}
\alias{explained.TSFmodel}
\alias{explained.FAmodel}

\title{Calculate Explained Portion of Data}
\usage{
    explained(object, ...)
    \method{explained}{TSFmodel}(object, f=factors(object),
                  names=seriesNames(object), ...)
    \method{explained}{FAmodel}(object, f=factors(object),
                  names=dimnames(loadings(object))[[1]], ...)
}
\arguments{
    \item{object}{A TSFmodel or TSFestModel.}
    \item{f}{Factor values to use with the model.}
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
   \code{\link[dse]{simulate}},
   \code{\link{tfplot.TSFmodel}}, 
}
\author{Paul Gilbert}
\concept{TSFA}
\keyword{ts}


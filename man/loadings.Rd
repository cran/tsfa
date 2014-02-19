\name{loadings}
\alias{DstandardizedLoadings}
\alias{DstandardizedLoadings.TSFmodel}

\title{Extract the Loadings Matrix from an Object}
\usage{
    DstandardizedLoadings(x)
    \method{DstandardizedLoadings}{TSFmodel}(x)
}
\arguments{
    \item{x}{an object.}
}
\description{
   Extract the loadings matrix from an object. 
}

\value{a loadings matrix.}
\details{
   The loadings matrix in \code{TSFmodel} and 
   \code{TSFestModel} objects is similar to that described for the default,
   but calculated for a \acronym{TSFA} model. More details are provided in 
   \code{\link{estTSF.ML}}

   The \code{DstandardizedLoadings} matrix is the standardized loadings 
   for the differenced data.
}

\seealso{
   \code{stats:::\link[stats]{loadings}},
   \code{\link{factors}},
   \code{\link{factorNames}},
   \code{\link{estTSF.ML}},
   \code{\link{TSFmodel}},
}
\author{Paul Gilbert}
\keyword{internal}


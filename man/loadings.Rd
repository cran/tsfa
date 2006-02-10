\name{loadings}
\alias{loadings}
\alias{loadings.default}
\alias{loadings.TSFmodel}
\alias{loadings.TSFestModel}
%%\alias{standardizedLoadings}
%%\alias{standardizedLoadings.TSFestModel}
\alias{DstandardizedLoadings}
\alias{DstandardizedLoadings.TSFestModel}

\title{Extract the Loadings Matrix from an Object}
\usage{
    loadings(x)
    \method{loadings}{default}(x)
    \method{loadings}{TSFmodel}(x)
    \method{loadings}{TSFestModel}(x)
%%    standardizedLoadings(x)
%%    \method{standardizedLoadings}{TSFestModel}(x)
    DstandardizedLoadings(x)
    \method{DstandardizedLoadings}{TSFestModel}(x)
}
\arguments{
    \item{x}{an object.}
}
\description{
   Extract the loadings matrix from an object. 
}

\value{a loadings matrix.}
\details{
    \code{stats:::loadings} is defined as the default method for the 
   generic which replaces it. (See  \code{help(loadings, package="stats")} 
   for more details.)  The loadings matrix in \code{TSFmodel} and 
   \code{TSFestModel} objects is similar to that described for the default,
   but calculated for a \acronym{TSFA} model. More details are provided in 
   \code{\link{estTSF.ML}}
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


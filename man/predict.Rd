\name{predict}
\alias{predict.TSFmodel}
\alias{predict.TSFestModel}

\title{Predict Factor Scores from an Object.}
\usage{
    \method{predict}{TSFmodel}(object, 
           newdata = NULL, factorNames.=factorNames(object), ...)
    \method{predict}{TSFestModel}(object, 
           newdata = NULL, factorNames.=factorNames(object), ...)
}
\arguments{
    \item{object}{an object from which a matrix (predictor) can be extracted 
            to apply to the data.}
    \item{newdata}{data to which the predictor should be applied.}
    \item{factorNames.}{names to be given to the calculated predicted factor 
         score series.}
}
\description{
   Predict factor scores using the predictor from object.
   }

\value{Predicted factor scores series.}
\details{
    If \code{newdata} is not supplied then it is extacted from object if possible 
    (which is normally the data the model was estimated with), and otherwise 
    an error is indicated. The 
    predicted factor scores are given by \code{newdata \%*\% t(LB)}, where LB is the 
    factor score coefficient matrix extracted from object. 
    This is the Barlett factor score coefficient matrix if \code{TSFmodel}
    or \code{TSFestModel} objects were estimated with \code{estTSF.ML}.
    }

\seealso{
   \code{\link[stats]{predict}},
   \code{\link{factors}},
   \code{\link{factorNames}},
   \code{\link{TSFmodel}}
   }
\author{Paul Gilbert}
\keyword{internal}


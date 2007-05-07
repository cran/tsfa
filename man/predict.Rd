\name{predict}
\alias{predict.FAmodel}
\alias{predict.TSFmodel}
%%\alias{predict.TSFestModel}

\title{Predict Factor Scores from an Object.}
\usage{
    \method{predict}{FAmodel}(object, 
           data = NULL, factorNames.=factorNames(object), ...)
    \method{predict}{TSFmodel}(object, 
           data = object$data, factorNames.=factorNames(object), ...)
%%    \method{predict}{TSFestModel}(object, 
%%           data = NULL, factorNames.=factorNames(object), ...)
}
\arguments{
    \item{object}{an object from which a matrix (predictor) can be extracted 
            to apply to the data.}
    \item{data}{data to which the predictor should be applied.}
    \item{factorNames.}{names to be given to the calculated predicted factor 
         scores.}
    \item{...}{additional arguments currently unused.}
}
\description{
   Predict factor scores using the predictor from object.
   }

\value{Predicted factor scores.}
\details{
    If \code{data} is not supplied then it is extacted from object if possible 
    (which is normally the data the model was estimated with), and otherwise 
    an error is indicated. The 
    predicted factor scores are given by \code{data \%*\% t(LB)}, where LB is the 
    factor score predictor matrix extracted from object. 
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
\concept{TSFA}
\keyword{ts}


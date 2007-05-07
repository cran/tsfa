\name{FAmodel}
\alias{FAmodel}
\alias{FAmodel.default}
\alias{FAmodel.FAmodel}

\title{Construct a Factor Model}
\usage{
    FAmodel(obj, ...)
    \method{FAmodel}{default}(obj,  Omega=NULL, Phi=NULL, LB=NULL, LB.std=NULL,
   stats=NULL,  ...)
    \method{FAmodel}{FAmodel}(obj, ...)
}
\arguments{
    \item{obj}{The loadings matrix (\eqn{B}{B}) in the default (constructor) method.
	       In other methods, an object from which the model should be
	       extracted.}
    \item{Omega}{Covariance of the idiosyncratic term.}
    \item{Phi}{Covariance of the factors.}
    \item{LB}{Factor score predictor matrix.}
    \item{LB.std}{The standardized factor score predictor matrix.}
    \item{stats}{An optional list of statistics from model estimation.}
    \item{...}{arguments passed to other methods or stored in the object.}
}
\description{
   The default method constructs a FAmodel. Other methods extract a 
   FAmodel from an object.
}

\details{
  The default method is the constructor for \code{FAmodel} objects. Other
  methods extract a \code{FAmodel} object from other objects that contain one.
  The loadings must be supplied to the default method. Omega, 
  Phi, and LB are included when the object comes from an estimation method, but
  are not necessary when the object is being specified in order to simulate.
  The model is defined by

      \deqn{y_t  = B f_t + \varepsilon_t,}{%
            y(t) = B f(t) + eps(t),}
 
  where the factors \eqn{f_t}{f(t)} have covariance \eqn{\Phi}{Phi} 
  and \eqn{\varepsilon_t}{eps(t)} have covariance \eqn{\Omega}{Omega}.
  The loadings matrix \eqn{B}{B} is \eqn{M \times k}{M x k}, where \eqn{M}{M} 
  is the number of indicator variables
  (the number of indicators in \eqn{y}{y}) and \eqn{k}{k} is the number of 
  factors.
}
\value{A FAmodel.}
\examples{
  B <- t(matrix(c(0.9, 0.1,
		  0.8, 0.2,
		  0.7, 0.3,
                  0.5, 0.5, 
		  0.3, 0.7,
 		  0.1, 0.9), 2,6))

  z <- FAmodel(B)
  z
  loadings(z)
}
\seealso{
   \code{\link{TSFmodel}},
   \code{\link{estFAmodel}}
}
\author{Paul Gilbert}
\concept{TSFA}
\keyword{ts}

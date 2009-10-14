\name{TSFmodel}
\alias{TSFmodel}
\alias{TSFmodel.default}
\alias{TSFmodel.TSFmodel}
\alias{TSFmodel.FAmodel}
\title{Construct a Time Series Factor Model}
\usage{
    TSFmodel(obj, ...)
    \method{TSFmodel}{default}(obj, f=NULL, Omega = NULL, Phi=NULL, LB = NULL, 
	 positive.data=FALSE, names=NULL, ...)
    \method{TSFmodel}{TSFmodel}(obj, ...)
    \method{TSFmodel}{FAmodel}(obj, f=NULL, positive.data=FALSE, names=NULL, ...)
}
\arguments{
    \item{obj}{The loadings matrix (\eqn{B}{B}) in the default (constructor) method.
	       In other methods, an object from which the model should be
	       extracted.}
    \item{f}{matrix of factor series.}
    \item{Omega}{Covariance of the idiosyncratic term.}
    \item{Phi}{Covariance of the factors.}
    \item{LB}{Factor score coefficient matrix.}
%%    \item{break.points}{a list of dates, one shorter than the list of loadingss,
%%       indicating the start date of each successive loading matrix.}
    \item{positive.data}{logical indicating if any resulting negative values
       should be set to zero.}
    \item{names}{vector of strings indicating names to be given to output
       series.}
    \item{...}{arguments passed to other methods or stored in the object.}
}
\description{
   The default method constructs a TSFmodel. Other methods extract a 
   TSFmodel from an object.
}

\details{
  The default method is the constructor for \code{TSFmodel} objects. Other
  methods extract a \code{TSFmodel} object from other objects that contain one.
  The loadings and the factors must be supplied to the default method. Omega, 
  Phi, and LB are included when the object comes from an estimation method, but
  are not necessary when the object is being specified in order to simulate.
  The model is defined by

      \deqn{y_t  = B f_t + \varepsilon_t,}{%
            y(t) = B f(t) + eps(t),}
 
  where the factors \eqn{f_t}{f(t)} have covariance \eqn{\Phi}{Phi} 
  and \eqn{\varepsilon_t}{eps(t)} have covariance \eqn{\Omega}{Omega}.
  The loadings matrix \eqn{B}{B} is \eqn{M \times k}{M x k}, where \eqn{M}{M} 
  is the number of indicator variables
  (the number of series in \eqn{y}{y}) and \eqn{k}{k} is the number of 
  factor series.

  The estimation method \code{estTSF.ML} returns a \code{TSFmodel} as part of a 
  \code{TSFestModel} that has additional information about the estimation.
  
%%  WARNING: THIS PART OF THE CODE IS UNTESTED:
%%  In the more complicated case, where breakpoints are added,  a list of
%%  loadings are used with the change starting at each date in the list 
%%  \code{break.points}. The first loadings matrix is used beginning at 
%%  start(factors) and ending one period before
%%  \code{break.points[[1]]}.
}
\value{A TSFmodel.}
\examples{
  f <- matrix(c(2+sin(pi/100:1),5+3*sin(2*pi/5*(100:1))),100,2)
  B <- t(matrix(c(0.9, 0.1,
		  0.8, 0.2,
		  0.7, 0.3,
                  0.5, 0.5, 
		  0.3, 0.7,
 		  0.1, 0.9), 2,6))

  z <- TSFmodel(B, f=f)
  tfplot(z)
}
\seealso{
   \code{\link{simulate.TSFmodel}},
   \code{\link[dse]{simulate}},
   \code{\link{estTSF.ML}}
}
\author{Paul Gilbert}
\concept{TSFA}
\keyword{ts}

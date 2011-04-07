\name{simulate.TSFmodel}
\alias{simulate.TSFmodel}

\title{Simulate a Time Series Factor Model}
\usage{
    \method{simulate}{TSFmodel}(model, f=factors(model), Cov=model$Omega,
       sd=NULL, noise=NULL, rng=NULL, noise.model=NULL, ...)
}
\arguments{
    \item{model}{A TSFmodel.}
    \item{f}{Factors to use with the model.}
    \item{Cov}{covariance of the idiosyncratic term.}
    \item{sd}{see \code{\link[dse]{makeTSnoise}}.}
    \item{noise}{see \code{\link[dse]{makeTSnoise}}.}
    \item{rng}{see \code{\link[dse]{makeTSnoise}}.}
    \item{noise.model}{see \code{\link[dse]{makeTSnoise}}.}
    \item{...}{arguments passed to other methods.}
}
\description{
   Simulate a TSFmodel to generate time series data (indicators) using 
   factors and loadings from the model.
}
\details{
  \code{simulate.TSFmodel} generates artifical data (indicators or measures) 
  with a given \code{TSFmodel} (which has factors and loadings). The \code{obj} 
  should be a \code{TSFmodel}. This might be a model constructed with
  \code{\link{TSFmodel}} or as returned by \code{\link{estTSF.ML}}.
  
  The number of factor series is determined by the number of columns in 
  the time series matrix \code{f} (the factors in the model object). 
  This must also be the number of columns
  in the loadings matrix \eqn{B}{B} (in the model object). 
  The number of rows in the loadings matrix determines the number of 
  indicator series generated (the number of columns in the matrix result). 
  The number of rows in the time series factor matrix determines the number of
  time observations (periods) in the indicator series generated, that is, the
  number of rows in the matrix result.
  
  \code{simulate} passes \code{Cov}, \code{sd}, \code{noise}, \code{rng}, 
  and \code{noise.model} to \code{\link[dse]{makeTSnoise}} to generate 
  the random idiosyncratic term \eqn{\varepsilon_t}{eps(t)}, 
  which will have the same dimension as the generated indicator series that
  are returned. 
  \eqn{\varepsilon_t}{eps(t)} will have random distribution determined by other
  arguments passed to \code{\link[dse]{makeTSnoise}}. 
  Note that the covariance of the generated indicator series \eqn{y_t}{y(t)}
  is also influenced by the covariance of the factors \eqn{f}{f}.
  
  The calculation to give the generated artificial time series indicator
  data matrix \eqn{y}{y} is 

      \deqn{y_t  = B f_t + \varepsilon_t.}{%
            y(t) = B f(t) + eps(t).}

  \code{simulate.TSFmodel} can use a \code{\link{TSFmodel}} that has only 
  \code{B}  and \code{f} specified, but in this case one of \code{Cov}, 
  \code{sd}, \code{noise}, or \code{noise.model} must be specified as the
  default \code{Omega} from the model is not available.
}
\value{A time series matrix.}
\examples{
  f <- matrix(c(2+sin(pi/100:1),5+3*sin(2*pi/5*(100:1))),100,2)
  B <- t(matrix(c(0.9, 0.1,
		  0.8, 0.2,
		  0.7, 0.3,
                  0.5, 0.5, 
		  0.3, 0.7,
 		  0.1, 0.9), 2,6))

  z <- simulate(TSFmodel(B, f=f), sd=0.01)
  tfplot(z)
}

\seealso{
   \code{\link{TSFmodel}},
   \code{\link{estTSF.ML}},
   \code{\link[dse]{simulate}},
   \code{\link{tfplot.TSFmodel}}, 
   \code{\link{explained.TSFmodel}}
}
\author{Paul Gilbert}
\concept{TSFA}
\keyword{ts}


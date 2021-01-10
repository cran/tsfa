\name{estFAmodel}
\alias{estFAmodel}

\title{Estimate a Factor Model}
\usage{
     estFAmodel(Sigma, p, n.obs=NA,
                est="factanal", 
		estArgs=list(scores="none", control=list(opt=list(maxit=10000))),
		rotation=if(p==1) "none" else "quartimin", rotationArgs=NULL,
		GPFargs=list(Tmat=diag(p), normalize=TRUE, eps=1e-5, maxit=1000),
		BpermuteTarget=NULL,
                factorNames=paste("Factor", seq(p)),
                indicatorNames=NULL)
}
\arguments{
    \item{Sigma}{covariance of the data matrix.}
    \item{n.obs}{integer indication number of observations in the dataset.}
    \item{p}{integer indication number of factors to estimate.}
    \item{est}{name of the estimation function.}
    \item{estArgs}{list of aarguments passed to the estimation function.}
   \item{rotation}{character vector indicating the factor 
         rotation method (see \pkg{GPArotation} for many options).}
    \item{rotationArgs}{list of arguments passed to the rotation method, 
          specifying arguments for the rotation criteria. See \code{GPFoblq}.}
    \item{GPFargs}{list of arguments passed to \code{GPFoblq} or \code{GPForth}
           for rotation optimization}
    \item{BpermuteTarget}{matrix of loadings. If supplied, this is used to 
             permute the order of estimated factors and change signs. (It is
	     for comparison with other results.}
    \item{factorNames}{vector of strings indicating names of factor series.}
    \item{indicatorNames}{vector of strings indicating names of indicator series.}
}

\description{
   Estimate an \code{FAmodel}. 
}

\value{A \code{FAmodel} object (see details). 
    
    }

\details{
  The default \code{est} method and \code{quartimin} rotation give parameters 
  using standard 
  (quasi) ML factor analysis (on the correlation matrix and then scaled back).
  The function \code{factanal} with no rotation is used to find the initial
  (orthogonal) solution. Rotation is then done 
  (by default with quartimin using \code{GPFoblq} optimization).
  \code{factanal} always uses the correlation matrix, so standardizing does 
  not affect the solution. 
  
  If \code{rotation} is \code{"none"} the result of the \code{factanal} 
  estimation is not rotated. In this case, to avoid confusion with a rotated
  solution, the factor covariance matrix \code{Phi} is returned as \code{NULL}.
  Another possibility for its value would be the identity matrix, but this is
  not calculated so \code{NULL} avoids confusion.
  
  The arguments \code{rotation}, \code{rotationArgs} are used for rotation.
  The quartimin default uses GPArotation and its default
  \code{normalize=TRUE}, \code{eps=1e-5}, \code{maxit=1000}, and \code{Tmat=I} 
  are passed through the rotation method to \code{\link[GPArotation]{GPFoblq}}. 

  The estimated loadings, Bartlett predictor matrix, etc.,
  are put in the returned \code{FAmodel} (see below).
  The Bartlett factor score coefficient matrix can be calculated as 

  \deqn{(B' \Omega^{-1} B)^{-1} B' \Omega^{-1} x}{%
        (B' Omega exp(-1) B) exp(-1) B' Omega exp(-1) x}

  or equivalently as

   \deqn{(B' \Sigma^{-1} B)^{-1} B' \Sigma^{-1} x,}{%
         (B' Sigma exp(-1) B) exp(-1) B' Sigma exp(-1) x}

  The first is simpler because \eqn{\Omega}{Omega} is diagonal, but breaks down 
  with a Heywood case, because  \eqn{\Omega}{Omega} is then singular (one or 
  more of its diagonal elements are zero). The second only requires 
  nonsingularity of \eqn{\Sigma}{Sigma}. Typically, \eqn{\Sigma}{Sigma} is not singular 
  even if \eqn{\Omega}{Omega} is singular.
  \eqn{\Sigma}{Sigma} is calculated  from \eqn{B  \Phi B' + \Omega}{B  Phi B' + Omega}, 
  where \eqn{B, \Phi,}{B, Phi,} and \eqn{\Omega}{Omega} are the
  estimated values returned from \code{factanal} and rotated. 
  The data covariance could also be used for \eqn{\Sigma}{Sigma}. 
  (It returns the same result with this estimation method.)
  
  The returned \code{FAmodel} object is a list containing
  \describe{
        \item{loadings}{the estimated loadings matrix.}
        \item{Omega}{the covariance of the idiosyncratic component (residuals).} 
        \item{Phi}{the covariance of the factors.}
        \item{LB}{the Bartlett predictor matrix.} 
        \item{LB.std}{the standardized Bartlett predictor matrix.}  
        \item{estConverged}{a logical indicating if estimation converged.}
        \item{rotationConverged}{a logical indicating if rotation converged.}
        \item{orthogonal}{a logical indicating if the rotation is orthogonal.}
        \item{uniquenesses}{the uniquenesses.}
        \item{call}{thearguments of the function call.} 
        }
  }
\references{
  Gilbert, Paul D. and Meijer, Erik (2005)
  Time Series Factor Analaysis with an Application to Measuring Money.
  Research Report 05F10, University of Groningen, SOM Research School.
  Available from \url{https://hdl.handle.net/11370/d7d4ea3d-af1d-487a-b9b6-c0816994ef5a}.
  }
\examples{
  data("WansbeekMeijer", package="GPArotation")
  fa.unrotated  <-  estFAmodel(NetherlandsTV, 2, n.obs=2150, rotation="none" )
  fa.varimax <- estFAmodel(NetherlandsTV, 2, n.obs=2150, rotation="Varimax" )
  fa.eiv     <- estFAmodel(NetherlandsTV, 2, n.obs=2150, rotation="eiv" )
  fa.oblimin <- estFAmodel(NetherlandsTV, 2, n.obs=2150, rotation="oblimin" )
  
  cbind(loadings(fa.unrotated), loadings(fa.varimax), loadings(fa.oblimin), loadings(fa.eiv))
  }

\seealso{
   \code{\link{estTSF.ML}},
   \code{\link[GPArotation]{rotations}},
   \code{\link{factanal}}
}
\author{Paul Gilbert and Erik Meijer}
\concept{TSFA}
\keyword{ts}


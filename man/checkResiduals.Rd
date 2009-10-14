\name{checkResiduals.TSFmodel}
\alias{checkResiduals.TSFmodel}
\title{Check Time Series Idiosyncratic Component}
\usage{
    \method{checkResiduals}{TSFmodel}(obj, data=obj$data, diff.=TRUE, ...)
}
\arguments{
    \item{obj}{TSFmodel object for which the idiosyncratic component should be
       examined (as if it were a residual).}
    \item{data}{data from which the idiosyncratic component should be
       calculated.}
    \item{diff.}{logical indicating if data and explained should be differenced.}
    \item{...}{arguments to be passed to checkResiduals default methods.}
}
\description{
   The data is subtracted from the explained data (after differencing if 
   \code{diff} is \code{TRUE}, the default) and the result is treated as a
   residual. Its covariance, the sum of the diagonal elements of the covariance,
   and the sum of the off-diagonal elements of the covariance are printed.
   The residual is then passed to the default method for \code{checkResiduals}
   which produces several diagonistic plots and (invisibly) returns statistics.
   See \code{\link[dse]{checkResiduals}} for more details. 
   Calculation of partial autocorrelations can be problematic.
   
   Some care should be taken interpreting the results. Factor estimation does
   not minimize residuals, it extracts common factors.
}


\examples{
  if (require("CDNmoney")){
    data("CanadianMoneyData.asof.28Jan2005", package="CDNmoney")
    data("CanadianCreditData.asof.28Jan2005", package="CDNmoney")
 
    z <- tframed(tbind(
  	MB2001,
  	MB486 + MB452 + MB453 ,
  	NonbankCheq,
  	MB472 + MB473 + MB487p,
  	MB475,
  	NonbankNonCheq + MB454 + NonbankTerm + MB2046 + MB2047 + MB2048 +
  	MB2057 + MB2058 + MB482),
  	names=c("currency", "personal cheq.", "NonbankCheq",
  	"N-P demand & notice", "N-P term", "Investment" )
      )
 
    z <- tfwindow(tbind (z, ConsumerCredit, ResidentialMortgage,
  			    ShortTermBusinessCredit, OtherBusinessCredit),
  	 start=c(1981,11), end=c(2004,11))
 
    cpi <- 100 * M1total / M1real
    popm <- M1total / M1PerCapita
    scale <- tfwindow(1e8 /(popm * cpi), tf=tframe(z))

    MBandCredit <- sweep(z, 1, scale, "*")
    c4withML  <- estTSF.ML(MBandCredit, 4)

    checkResiduals(c4withML, pac=FALSE)
  }
}

\seealso{
   \code{\link[dse]{checkResiduals}},
   \code{\link{TSFmodel}},
   \code{\link{estTSF.ML}}
}
\author{Paul Gilbert}
\concept{TSFA}
\keyword{ts}

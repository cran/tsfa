\name{tsfa-package}
\alias{tsfa-package}
\alias{tsfa.Intro}
\docType{package}

\title{Time Series Factor Analysis (TSFA)}

\description{TSFA extends standard factor analysis (FA) to 
time series data. Rotations methods can be applied as in FA.
A dynamic model of the factors is not assumed, but could be 
estimated separately using the extracted factors.
}

\details{
\tabular{ll}{
Package: \tab tsfa\cr
Depends: \tab R (>= 2.0.0), GPArotation, setRNG (>= 2004.4-1), tframe (>=
2006.1-1),\cr
         \tab dse (>= 2006.1-1), EvalEst (>= 2006.1-1)\cr
Suggests: \tab CDNmoney\cr
License: \tab GPL Version 2.\cr
URL: \tab http://tsanalysis.r-forge.r-project.org/\cr
}

The main functions are:
\preformatted{
DstandardizedLoadings   Extract standardized loadings from an object
loadings                Extractloadings from an object
estTSF.ML               Estimate a time series factor model
factors                 Extract time series factors from an object
FAmodelFitStats         Various fit statistics.
simulate                Simulate a time series factor model
summary                 Summary methods for \pkg{tsfa} objects
tfplot                  Plot methods for \pkg{tsfa} objects
TSFmodel                Construct a time series factor model
}

An overview of how to use the package is available in the vignette
\code{tsfa} (source, pdf).
}

\author{
Paul Gilbert <pgilbert.ttv9z@ncf.ca> and
Erik Meijer  <meijer@rand.org>

Maintainer: Paul Gilbert <pgilbert.ttv9z@ncf.ca>  
}
\references{
  Gilbert, Paul D. and Meijer, Erik (2005)
  Time Series Factor Analaysis with an Application to Measuring Money.
  Research Report 05F10, University of Groningen, SOM Research School.
  Available from \url{https://hdl.handle.net/11370/d7d4ea3d-af1d-487a-b9b6-c0816994ef5a}.

  Gilbert, Paul D. and Meijer, Erik (2006) Money and Credit Factors.
  Bank of Canada Working Paper 2006-3, available at
\url{https://www.bankofcanada.ca/2006/03/working-paper-2006-3/}.
}
\keyword{ package }
\keyword{ ts }

\seealso{
   \code{\link{estTSF.ML}}, 
   \code{\link[GPArotation:GPArotation-package]{GPArotation}}, 
   \code{\link[tframe:tframe-package]{tframe}}, 
   \code{\link[dse:dse-package]{dse}} 
}

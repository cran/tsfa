\name{permusign}
\alias{permusign}


\title{Internal Utility to Permute the Loadings Matrix.}
\usage{
    permusign(B, Btarget, Phi=NULL)
}
\arguments{
    \item{B}{proposed loadings matrix.}
    \item{Btarget}{target loadings matrix.}
    \item{Phi}{proposed Phi matrix.}
}
\description{
   Internal utility to permute the loadings matrix.
}

\value{list with a permuted and sign changed loadings matrix and the
    corresponding \code{Phi} matrix.}

\seealso{
   \code{\link{factors}},
   \code{\link{factorNames}},
   \code{\link{TSFmodel}},
}
\author{Paul Gilbert and Erik Meijer}
\keyword{internal}



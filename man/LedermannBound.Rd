\name{LedermannBound}
\alias{LedermannBound}

\title{Ledermann Bound for Number of Indicators}
\usage{
    LedermannBound(M)
}
\arguments{
    \item{M}{an integer indicating the number of indicator variables or a matrix
    of data, in which case ncol(M) is used as the number of indicator variables.}
}

\value{The Ledermann bound, a positive real number.}
\description{
   The Ledermann bound is given by the solution \eqn{k}{k} for
   \eqn{(M - k)^2 \geq M + k}{(M - k)^2 >= M + k}, where \eqn{M}{M} is 
   the number of
   indicator variables. The maximum possible number of factors is the largest
   integer smaller than or equal \eqn{k}{k}.
}
\seealso{
   \code{\link{FAfitStats}}
}
\references{
  Tom Wansbeek and Erik Meijer (2000) \emph{Measurement 
   Error and Latent Variables in Econometrics}, Amsterdam: North-Holland. 
   (note p169.)
  }

\author{Paul Gilbert and Erik Meijer}
\concept{TSFA}
\keyword{ts}


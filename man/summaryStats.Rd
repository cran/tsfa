
\name{summaryStats}
\alias{summaryStats}
\alias{summaryStats.TSFmodelEstEval}

\title{Summary Statistics Calculations}
\usage{
    summaryStats(object, ...)
    \method{summaryStats}{TSFmodelEstEval}(object, ...)
}
\arguments{
    \item{object}{a TSFestModel object to summarize.}
    \item{...}{further arguments passed to other methods.}
}

\value{a list passed of statistics.}
\description{
   Calculates various statistics from a TSFmodelEstEval object returned 
   by EstEval. This
   function is for use by the summary and tfplot methods and would not typically
   be called by a user.
}
\seealso{
   \code{\link[EvalEst]{EstEval}},
   \code{\link{summary.TSFmodelEstEval}},
   \code{\link{tfplot.TSFmodelEstEval}}
}
\author{Paul Gilbert and Erik Meijer}
\concept{TSFA}
\keyword{ts}
\keyword{internal}


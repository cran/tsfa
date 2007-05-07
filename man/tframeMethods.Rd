\name{tframeMethods}
\alias{tframe.TSFmodel}
\alias{tfplot.TSFmodel}
\alias{tfplot.TSFfactors}
\alias{tfplot.TSFexplained}
\alias{tfplot.TSFmodelEstEval}

\alias{diff.TSFmodel}
\alias{diff.TSFexplained}
\alias{diff.TSFfactors}

\title{Time Series Factor Methods for tframe Generics}
\usage{
   \method{tframe}{TSFmodel}(x) 

   \method{tfplot}{TSFmodel}(x, ..., tf=tfspan(x, ...), start=tfstart(tf), end=tfend(tf), 
		 series = seq(nfactors(x)), 
		 Title = "Model factors", 
		 lty = 1:5, lwd = 1, pch = NULL, col = 1:6, cex = NULL,
		 xlab = NULL, ylab = factorNames(x), xlim = NULL, ylim = NULL, 
		 graphs.per.page = 5,
                 par=NULL, mar = par()$mar, reset.screen = TRUE)
   \method{tfplot}{TSFfactors}(x,..., tf=tfspan(x, ...), start=tfstart(tf), end=tfend(tf), 
		 series=seq(nfactors(x)),
		 Title="Estimated factors (dashed) and true (solid)", 
		 lty = c("dashed", "solid"), lwd = 1, pch = NULL, col = 1:6, cex = NULL,
		 xlab=NULL, ylab=factorNames(x), xlim = NULL, ylim = NULL, 
		 graphs.per.page=5, par=NULL, mar=par()$mar, reset.screen=TRUE)
   \method{tfplot}{TSFexplained}(x,..., tf=tfspan(x, ...), start=tfstart(tf), end=tfend(tf), 
		 series=seq(nseries(x)),
		 Title="Explained (dashed) and actual data (solid)", 
		 lty = c("dashed", "solid"), lwd = 1, pch = NULL, col = 1:6, cex = NULL,
		 xlab=NULL, 
		 ylab=seriesNames(x), 
		 xlim = NULL, ylim = NULL,
		 graphs.per.page=5, par=NULL, mar=par()$mar, reset.screen=TRUE)
   \method{tfplot}{TSFmodelEstEval}(x, ...,  tf=NULL, start=tfstart(tf), end=tfend(tf), 
		 series=seq(nseries(factors(x))),
		 Title="Monte Carlo Results", 
		 lty = c("solid", "dotdash", "dashed",  "dashed"), lwd = 1, pch = NULL, 
		 col = c("black", "red", "red", "red"), cex = NULL,
		 xlab=NULL, 
		 ylab=seriesNames(factors(x$truth)), 
		 xlim = NULL, ylim = NULL,
		 graphs.per.page=5, par=NULL, mar=par()$mar, reset.screen=TRUE,
		 diff.=FALSE,  percentChange.=FALSE, 
		 PCcentered.=FALSE, summary.=TRUE)

    \method{diff}{TSFmodel}(x, ...)
    \method{diff}{TSFexplained}(x, ...)
    \method{diff}{TSFfactors}(x, ...)
    \method{diff}{factorsEstEval}(x, ...)
}
\arguments{
    \item{x}{an object.}
    \item{x}{a TSFmodel, TSFestModel, TSFexplained, or TSFfactors object 
        for plotting or differencing.}
    \item{diff.}{logical indicating if differenced data should be plotted.}
    \item{percentChange.}{logical indicating if percent change data 
          should be plotted.}
    \item{PCcentered.}{logical indicating if centered percent change data 
          should be plotted.}
    \item{summary.}{logical indicating if mean and 1 SD bounds should be 
         plotted in place of all estimates.}
    \item{tf}{See generic tfplot method}
    \item{start}{See generic tfplot method}
    \item{end}{See generic tfplot method}
    \item{series}{See generic tfplot method}
    \item{Title}{string to use for title of factors plot (but see details).}
    \item{lty}{See generic tfplot method}
    \item{lwd}{See generic tfplot method}
    \item{pch}{See generic tfplot method}
    \item{col}{See generic tfplot method}
    \item{cex}{See generic tfplot method}
    \item{xlab}{See generic tfplot method}
    \item{ylab}{See generic tfplot method}
    \item{xlim}{See generic tfplot method}
    \item{ylim}{See generic tfplot method}
    \item{graphs.per.page}{See generic tfplot method}
    \item{par}{See generic tfplot method}
    \item{mar}{See generic tfplot method}
    \item{reset.screen}{See generic tfplot method}
    \item{...}{other objects to plot (currently unused).}
}
\description{Plot or difference objects. See the generic descriptions.}

\value{diff returns an object in which the time series data has been 
   differenced. tfplot returns an invisible value but is executed mainly
   for the side-effect (plot).}
\details{
    The \code{Title} is not put on the plot if the global option
    PlotTitles is FALSE. This can be set with \code{options(PlotTitles=FALSE)}. 
    This provides a convenient  mechanism to omit all titles when the title
    may be added separately (e.g. in Latex).
}

\seealso{
   \code{\link{TSFmodel}},
   \code{\link{estTSF.ML}},
   \code{\link{simulate.TSFmodel}},
   \code{\link[tframe]{tfplot}},
   \code{\link{diff}},
   \code{\link{factors}},
   \code{\link{explained}},
   \code{\link{factorNames}},
   \code{\link{TSFmodel}}
}
\author{Paul Gilbert}
\concept{TSFA}
\keyword{internal}

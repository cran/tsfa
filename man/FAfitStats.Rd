\name{FAfitStats}
\alias{FAfitStats}
\alias{FAfitStats.default}
\alias{FAfitStats.TSFmodel}

\encoding{latin1}
%%\encoding{UTF-8}
\title{Summary Statistics for a TSFA Models}
\usage{
    FAfitStats(object, ...)
    \method{FAfitStats}{default}(object, diff. = TRUE,
               N=(nrow(object) - diff.), 
	       control=list(lower = 0.0001, opt=list(maxit=1000)), ...)
    \method{FAfitStats}{TSFmodel}(object, diff. = TRUE,
               N=(nrow(object$data) - diff.), ...)
}
\arguments{
    \item{object}{a time series matrix or TSFestModel.}
    \item{diff.}{logical indicating if data should be differenced.}
    \item{N}{sample size.}
    \item{control}{a list of arguments passed to \code{factanal}.}
    \item{...}{further arguments passed to other methods.}
}

\value{a vector or list of various fit statistics. See details.}
\description{
   \code{FAfitStats} calculates various statistics for a \code{TSFestModel} 
   or all possible (unrotated factanal) models for a data matrix. This
   function is also used by the summary method for a \code{TSFestModel}.
}
\details{
  In the case of the method for a \code{TSFmodel} the model parameters are 
  extracted from the model and the 
  result is a vector of various fit statistics (see below). 
  (Calculations are done by the internal function \code{FAmodelFitStats}.) 

  Most of these statistics are described in 
  \cite{Wansbeek and Meijer} (2000, WM below).
  The sample size \eqn{N} is used in the calculation of these statistics.
  The default is the number of number of observations, as in WM. That is, the
  number of rows in the data matrix, minus one if the data is differenced. Many
  authors use \eqn{N - 1}{N - 1}, which would be \eqn{N-2}{N-2} 
  if the data is differenced.    
  The exact calculations can be determined by examining the code: 
  \code{print(tsfa:::FAmodelFitStats)}. The vector of statistics is:
  \describe{
    \item{chisq}{Chi-square statistic (see, for example, WM p298).}
    \item{df}{degrees of freedom, which takes the rotational freedom 
  	 into account (WM p169).}
    \item{pval}{p-value}
    \item{delta}{delta}
    \item{RMSEA}{Root mean square error of approximation (WM p309).}
    \item{RNI}{Relative noncentrality index (WM p307).}
    \item{CFI}{Comparative fit index (WM p307).}
%%    \item{NFI}{Normed fit index (WM, p306).}
%%    \item{TLI}{Tucker-Lewis Index, also called Nonnormed fit index (WM p308).}
%%    \item{BL86}{Bollen (1986)}
%%    \item{BL89}{Bollen (1989)}
    \item{MCI}{McDonald's centrality index.}
    \item{GFI}{Goodness of fit index ( \enc{Jöreskog}{Joreskog} and
    \enc{Sörbom}{Sorbom}, 1981, 1986, WM p305).}
    \item{AGFI}{Adjusted GFI (\enc{Jöreskog}{Joreskog} and
    \enc{Sörbom}{Sorbom}, 1981, 1986).}
%%   Hoelter's (1983) Critical N; note that the significance level
%%    alpha must be given. E.g.: alpha 0.05
    \item{AIC}{Akaike's information criterion (WM p309).}
    \item{CAIC}{Consistent AIC(WM p310).}
    \item{SIC}{Schwarz's Bayesian information criterion.}
    \item{CAK}{Cudeck & Browne's rescaled AIC.}
    \item{CK}{Cudeck & Browne's cross-validation index.}
  }

  The information criteria account for rotational freedom.
  Some of these goodness of fit statistics should be used with caution, because
  they are not yet based on sound statistical theory. Future versions of tsfa
  will probably provide improved versions of these goodness-of-fit statistics.

  In the case of the default method, which expects a matrix of data with columns
  for each indicator series, models are calculated with \code{factanal} for 
  factors up to the Ledermann bound. No rotation is needed, since rotation does 
  not affect the fit statistics. Values for the saturated model are also
  appended to facilitate a sequential comparison. 
  
  If \code{factanal} does not obtain a satisfactory solution  it may produce an
  error "unable to optimize from these starting value(s)." This can sometimes be
  fixed by increasing the \code{opt, maxit} value in the \code{control} list.
  
  The result for the default method is a list with elements
  \describe{
    \item{fitStats}{a matrix with rows as for a single model above, 
  	and a column for each possible number of factors.}
    \item{seqfitStats}{a matrix with rows \code{chisq}, \code{df}, and
  	 \code{pval}, and columns indicating the comparative fit for an additional
  	 factor starting with the null (zero factor) model. 
  	 (See also independence model, WM, p305)
  	 }
    }

  The largest model can correspond to the saturated model, but will not if the 
  Ledermann bound is not an integer, or even in the case of an integer bound
  but implicit contraints resulting in a Heywood case (see Dijkstra, 1992).
  In these situations it might make sense to remove the model
  corresponding to the largest integer, and make the last sequential comparison
  between the second to largest integer and the saturated solution. The code
  does not do this automatically. 
  }

\seealso{
   \code{\link{FAmodelFitStats}},
   \code{\link{summary}},
   \code{\link{summary.TSFmodel}},
   \code{\link{summaryStats}},
   \code{\link{LedermannBound}}
}
\references{
  Dijkstra, T. K. (1992) On Statistical Inference with Parameter Estimates
    on the Boundary of the Parameter Space, 
    \emph{British Journal of Mathematical and Statistical
    Psychology}, \bold{45}, 289--309.

  Hu, L.-t., and Bentler, P. (1995) Evaluating model fit. In R. H. Hoyle
  (Ed.), \emph{Structural equation modeling: Concepts, issues, and
  applications} (pp. 76--99). Thousand Oaks, CA: Sage.   

  \enc{Jöreskog}{Joreskog}, K. G., and \enc{Sörbom}{Sorbom}, D. (1981) 
  \emph{LISREL V user's guide}. Chicago: National  Educational Resources.

  \enc{Jöreskog}{Joreskog}, K. G., and \enc{Sörbom}{Sorbom}, D. (1986) 
  LISREL VI: Analysis of linear structural relationships by maximum likelihood, instrumental
  variables, and least squares methods (User's Guide, 4th ed.).
  Mooresville, IN: Scientific Software.

  Ogasawara, Haruhiko. (2001). Approximations to the Distributions of Fit
  Indexes for Misspecified Structural Equation Models.
  \emph{Structural Equation Modeling}, \bold{8}, 556--574.

  Wansbeek, Tom and Meijer, Erik (2000) \emph{Measurement Error and 
   Latent Variables in Econometrics}, Amsterdam: North-Holland.
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

    FAfitStats(MBandCredit)

    c4withML  <- estTSF.ML(MBandCredit, 4)
    FAfitStats(c4withML)
  }
}
\author{Paul Gilbert and Erik Meijer}
\concept{TSFA}
\keyword{ts}



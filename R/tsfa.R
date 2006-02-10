###################################################

#   Time Series Factor Analysis (Latent Variables)

###################################################

# I don't think this exists really works for namespaces
# The if also seems to cause some problems for codoc
#if (!exists("loadings.default", mode="function")){
  loadings.default  <- stats::loadings
  loadings <- function(x)UseMethod("loadings")
#  }

loadings.TSFmodel <- function(x)x$B
loadings.TSFestModel <- function(x) loadings(TSFmodel(x))

DstandardizedLoadings <- function(x)UseMethod("DstandardizedLoadings")
DstandardizedLoadings.TSFestModel <- function(x){
    r <- diag(1/sqrt(diag(cov(diff(x$data))))) %*% loadings(x)
    dimnames(r) <- dimnames(loadings(x))
    r
    }


# standardizedLoadings <- function(x)UseMethod("standardizedLoadings")
# standardizedLoadings.TSFestModel <- function(x){
#   To transform everything back to
#   the undifferenced scales and then standardize, then somehow a
#   pseudo-Phi and pseudo-Omega (called \Gamma_t and \Psi_t in the TSFA
#   paper) must be computed first.    
#     r <- diag(1/sqrt(diag(cov(x$data)))) %*% loadings(x)
#     dimnames(r) <- dimnames(loadings(x))
#     r
#     }

TSFmodel <- function(obj, ...)UseMethod("TSFmodel")
TSFmodel.TSFmodel <- function(obj, ...) obj #extractor
TSFmodel.TSFestModel <- function(obj, ...) obj$model  #extractor

TSFmodel.default <- function(obj, f=NULL, Omega=NULL, Phi=NULL, LB=NULL,
        positive.data=FALSE, names=NULL, ...)
  {#  (... further arguments, currently disregarded)
   #  arg break.points=NULL, not yet supported
   #  obj should be (hat)B
   #     x =  B f + e  # NB no mean added to give x
   #   vector processes x, f, and e have times series matrices of data so
   #   calculation is eg t(B %*% t(f))

   if(is.null(f)) stop(" f must be specified.")

   if(!is.matrix(obj))
     stop("TSFmodel.default requires a loadings matrix (factor loadings) as the first argument.")

   if(ncol(obj)!=nseries(f))  stop("dimensions of obj (B) and  f do not agree.")
   if(is.null(names)) names <- paste("Series", seq(nrow(obj)))
   classed(list(B=obj, f=f, Omega=Omega, Phi=Phi, LB=LB, 
	   positive.data=positive.data,
           names=names,  #break.points=break.points, 
	   dots=list(...)), "TSFmodel") # constructor
  }

simulate.TSFestModel <- function(model, Cov=TSFmodel(model)$Omega, sd=NULL, noise=NULL, 
	rng=NULL, noise.model=NULL, ...) {
   simulate(TSFmodel(model), Cov=Cov, sd=sd, noise=noise, 
	rng=rng, noise.model=noise.model, ...)
   }

simulate.TSFmodel <- function(model, Cov=model$Omega, sd=NULL, noise=NULL,  
	rng=NULL, noise.model=NULL, ...)
   {#  (... further arguments, currently disregarded)
    # tframe and periods are taken from factors (f) 
    
    if ( is.null(Cov) & is.null(sd) & is.null(noise) & is.null(noise.model))
      stop("One of Cov, sd, noise, or noise.model, must be specified.")

    p <- if (is.matrix(model$B)) nrow(model$B) else nrow(model$B[[1]])
    sampleT <- periods(factors(model))
     noise <- makeTSnoise(periods(model$f), p, 1, noise=noise, rng=rng,
                        Cov=Cov, sd=sd, noise.model=noise.model)

    #use the calculation in explained but discard the class setting
    x <- unclass(explained(model)) + noise$w  
    if (model$positive.data && any( x < 0 )) {
        warning("negative simulated data values set to zero.")
        x[ x < 0 ] <- 0
        }
    attr(x, "noise") <- noise
    attr(x, "TSFmodel") <- model
    tframed(x, tf=factors(model), names=seriesNames(model))
}


factors <- function(x)UseMethod("factors")
factors.TSFmodel <- function(x) classed(x$f, c("TSFfactors", class(x$f)))
factors.TSFestModel <- function(x) {
   r <- factors(TSFmodel(x))
   trueM <- attr(x$data, "TSFmodel")
   if(!is.null(trueM)) attr(r, "true") <- factors(trueM)
   r
   }


factors.EstEval  <- function(x)
   {N <- length(x$result)
    r <- vector(mode="list", length=N)
    for (i in 1:N) r[[i]] <- factors(x$result[[i]])
    classed(list(result=r, truth=factors(x$truth)),
            c("factorsEstEval", "EstEval"))
  }

diff.TSFmodel  <- function (x, ...){
  x$f <- diff(x$f)
  x 
  }

diff.TSFestModel  <- function (x, ...){
  x$model$f <- diff(x$model$f)
  trueM <- attr(x$data, "TSFmodel")
  x$data <- diff(x$data)
  if(!is.null(trueM)){
     trueM$f <- diff(factors(trueM))
     attr(x$data, "TSFmodel") <- trueM
     }
  x 
  }

diff.TSFexplained  <- function (x, ...){
  tf <- tfdiff(tframe(x))
  r  <- tframed(diff(unclass(x)), tf=tf)
  d  <- attr(x, "data")
  if(!is.null(d)) attr(r, "data") <- tframed(diff(d), tf=tf)
  classed(r, c("TSFexplained", class(r))) 
  }

diff.TSFfactors  <- function (x, ...){
  tf <- tfdiff(tframe(x))
  r <- tframed(diff(unclass(x)), tf=tf)
  truef <- attr(x, "true")
  if(!is.null(truef)) attr(r, "true") <- tframed(diff(unclass(truef)), tf=tf)
  classed(r, c("TSFfactors", class(r))) 
  }

diff.factorsEstEval  <- function (x, ...){
  N <- length(x$result)
    r <- vector(mode="list", length=N)
    for (i in 1:N) r[[i]] <- diff(x$result[[i]])
    classed(list(result=r, truth=diff(x$truth)),
            c("factorsEstEval", "EstEval"))
  }

coef.TSFmodel    <- function(object, ...) {c(object$B, diag(object$Omega))}
coef.TSFestModel <- function(object, ...) {coef(TSFmodel(object))}

nfactors <- function(x) UseMethod("nfactors")
nfactors.TSFmodel <- function(x) {nseries(factors(x))}
nfactors.TSFestModel <- function(x) {nfactors(TSFmodel(x))}
nfactors.EstEval <- function(x) {nfactors(x$truth)}
nfactors.TSFfactors <- function(x) {nseries(x)}

seriesNames.TSFestModel <- function(x) {seriesNames(x$data)}
seriesNames.TSFmodel <- function(x) {x$names}

factorNames <- function(x) UseMethod("factorNames")
factorNames.TSFmodel    <- function(x) seriesNames(factors(x))
factorNames.TSFestModel <- function(x) factorNames(TSFmodel(x))
factorNames.EstEval     <- function(x) factorNames(x$truth)
factorNames.TSFfactors     <- function(x) seriesNames(x)



tfplot.TSFmodel <- function(x,..., tf=tfspan(x), start=tfstart(tf), end=tfend(tf), 
      series=seq(nfactors(x)),
      Title="Model factors", 
      lty = 1:5, lwd = 1, pch = NULL, col = 1:6, cex = NULL,
      xlab=NULL, ylab=factorNames(x), xlim = NULL, ylim = NULL, 
      graphs.per.page=5, par=NULL, mar=par()$mar, reset.screen=TRUE) {

   tfplot(factors(x),...,  tf=tf, start=start, end=end, 
	series=series, Title=Title, 
        lty=lty, lwd=lwd, pch=pch, col=col, cex=cex,
        xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
	graphs.per.page=graphs.per.page, 
	par=par, mar=mar, reset.screen=reset.screen)
   }

tfplot.TSFestModel <- function(x,...)  tfplot(factors(x), ...)

tfplot.TSFfactors <- function(x,..., tf=tfspan(x), start=tfstart(tf), end=tfend(tf), 
      series=seq(nfactors(x)),
      Title="Estimated factors (dashed) and true (solid)", 
      lty = c("dashed", "solid"), lwd = 1, pch = NULL, col = 1:6, cex = NULL,
      xlab=NULL, ylab=factorNames(x), xlim = NULL, ylim = NULL, 
      graphs.per.page=5, par=NULL, mar=par()$mar, reset.screen=TRUE) {

   tfplot(unclass(x), attr(x, "true"),  ..., 
        tf=tf, start=start, end=end, 
	series=series, Title=Title, 
        lty=lty, lwd=lwd, pch=pch, col=col, cex=cex,
        xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
	graphs.per.page=graphs.per.page, 
	par=par, mar=mar, reset.screen=reset.screen)
   }


tfplot.TSFexplained <- function(x,..., tf=tfspan(x), start=tfstart(tf), end=tfend(tf), 
      series=seq(nseries(x)),
      Title="Explained (dashed) and actual data (solid)", 
      lty = c("dashed", "solid"), lwd = 1, pch = NULL, col = 1:6, cex = NULL,
      xlab=NULL, 
      ylab=seriesNames(x), 
      xlim = NULL, ylim = NULL,
      graphs.per.page=5, par=NULL, mar=par()$mar, reset.screen=TRUE) {
    tfplot( unclass(x), attr(x, "data"), 
                Title=Title,
                tf=tf, start=start, end=end, series=series,  
                lty=lty, lwd=lwd, pch=pch, col=col, cex=cex,
                xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
		graphs.per.page=graphs.per.page, 
		par=par, mar=mar, reset.screen=reset.screen)
  }


LedermannBound  <- function(M) {
   if (is.matrix(M)) M <- ncol(M)
   if (1 != length(M)) stop("M must be an integer number of indicator variables.")
   
   # this works if only the biggest integer solution is needed.
   #max(seq(M)[(M - seq(M))^2 >= M + seq(M)])
   
   # solve (M^2-M) - (2M+1)k + k^2 = 0 for k
   r <- polyroot(c(M^2-M, -(2*M+1), 1))
   if (any(Im(r) > .Machine$double.eps))
            warning("LedermannBound has complex solution")
   r <- Re(r)
   r[(0 <=r) & (r <= M)]
   }

FAfitStats <- function(object, ...)UseMethod("FAfitStats")

FAfitStats.default <- function(object, diff.=TRUE, 
                           N=(nrow(object) - diff.), 
			   control=list(lower = 0.0001), ...) {
    if (!is.matrix(object)) stop("FAfitStats.default expects a matrix object.")
    corDX  <- cor(if (diff.) diff(object) else object)
    #covDX  <- cov(if (diff.) diff(object) else object)
    nvar   <- ncol(object)
    maxfact <- floor(LedermannBound(nvar))
    OmegaTot <- matrix(NA,nvar, maxfact+1)
    
    # Fit statistics could be calculate with either corDX or covDX. Factanal
    # transforms a covariance matrix to a correlation matrix: the loadings and 
    # uniquenesses it returns are for the standardized solution (correlation 
    #  matrix) in whether it is given a cov or cor matrix. This means the
    # loadings an uniquenesses must be converted bac to the unstandardized
    # scale, if the fit staistics are to be calculated with the cov.
    # Commented code below does the calculation using covariances.

    # zero factors
    fitStats   <- FAmodelFitStats(NULL, NULL, diag(corDX), corDX, N)
    OmegaTot[,1] <- diag(corDX)
    nm <- list(names(fitStats), c(0, seq(maxfact), "saturated"))
    
    for (j in 1:maxfact) { # loop over number of factors
     	# estimate parameters
    	FA <- factanal(factors = j, covmat = corDX, n.obs = N,
    			scores = "none", rotation = "none", control=control)
    	#FA <- factanal(factors = j, covmat = covDX, n.obs = N,
    	#		scores = "none", rotation = "none", control=control)
     	OmegaTot[,j+1] <- FA$uniquenesses
        
        # compute fit statistics

 	fitStats <- cbind(fitStats, 
            FAmodelFitStats(FA$loadings, diag(1,j,j), FA$uniquenesses, corDX, N))
	#D <- diag(sqrt(diag(covDX)))
 	#fitStats <- cbind(fitStats, 
	#   FAmodelFitStats(D %*% FA$loadings, diag(1,j,j),
	#                               diag(covDX) * FA$uniquenesses, covDX, N))
        }

    Hey <- apply(control$lower == OmegaTot, 2, any)
    if(any(Hey)) warning("Heywood cases: ", 
                   paste((0:maxfact)[Hey], collapse=" "), " factor model(s)")

    # saturated model
    k <- maxfact
    M <- nvar
    nparc <- M * k + M - (k *(k-1))/2 # no. of param's corrected for rotation
    fitStats <- cbind(fitStats, c(
        0,0,NA,0,         # chisq=0, df=0, pval=na, delta=0
        NA,1,1,  NA,      # RMSEA=na, RNI=1, CFI=1, MCI=na
        1,1,0,            # GFI=1, AGFI=1, AIC=0
        (1 + log(N)) * nparc,  # CAIC
        log(N) * nparc,        # SIC
        (2 * nparc)/N,         # CAK
        2 * nparc/(N-M-1) ))    # CK


    dimnames(fitStats) <- nm
    
    # differences between consecutive models

    seqfitStats <- NULL
    for (j in 1:(maxfact+1))  seqfitStats <- cbind(seqfitStats, 
        c( fitStats["chisq",j] - fitStats["chisq",j+1],
	   fitStats["df",   j] - fitStats["df",   j+1],
	   pchisq(fitStats["chisq",j] - fitStats["chisq",j+1], 
	          fitStats["df",   j] - fitStats["df",   j+1], lower.tail=FALSE)))
    
    nm <- dimnames(fitStats)[[2]]
    dimnames(seqfitStats) <- list(c("chisq", "df", "pval"), 
                                  paste(nm[-1], "vs", nm[-length(nm)]))
        
    list(fitStats=fitStats, seqfitStats=seqfitStats)  #, OmegaTot= OmegaTot)
    }

FAfitStats.TSFestModel <- function(object, diff.=TRUE,
                             N=(nrow(object$data) - diff.), ...) {
    # This uses unstandardized B and Omega (and covDX rather than corDX).
    # This should be the same as standardized for MLE, but not for other
    #  estimation methods. Standardized may be better conditioned numerically,
    #  so might be perferred for ML, but probably does not seem to make much
    #  differene in simple tests. Might consider using both.
    
    X <- if(diff.) diff(object$data) else object$data
    FAmodelFitStats(loadings(object), object$model$Phi, diag(object$model$Omega),
                   cov(X), N)
    } 


FAmodelFitStats <- function(B, Phi, omega, S, N) {
  tr <- function(A) {sum(diag(A))} # local function, trace of a matrix

  # Fit statistics of FA model, based on standard likelihood.
  # No consistency checks are made.
  #
  # See, e.g., Wansbeek, T., & Meijer, E. (2000). Measurement error and latent
  #   variables in econometrics. Amsterdam: North-Holland. (W&M below)
  #
  # B     = loadings
  # Phi   = cov. matrix of factors
  # omega = vector of error variances
  # S     = sample covariance matrix, or correlation matrix if B and omega are
  #         standardized. (see the notes in FAfitStats.default)
  # N     = sample size (Many authors prefer sample size - 1.)
  # k     = number of factors (may be zero)

  # numbers of variables and factors

  M <- nrow(S)
  k <- if (is.null(B)) 0 else ncol(B)

  # Saturated model: all elements of cov. matrix are free parameters.

  const <- log(det(S)) + M

  # Null model: independence model (W&M, p. 305).

  Sigma0 <- diag(c(diag(S)))
  chisq0 <- N * (log(det(Sigma0)) + tr(solve(Sigma0) %*% S) - const)
  df0    <- 0.5 * (M^2 - M)
  delta0 <- max(0, chisq0 - df0)

  # Target model, i.e., the one from which B, Phi, and omega are obtained.

  # Model-implied covariance matrix
  if (k == 0)  Sigma <- diag(c(omega))
  else if (is.null(Phi)) Sigma <- B %*% t(B) + diag(c(omega))
  else                   Sigma <- B %*% Phi %*% t(B) + diag(c(omega))


  # Chi-square statistic, its degrees of freedom, and its p-value (W&M, p. 298).
  # Note: the df takes the rotational freedom into account (cf. W&M, p. 169).
  chisq <- N * (log(det(Sigma)) + tr(solve(Sigma) %*% S) - const)
  df    <- 0.5 * ( (M - k)^2  -  (M + k) )
  pval  <- pchisq(chisq, df, lower.tail=FALSE)

  # Estimate of noncentrality parameter (W&M, p. 307).
  delta <- max(0, chisq - df)

  # Comparative fit index (W&M, p. 307).
  if(chisq0 <= df0)       CFI <- 0 # Null model fits very well: extremely unlikely.
  else if(chisq <= df)    CFI <- 1 # Target model fits very well.
  else if(delta0 < delta) CFI <- 0 # Null model fits better than target model: also extremely unlikely.
  else        CFI <- 1 - delta/delta0  # The most common situation: null model fits very badly,
                             # target model fits better, but not perfectly.

  # Root mean square error of approximation (W&M, p. 309).
  RMSEA <- if (df > 0)  sqrt(delta/(N * df)) else  Inf


  # Dozens of other fit indexes possible. See output of LISREK, EQS,
  # Amos, Mplus, and Mx, and the paper by Ogasawara (SEM, ca. 2001).
  # Most are very bad. Here are some possibilities, from W&M (chap. 10),
  # Hu & Bentler (1995), and Bollen (1989, pp. 256-281):

   # NFI  <- 1 - chisq/chisq0 # Normed fit index
   # TLI  <- (chisq0/df0 - chisq/df)/(chisq0/df0 - 1) # Tucker-Lewis Index,
  	       # also called Nonnormed fit index
   # BL86 <- 1 - (chisq/df)/(chisq0/df0) # Bollen (1986)
   # BL89 <- (chisq0 - chisq)/(chisq0 - df) # Bollen (1989)
   RNI  <- 1 - delta/delta0 # Relative noncentrality index
   MCI  <- exp(-0.5 * (chisq - df)/N) # McDonald's centrality index
  
  # # Some indexes by Joreskog & Sorbom (1981, 1986):
  
   T1	<- solve(Sigma) %*% S  # Temporary matrix
   T2	<- T1 - diag(1,M,M)    # Temporary matrix
   GFI  <- 1 - tr(T2 %*% T2)/tr(T1 %*% T1) # Goodness of fit index
   AGFI <- 1 - ((M * (M+1))/(2 * df)) * (1 - GFI) # Adjusted GFI
   #RMR <- sqrt(((sum( (   c(S - Sigma))^2 ) +
   #		 sum( (diag(S - Sigma))^2 ))/(M *(M+1))) #Root mean-square residual

  # # Hoelter's (1983) Critical N; note that the significance level
  # # alpha must be given. E.g.:
  # # alpha <- 0.05
  # # Apparently, this is the definition of Hoelter:
  # CN <- (qnorm(alpha/2, lower.tail=FALSE) + sqrt(2*df-1))^2/(2*chisq/N) + 1
  # # but I've also seen the following, which makes more sense:
  # # CN <- qchisq(alpha, df, lower.tail=FALSE)/(chisq/N) + 1
  # # although both are not very useful.

  # # Some information criteria; accounting for rotational freedom
   nparc <- M * k + M - (k *(k-1))/2 # no. of param's corrected for rotation
   AIC   <- chisq - 2 * df # or chisq + 2 * nparc # Akaike's info. crit.
   CAIC  <- chisq + (1 + log(N)) * nparc # Consistent AIC
   SIC   <- chisq + log(N) * nparc # Schwarz's Bayesian info. crit.
   CAK   <- (chisq + 2 * nparc)/N  # Cudeck & Browne's rescaled AIC
   CK	 <- chisq/N + 2 * nparc/(N-M-1) # Cudeck & Browne's cross-val. index

  r <- c(chisq, df, pval, delta, RMSEA, RNI, CFI,
         MCI, GFI, AGFI, AIC, CAIC, SIC, CAK, CK) #NFI, TLI, BL86, BL89

  names(r) <- c("chisq", "df", "pval", "delta", "RMSEA", "RNI", "CFI",
         "MCI", "GFI", "AGFI", "AIC", "CAIC", "SIC", "CAK", "CK")
  r
}


summary.TSFmodel <- function(object, ...)
 {classed(list(k=nfactors(object), M=nrow(object$B),
      Snames=seriesNames(object), Fnames=factorNames(object),
      N=periods(factors(object)), S=start(factors(object)), E=end(factors(object)),
      Omega=!is.null(object$Omega), 
      Phi=!is.null(object$Phi), 
      LB=!is.null(object$LB), 
      positive.data=object$positive.data
      ), "summary.TSFmodel")

 }

print.summary.TSFmodel <- function (x, ...)
  {cat(x$k, " factors: ",     x$Fnames, "\n")
   cat(x$M, " indicators: ",  x$Snames, "\n")
   cat("factors have ", x$N, " periods from:", x$S, " to ", x$E, "\n")
   cat("Omega ", (if(x$Omega) "is" else "is not"), " specified.\n")
   cat("Phi   ", (if(x$Phi) "is" else "is not"), " specified.\n")
   cat("LB    ", (if(x$LB) "is" else "is not"), " specified.\n")
   cat("positive.data ", (if(x$positive.data) "is" else "is not"), " specified.\n")
  }

summary.TSFestModel <- function(object, ...)
 {
  fitStats <- FAfitStats(object)
  est  <- TSFmodel(object)
   
  barx     <- colMeans(object$data)
  barx.est <- colMeans(explained(est))
  hatk     <- colMeans(factors(est))

  barDx     <-  colMeans(diff(object$data ))
  barDx.est <-  colMeans(diff(explained(est)))
  hatDk     <-  colMeans(diff(factors(est)))

  true <- attr(object$data, "TSFmodel")
  if (is.null(true)) 
      B.true <- hatk.true <- hatDk.true <- NULL 
  else  {
      B.true     <- true$B
      hatk.true  <- colMeans(true$f)
      hatDk.true <- colMeans(diff(factors(true)))
      }
  
  classed(list(Snames=seriesNames(object),Fnames=factorNames(est),
      fitStats=fitStats, B.estimate=est$B,   B.true=B.true,
      #stdB.estimate=standardizedLoadings(object), 
      DstdB.estimate=DstandardizedLoadings(object),
      barDx=barDx, barDx.est=barDx.est,   barx=barx,  barx.est=barx.est,
      hatDk=hatDk, hatDk.true=hatDk.true, hatk=hatk,  hatk.true=hatk.true),
	  "summary.TSFestModel")
  }

print.summary.TSFestModel <- function (x, ...)
  {cat("     Estimated loadings:\n"); print(x$B.estimate)
   cat("\n     Standardized (using differenced data covariance):\n")
   print(x$DstdB.estimate)
   #cat("\n     Standardized (using undifferenced data covariance):\n")
   #print(x$stdB.estimate) 
   if (!is.null(x$B.true))
     {cat("\n   true loadings:\n"); print(x$B.true)
      cat("\n   loadings estimation error:\n"); print(x$B.estimate - x$B.true)
     }


   z <- rbind(x$barx.est,x$barx,  x$barx.est - x$barx)
   dimnames(z) <- list(c("explained","actual","error"), x$Snames)
   cat("\n                 Mean of data:\n"); print(z)

   z <- rbind(x$barDx.est,x$barDx,  x$barDx.est - x$barDx)
   dimnames(z) <- list(c("explained","actual","error"), x$Snames)
   cat("\n		  Mean of differenced data:\n"); print(z) 

   if (!is.null(x$hatk.true))
     {z <- rbind(x$hatk, x$hatk.true,  x$hatk - x$hatk.true)
      dimnames(z) <- list(c("estimated","true","error"), x$Fnames)
     }
   else
     {z <- x$hatk
      names(z) <- x$Fnames
     }
   cat("\n     Mean of factors:\n"); print(z)
   
   if (!is.null(x$hatDk.true))
     {z <- rbind(x$hatDk, x$hatDk.true,  x$hatDk - x$hatDk.true)
      dimnames(z) <- list(c("estimated","true","error"), x$Fnames)
     }
   else 
     {z <- x$hatDk
      names(z) <- x$Fnames
     }
   cat("\n    Mean of differenced factors:\n"); print(z)

   cat("\n   Fit statistics:\n")
   print(x$fitStats)
   invisible(x)
  }



distribution.factorsEstEval <- function (obj, ..., bandwidth = "nrd0",
        cumulate=TRUE, graphs.per.page = 5, Title=NULL)
  {# if cumulate is true then a distribution is plotted, otherwise,
   # a time series graph of the true and one 1 sd bands
    truth <- obj$truth
    r <- array(NA, c(length(obj$result), dim(truth)))
    otherobj <- list(...)
    obr <- list()
    for (ob in otherobj)
      {if (! testEqual(truth, ob$truth))
                    warning("object true values do not correspond.")
       rx <- r
       for (i in 1:length(ob$result)) rx[i,,] <- ob$result[[i]] - truth
       obr <- append(obr, list(rx))
      }
    for (i in 1:length(obj$result)) r[i,,] <- obj$result[[i]] - truth
    xlab <- "factor "
    old.par <- par(mfcol = c(min(graphs.per.page, ncol(truth)), 1),
        mar = c(5.1, 4.1, 4.1, 2.1), no.readonly = TRUE)
    on.exit(par(old.par))
    if (cumulate)
       for (i in 1:ncol(truth))
         {if (is.R())
	    {rd <- density(c(r[,,i]), bw = bandwidth)
             rdy <- rd$y
	     for (rx in obr)
	         rdy <- cbind(rdy, density(c(rx[,,i]), bw = bandwidth)$y)
            }
	  else 
	    {rd <- ksmooth(r[,, i], bandwidth = bandwidth)
	     rdy <- rd$y
	    }
	  matplot(rd$x, rdy, type = "l",
	     ylab = "density", xlab = paste(xlab, i), main = "")
          if (i == 1) title(main = Title)
	 }
    else
      {rd <- apply(r,c(2,3), FUN="var")^0.5
       tfplot(truth, truth + rd, truth - rd,
                   Title=Title, graphs.per.page = graphs.per.page)
      }
    invisible()
}

checkResiduals.TSFestModel <- function (obj, diff.=TRUE, ...) {
	res <- if (diff.) diff(explained(obj)) - diff(obj$data)
	           else explained(obj) - obj$data
	seriesNames(res) <- seriesNames(obj$data)
	cat("residual covariance matrix\n")
	cv <- cov(res)
	print(cv)
	cat("\nsum of trace cov: ", sum(diag(cv)), "\n")
	cat("sum of abs (off-diag of cov): ", sum(abs(cv - diag(cv))), "\n")
	checkResiduals(res, ...)
}


permusign <- function(B, Btarget, Phi=diag(1,ncol(B))) {

# Selects the permutation and signs of the columns of the factor loadings B
# that resembles the Btarget matrix the most.
# Phi matrix (cov of factors) may need to be reordered for the permutation.
  
  ############################## local functions
  permute <- function(x){
     # with thanks to Bill Venables
     if (length(x) <= 1) as.matrix(x) else{
  	 M <- NULL
  	 for (i in seq(length(x))) M <- rbind(M, cbind(x[i], Recall(x[-i])))
  	 M
  	 }
     }
  signsw <- function(Bprop, Bnew, Btarget){
     # compare distance of Bprop from Btarget and also Bprop with column 
     # signs switched. If the best of these is better than Bnew to Btarget
     #  return the column signs (1 or -1), otherwise return NULL
     signs <- rep(1, ncol(Bprop))
     d1 <- colSums((Btarget - Bprop)^2)
     d2 <- colSums((Btarget + Bprop)^2) # all col signs reversed
     if ( sum(pmin(d1, d2)) < sum((Btarget - Bprop %*% diag(signs))^2))
     	signs <- 2 * (d1 < d2) - 1
     # the fuzz (1e-12) seems to be necessary to avoid rounding error causing
     # T when things should be equal,with the result that random 
     #  permutations occur.
     if (sum((Btarget - Bnew)^2) > 1e-12 +
         sum((Btarget - Bprop %*% diag(signs))^2) ) signs else NULL
     }

  ############################## end local functions

  P    <- permute(seq(ncol(B))) # permutation matrix
  Bnew   <- B
  PhiNew <- Phi
  for (j in seq(nrow(P))) {
    Bprop <- B[,P[j,]]
    signs <- signsw(Bprop, Bnew, Btarget)
    if(!is.null(signs)){
       #cat(j, ":", P[j,], signs)
       Bnew   <-  Bprop %*% diag(signs)
       PhiNew <-  (Phi[P[j,],P[j,]]) * outer(signs,signs)
       }
    }
  list(loadings=Bnew,Phi=PhiNew)
}



summaryStats <- function(object, ...) UseMethod("summaryStats")

summaryStats.TSFmodelEstEval <- function(object, ...) {

  N <- length(object$result)
  if (N <2 ) stop("This requires more than one replication.")

  meanhatf <- sdhatf <- meanhatDf <- sdhatDf <- meanhatPCf <- 
              sdhatPCf <- meanhatB <- sdhatB <- 
	      factanalConverged <- rotationConverged <- 0

  for (m in object$result) { 
      meanhatf  <- meanhatf   + factors(m)
      sdhatf	<- sdhatf     + factors(m)^2
      meanhatDf <- meanhatDf  + diff(factors(m))
      sdhatDf   <- sdhatDf    + diff(factors(m))^2
      meanhatPCf<- meanhatPCf + percentChange(factors(m))
      sdhatPCf  <- sdhatPCf   + percentChange(factors(m))^2
      meanhatB  <- meanhatB   + TSFmodel(m)$B
      sdhatB	<- sdhatB     + TSFmodel(m)$B^2
      if (!is.null(TSFmodel(m)$dots$factanalConverged) && !TSFmodel(m)$dots$factanalConverged)
          factanalConverged <- factanalConverged + 1
      if (!is.null(TSFmodel(m)$dots$rotationConverged) && !TSFmodel(m)$dots$rotationConverged)
          rotationConverged <- rotationConverged + 1
      }
 
  true <- factors(object$truth)
  tf <- tframe(true)
  
  meanhatf   <- tframed(meanhatf /N, tf)
  meanhatDf  <- meanhatDf /N
  meanhatPCf <- meanhatPCf /N
  meanhatB   <- meanhatB /N
 
  sdhatf   <- tframed(sqrt(sdhatf /N - meanhatf^2), tf)
  sdhatDf  <- sqrt(sdhatDf  /N - meanhatDf^2)
  sdhatPCf <- sqrt(sdhatPCf /N - meanhatPCf^2)
  sdhatB   <- sqrt(sdhatB   /N - meanhatB^2)
  
  list(true=true, Btrue=TSFmodel(object$truth)$B,
  	meanhatf   =  meanhatf, 
  	meanhatDf  =  meanhatDf, 
  	meanhatPCf =  meanhatPCf, 
  	meanhatB   =  meanhatB,   
  	sdhatf     =  sdhatf,	
  	sdhatDf    =  sdhatDf,  
  	sdhatPCf   =  sdhatPCf, 
  	sdhatB     =  sdhatB,
        factanalConverged = factanalConverged,
        rotationConverged = rotationConverged)
  }


summary.TSFmodelEstEval <- function(object, ...) {
  sm <- summaryStats(object, ...)
  classed(list(
      meanhatf.error  = colMeans(sm$meanhatf  - sm$true),
      meanSDhatf      = colMeans(sm$sdhatf),
      meanhatDf.error = colMeans(sm$meanhatDf - diff(sm$true)),
      meanSDhatDf     = colMeans(sm$sdhatDf),
      meanhatPCf.error= colMeans(sm$meanhatPCf -
                               percentChange(sm$true)),
      meanSDhatPCf    = colMeans(sm$sdhatPCf),
      meanhatB.error  = sm$meanhatB - sm$Btrue, 
      SDhatB          = sm$sdhatB,
      factanalConverged = sm$factanalConverged,
      rotationConverged = sm$rotationConverged), "summary.TSFmodelEstEval")
  }

print.summary.TSFmodelEstEval <- function(x, digits = options()$digits, ...) {
  cat("    mean hat{f} error\n") ; print(x$meanhatf.error, digits=digits)
  cat("\n    mean  SD hat{f}\n")   ; print(x$meanSDhatf, digits=digits)
  cat("\n    mean diff hat{f} error\n");  print(x$meanhatDf.error,  digits=digits)
  cat("\n    mean %change hat{f} error\n");  print(x$meanhatPCf.error,  digits=digits)
  cat("\n    mean hat{B} error") ; print(x$meanhatB.error, digits=digits) 
  cat("\n      SD hat{B}")       ; print(x$SDhatB, digits=digits) 
  cat("\n    Estimates NOT converged: ", x$factanalConverged) 
  cat("\n    Rotations NOT converged: ", x$rotationConverged) 
  cat("\n") 
  invisible(x)
  }


tfplot.TSFmodelEstEval <- function(x, diff.=FALSE,  percentChange.=FALSE,
        PCcentered.=FALSE, summary.=TRUE, ...) {

   #if summary. is FALSE then all of the factors are plotted,
   # otherwise the mean and 1 SD bounds are plotted as follows:
   #if diff. is TRUE then the differenced factors are plotted
   #if percentChange. is TRUE then the PC factors are plotted
   #if PCcentered. is TRUE then the PC factors less means are plotted
   # otherwise the undifferenced factors are plotted

   true <- factors(x$truth)
   
   if(!summary.) {
     tfplot(factors(x), tf = tframe(true),
        truth = true,
        Title = "Estimated (and true) results",
      ylab = seriesNames(truth), remove.mean = FALSE, graphs.per.page = 5,
       mar = par()$mar, reset.screen = TRUE, ...)
      } 
   else {
      sm <- summaryStats(x)

      if(diff.) { # factor difference
    	  tfplot(diff(true),
    		 sm$meanhatDf, 
    		 sm$meanhatDf  + 1.96 * sm$sdhatDf,  
    		 sm$meanhatDf  - 1.96 * sm$sdhatDf)
    	  }
      else if(percentChange.){ # factor growth rates
    	  tfplot(percentChange(true),
    		 sm$meanhatPCf, 
    		 sm$meanhatPCf  + 1.96 * sm$sdhatPCf,  
    		 sm$meanhatPCf  - 1.96 * sm$sdhatPCf)
    	  }
      else if(PCcentered.){ # factor growth rates: bias (clearer picture?)
    	  growth <- percentChange(true)
	  tfplot(sm$meanhatPCf - growth, 
    		sm$meanhatPCf  + 1.96 * sm$sdhatPCf - growth, 
    		sm$meanhatPCf  - 1.96 * sm$sdhatPCf - growth)
          }
      else { # factors
    	  tfplot(true,
    		 sm$meanhatf, 
    		 sm$meanhatf  + 1.96 * sm$sdhatf, 
    		 sm$meanhatf  - 1.96 * sm$sdhatf)
    	  }
      }
   invisible(x)
   }


predict.TSFmodel <- function(object, newdata = NULL, factorNames.=factorNames(object), ...){
      # prediction of factors with new data
      if (is.null(newdata)) stop("newdata must be supplied.")
      LB <- object$LB
      tframed(newdata %*% t(LB), tframe(newdata), names=factorNames.) #hatf
      }

predict.TSFestModel <- function(object, newdata=NULL, factorNames.=factorNames(object), ...){
      if (is.null(newdata)) newdata <- object$data
      LB <- object$model$LB
      tframed(newdata %*% t(LB), tframe(newdata), names=factorNames.) #hatf
      }

explained <- function(object, ...)UseMethod("explained")
explained.TSFestModel <- function (object, ...)
 {r <- explained(TSFmodel(object), names=seriesNames(object$data), ...) 
  attr(r, "data") <- object$data
  r
 }

explained.TSFmodel <- function (object, names=object$names, ...) {
  # portion of data explained by factors
  f <- object$f
  if (is.matrix(object$B)) x <- t(object$B %*% t(f))
#  else {warning("this part of explained.TSFmodel is not tested.")
#	if (is.list(B)  &&  !is.list(break.points))
#	     stop("a list of break.points must be specified for a list of B.")
#	if (length(break.points)+1 != length(B))
#	     stop("the list of break.points must be one shorter than the list of B.")
#	x <- NULL
#	starts <- list(tfstart(f))
#	for (i in break.points) starts <- append(starts, list(i))
#	ends <- list()
#	for (i in break.points)
#	   ends <- append(ends, list(addDate(i, -1, tffrequency(f))))
#	ends <- append(ends, list(tfend(f)))
#	for (i in seq(length(B))) x <- rbind(x,
#	 t(B[[i]] %*%  t(tfwindow(f, start=starts[[i]], end=ends[[i]]))) )
#	}
  x <- tframed(x, tf=tframe(f), names=names) 
  classed(x, c("TSFexplained", class(x))) 
  }



#######################################

estTSF.ML <- function(y, p, diff.=TRUE, 
                      rotation=if(p==1) "none" else "quartimin", 
		      methodArgs=NULL, 
		      normalize=TRUE, eps=1e-5, maxit=1000, Tmat=diag(p),
		      BpermuteTarget=NULL,
                      factorNames=paste("Factor", seq(p))) {

      # Estimate parameters using standard (quasi) ML factor analysis
      # (on the correlation matrix and then scaled back).
      # factanal always uses the cor matrix, so standardizing does not affect 
      # the solution. Both standardized and not can be calculated after.
      # With non ML methods this solutions may differ (and working with cov  
      # rather than cor is probabably better.

      if (p < 1) stop("number of factors must be a positive integer.")
     
      zz <- if (diff.) diff(y) else y
      zz <- sweep(zz,2,colMeans(zz), "-")
      #Sigma <- cov(zz) which is the same as
      Sigma  <- crossprod(zz)/(periods(zz) - 1)
 
      z <- factanal(covmat = Sigma, factors=p, scores="none",
    			  rotation="none", n.obs=(periods(y) - diff.))
      #  above should be the same as
      #  z <- factanal(if(diff.) diff(y) else y, factors=p,
      #	                 scores="none", rotation="none")

      factanalConverged <- z$converged
      
      stds	<- sd(if (diff.) diff(y) else y)
      uniquenesses <- z$uniquenesses
      hatOmega  <- stds * uniquenesses * stds
       
      # for debugging compare:  hatOmega - Omega, hatOmega, Omega

      # z$loadings is orth solution
      if (rotation == "none") {
         hatB <- diag(stds) %*% z$loadings              
	 Phi  <- NULL
	 rotationConverged <- NULL
 	 }
      else {	 
    	 rotB <- GPFoblq(z$loadings, Tmat=Tmat, 
			 normalize=normalize, eps=1e-5, maxit=1000,
    	                 method=rotation, methodArgs=methodArgs)
    	 hatB <- diag(stds) %*% rotB$Lh
	 rotationConverged <- rotB$convergence
     
    	 # Make sure columns are ordered properly and have the correct signs.
    	 Phi <- rotB$Phi
    	 if (! is.null(BpermuteTarget)) {
    	    z  <- permusign(hatB, BpermuteTarget, Phi=Phi)
    	    hatB <- z$loadings
    	    Phi  <- z$Phi
    	    }
	 }
      dimnames(hatB) <- list(seriesNames(y), factorNames)

      ### Compute Bartlett-predicted factor scores 
      #BO <- t(hatB) %*% diag(1/hatOmega) or
      #BO <- crossprod(hatB, diag(1/hatOmega) )
      #LB   <- solve(BO %*% hatB, BO)
      # above is the same as
      #LB   <- solve(t(hatB) %*% solve(hatOmega) %*% hatB) %*%
      #              t(hatB) %*% solve(hatOmega) 
      
      Sigma <- if(is.null(Phi)) hatB %*% t(hatB) + diag(hatOmega) else
                        hatB %*% Phi %*% t(hatB) + diag(hatOmega)
      BO <- crossprod(hatB, solve(Sigma) )
      LB   <- solve(BO %*% hatB, BO)

      dimnames(LB) <- list(factorNames, seriesNames(y))

      model <- TSFmodel(hatB, 
                        f=tframed(y %*% t(LB), tframe(y), names=factorNames), #hatf
			Omega=diag(hatOmega), Phi=Phi,
			LB=LB,   
	                positive.data=all(0<y), 
		        factanalConverged=factanalConverged,
		        rotationConverged=rotationConverged )
      # factanalConverged and rotationConverged should really be in est rather
      #  than model, but then that has to be passed back from EstEval
      classed(list(model=model, data=y, 
             estimates=list(estimation="estTSF.ML", diff.=diff., rotation=rotation,
	               uniquenesses = uniquenesses, 
	               BpermuteTarget=BpermuteTarget)), "TSFestModel")
      }

#######################################
#######################################

# this file has (tentative, experimental) code for information tests

################  tests#########
if (FALSE) {
require("tsfa")
require("numDeriv")
require("sandwich") # To computes HAC estimators.

### Construct data
data("CanadianMoneyData.asof.6Feb2004", package="CDNmoney")

cpi <- 100 * M1total / M1real
seriesNames(cpi) <- "CPI"
popm <- M1total / M1PerCapita
seriesNames(popm) <- "Population of Canada"

z <- tframed(tbind(
    MB2001,
    MB486 + MB452 + MB453 ,
    NonbankCheq,
    MB472 + MB473 + MB487p,
    MB475,
    NonbankNonCheq + MB454 + NonbankTerm + MB2046 + MB2047 + MB2048 +
    MB2057 + MB2058 + MB482),
    names=c("currency", "personal cheq.", "NonbankCheq",
    "N-P demand & notice", "N-P term", "Investment")
    )

z <- tfwindow(z, start=c(1986,1))
if( all(c(2003,12) ==end(z))) z <-tfwindow(z, end=c(2003,11))
MBcomponents <- 1e8 * z/matrix(tfwindow(popm * cpi,tf=tframe(z)),periods(z),6)

   z <- tsfa.with.seasons(MBcomponents, p=2, nseason=1,  diff.=TRUE,numDeriv.method="simple")
   zz <- tsfa.with.seasons(MBcomponents, p=2, nseason=1,  diff.=TRUE,numDeriv.method="Richardson")
   z5 <- tsfa.with.seasons(MBcomponents, p=2, nseason=1,  diff.=TRUE,
        numDeriv.method="Richardson",
	numDeriv.method.args=list(eps=1e-4, d=0.0001, r=8, v=2))
   z6 <- tsfa.with.seasons(MBcomponents, p=2, nseason=1,  diff.=TRUE,
        numDeriv.method="Richardson",
	numDeriv.method.args=list(eps=1e-4, d=0.0001, r=10, v=2))
  max(abs(z5$Acov.theta - z6$Acov.theta))
#4092040
  max(abs(zz$Acov.theta - z5$Acov.theta))

   z5e <- tsfa.with.seasons(MBcomponents, p=2, nseason=1,  diff.=TRUE,
        rotation="eiv",
	numDeriv.method="Richardson",
	numDeriv.method.args=list(eps=1e-4, d=0.0001, r=8, v=2))
   z6e <- tsfa.with.seasons(MBcomponents, p=2, nseason=1,  diff.=TRUE,
        rotation="eiv",
        numDeriv.method="Richardson",
	numDeriv.method.args=list(eps=1e-4, d=0.0001, r=10, v=2))
  max(abs(z5$Acov.theta - z5e$Acov.theta))
# fairly big as expected
  max(abs(z5e$Acov.theta - z6e$Acov.theta))
#4092040  why is this the same as above ???

   zzz <- tsfa.with.seasons(MBcomponents, p=2, nseason=12, diff.=TRUE,numDeriv.method="simple")
   zzzz <- tsfa.with.seasons(MBcomponents, p=2, nseason=12, diff.=TRUE,numDeriv.method="Richardson")

 # there are some big differences here, but they don't amount to much in the criteria
  max(abs(100*(z$Acov.theta - zz$Acov.theta)/(abs(zz$Acov.theta)+1e-6)))
  max(abs(z$Acov.theta - zz$Acov.theta))
  max(abs(diag(z$Acov.theta - zz$Acov.theta)))

 # NFI 
   (z$chisq.null - z$chisq)/z$chisq.null
   (zz$chisq.null - zz$chisq)/zz$chisq.null
   (zzz$chisq.null - zzz$chisq)/zzz$chisq.null
   (zzzz$chisq.null - zzzz$chisq)/zzzz$chisq.null

 # RNI 
  M <- NCOL(MBcomponents) 
  k <- 2
  T <- periods(MBcomponents) - 1  # -1 because diff=TRUE
  nu <- ((M-k)^2 - (M+k))/2
  nu.null <- (M * (M+1))/2 - M  # = ((M-k)2 - (M+k))/2 with k = 0
   ((z$chisq.null - nu.null) - (z$chisq- nu))/(z$chisq.null - nu.null)
  ((zz$chisq.null - nu.null) - (zz$chisq- nu))/(zz$chisq.null - nu.null)
  ((zzz$chisq.null - nu.null) - (zzz$chisq- nu))/(zzz$chisq.null - nu.null)
  ((zzzz$chisq.null - nu.null) - (zzzz$chisq- nu))/(zzzz$chisq.null - nu.null)

 # RMSEA 
   sqrt(max((z$chisq - nu)/(nu * T), 0))
   sqrt(max((zz$chisq - nu)/(nu * T), 0))
   sqrt(max((zzz$chisq - nu)/(nu * T), 0))
   sqrt(max((zzzz$chisq - nu)/(nu * T), 0))

 # AIC 
      z$chisq + 2*(M*(M+1)/2 - nu)
     zz$chisq + 2*(M*(M+1)/2 - nu)
    zzz$chisq + 2*(M*(M+1)/2 - nu)
   zzzz$chisq + 2*(M*(M+1)/2 - nu)

 # CAIC 
      z$chisq + (M*(M+1)/2 - nu) * (1+log(T))
     zz$chisq + (M*(M+1)/2 - nu) * (1+log(T))
    zzz$chisq + (M*(M+1)/2 - nu) * (1+log(T))
   zzzz$chisq + (M*(M+1)/2 - nu) * (1+log(T))

 # BIC 
      z$chisq + (M*(M+1)/2 - nu) * log(T)
     zz$chisq + (M*(M+1)/2 - nu) * log(T)
    zzz$chisq + (M*(M+1)/2 - nu) * log(T)
   zzzz$chisq + (M*(M+1)/2 - nu) * log(T)
}  # end if FALSE

################  from cov02#########

HAC.of.s.with.seasons <- function(y, nseason=1, diff.=TRUE, meat="meatHAC") {

    # Compute a heteroskedasticity and autocorrelation consistent (HAC)
    # estimator of the asymptotic covariance matrrix of the sample covariance
    # matrix, for a seasonal time series.
    # If diff. == TRUE then the first differences of y are used instead of y.
    #
    # nseason = the number of seasons. If <= 1, then no seasons, and the
    # means and covariances are assumed to be the same at each time point.
    # If > 1, then means and covariances are to be computed for each
    # season separately. The return matrix zz is a matrix consisting of
    # dummy variables, with as many columns as there are seasons.

    # First, assemble the data.

    dd <- data.for.HAC.with.seasons(y, nseason, diff.)

    # Compute the residuals from the regression of yy on zz,
    # i.e., the deviations from the seasonal means.

    # lm crashes on large datasets, so do it the old-fashioned way

    # mvr <- lm(dd$yy ~ dd$zz - 1) # -1 removes the superfluous intercept
    # mvm <- mvmean(resid(mvr))

    Gam <- solve(crossprod(dd$zz), crossprod(dd$zz,dd$yy)) # (Z'Z)^{-1}Z'Y
    EE  <- dd$yy - dd$zz %*% Gam
    mvm <- mvmean(EE)

    # Compute a HAC estimator for the covariance matrix of the means
    # of the columns of dd$yy. At some point, allow the user to pass
    # options to this function. (Paul?)

    if (meat == "meat") {
        A <- meat(mvm)    # outer product estimate
    } else if (meat == "meatHC") {
        A <- meatHC(mvm)  # HC estimate
    } else if (meat == "meatHAC") {
        A <- meatHAC(mvm) # HAC estimate
    } else {
        stop("Unknown meat function")
    }

    # Transform to a HAC estimator for the covariance matrix of s, using
    # the delta method, and return the result.

    return(HAC.s.from.HAC.m(A, colMeans(y)))
}


### The Jacobian of the theta.tsfa function.

jacobian.tsfa <- function(s, p, est="factanal", n.obs=NA,
        estArgs=list(scores="none", control=list(opt=list(maxit=10000))),
        rotation=if(p==1) "none" else "quartimin", rotationArgs=NULL,
        GPFargs=list(Tmat=diag(p),normalize=TRUE, eps=1e-5, maxit=1000),
        BpermuteTarget=NULL, factorNames=paste("Factor", seq(p)),
        numDeriv.method="Richardson",
        numDeriv.method.args=list(eps=1e-4, d=0.0001, r=4, v=2, show.details=FALSE)) {

    # Compute the jacobian of the implicit function theta(s) estimated
    # by tsfa, i.e., compute the derivative d theta/ds'.
    # theta is the vector of parameters
    # s is the vector containing the unique elements of the sample covariance
    #   matrix
    # p is the number of factors
    # the following arguments are settings of tsfa, see estTSF,
    # the last 2 arguments are settings of numDeriv, see jacobian
    #
    # Typically, Tmat and BpermuteTarget are not given (although the latter
    # may be given in simulation studies) and factorNames are not relevant
    # here.

    if (is.null(BpermuteTarget)) {

        # Create a target, thereby also finding out whether an orthogonal
        # or oblique rotation has been performed.

        Sigma <- unvecs(s)

        #------- From estTSF ---------------------------------------

        z <- estFAmodel(Sigma, p=p, n.obs=n.obs, est="factanal",
               estArgs=estArgs, rotation=rotation, rotationArgs=rotationArgs,
               GPFargs=GPFargs, BpermuteTarget=BpermuteTarget,
               factorNames=factorNames)

        #cat("results of estFAmodel with eiv parameterization:\n\n")
        #print(z)

        #------- End adapted from estTSF.ML -----------------------------------

        ny    <- nrow(Sigma)
        npar  <- ny * p + ny              # loadings + uniquenesses
        if (!z$stats$orthogonal) {              # oblique rotation
            npar <- npar + (p*(p+1))/2    # factor covariances
        }

        if (!z$stats$estConverged || !z$stats$rotationConverged) {
           # No convergence: return missing and let the calling program
           # handle it from there.
           return(matrix(NA,nrow=npar,ncol=(ny*(ny+1))/2))
        }

        Bt <- z$loadings

    } else {

        Bt <- BpermuteTarget
    }

    # Return the numerically computed jacobian.

    return(jacobian(func=theta.tsfa, x=s, p=p, method=numDeriv.method,
             method.args=numDeriv.method.args, est=est,
             estArgs=estArgs, rotation=rotation, rotationArgs=rotationArgs,
             GPFargs=GPFargs, BpermuteTarget=Bt, factorNames=factorNames))
}

### This is the main routine; does all the estimation and computation
### of the basic statistics. Only for comparison of models and for
### computing derived fit indexes, other routines must be used.

tsfa.with.seasons <- function(y, p, nseason=1, diff.=TRUE,
        est="factanal",
        estArgs=list(scores="none", control=list(opt=list(maxit=10000))),
        rotation=if(p==1) "none" else "quartimin", rotationArgs=NULL,
        GPFargs=list(Tmat=diag(p),normalize=TRUE, eps=1e-5, maxit=1000),
        BpermuteTarget=NULL, factorNames=paste("Factor", seq(p)),
        numDeriv.method="Richardson",numDeriv.method.args=list()) {

    # Compute estimate of asymptotic covariance matrix of the parameters
    # of a TSFA analysis, assuming seasonal data, although nseason=1
    # can be used for yearly data or if no seasonal effects are present.
    #
    # y = raw data, p = number of factors, nseason = number of seasons
    # diff. = indicator whether analysis is done on differenced data
    # or not. The other arguments are arguments of the function
    # jacobian.tsfa.

    # In addition, compute estimate of asymptotic covariance matrix of the
    # residuals (s - \hat{\sigma}), which is used for testing and fit
    # statistics, and compute chi-square statistics.

    # Compute the sample covariance matrix.

    zz    <- if (diff.) diff(y) else y
    n.obs <- periods(zz)
    SmplCov <- cov(zz)

    # Store it in vector s. (I think jacobian wants a vector, not a matrix
    # with one column.)

    s <- c(vecs(SmplCov))

    # Compute asymptotic covariance matrix of s.

    Upsilon <- HAC.of.s.with.seasons(y, nseason, diff.)

    # Compute the estimates. This is done a couple of times too often,
    # but that will not be important relatively, given the number of
    # times the estimators must be computed for the numerical derivatives.

    #------- From estTSF ---------------------------------------

    #z <- estFAmodel(Sigma, p, n.obs=n.obs, est="factanal",
    z <- estFAmodel(SmplCov, p, n.obs=n.obs, est="factanal",
           estArgs=estArgs, rotation=rotation, rotationArgs=rotationArgs,
           GPFargs=GPFargs, BpermuteTarget=BpermuteTarget,
           factorNames=factorNames)

    #------- End adapted from estTSF.ML -----------------------------------

    B     <- z$loadings
    omega <- diag(z$Omega)
    if (z$stats$orthogonal) { # orthogonal rotation
        theta <- mat.to.theta(loadings=B, uniquenesses=omega)
    } else { # oblique rotation
        Phi   <- z$Phi
        theta <- mat.to.theta(loadings=B, uniquenesses=omega, Phi=Phi)
    }

    # Compute the vector sigma of unique elements of the model-implied
    # covariance matrix Sigma.

    sig <- sigma.tsfa(theta,p,orthogonal=z$stats$orthogonal)

    if (!z$stats$estConverged || !z$stats$rotationConverged) {
        # No convergence: return whatever we currently have and a lot of
        # missings and let the calling program handle it from there.
        # We can do a little bit better than this, e.g. df can be computed.
        # Implement this later.

        if (z$stats$orthogonal) {
            return(list(loadings=B,uniquenesses=omega,theta=theta,
                Acov.theta=NULL,s=s,sigma=sig,chisq=NULL,df=NULL,
                chisq.null=NULL,df.null=NULL,orthogonal=TRUE,
                estConverged=z$stats$estConverged,
                rotationConverged=z$stats$rotationConverged,jacConverged=NULL))
        } else {
            return(list(loadings=B,Phi=Phi,uniquenesses=omega,theta=theta,
                Acov.theta=NULL,s=s,sigma=sig,chisq=NULL,df=NULL,
                chisq.null=NULL,df.null=NULL,orthogonal=FALSE,
                estConverged=z$stats$estConverged,
                rotationConverged=z$stats$rotationConverged,jacConverged=NULL))
        }
    }

    # Compute the jacobian J = d theta/ds'.

    J <- jacobian.tsfa(s=s, p=p, est=est, estArgs=estArgs,
            rotation=rotation, rotationArgs=rotationArgs,
            GPFargs=GPFargs, BpermuteTarget=BpermuteTarget,
            factorNames=factorNames, numDeriv.method=numDeriv.method,
            numDeriv.method.args=numDeriv.method.args)

    if (any(!is.finite(c(J)))) {

        # Infinite or missing values, so no proper convergence. Return
        # whatever we currently have and a lot of missings and let the
        # calling program handle it from there. We can do a little bit
        # better than this, e.g. df can be computed. Implement this later.

        if (z$stats$orthogonal) {
            return(list(loadings=B,uniquenesses=omega,theta=theta,
                Acov.theta=NULL,s=s,sigma=sig,chisq=NULL,df=NULL,
                chisq.null=NULL,df.null=NULL,orthogonal=TRUE,
                estConverged=z$stats$estConverged,
                rotationConverged=z$stats$rotationConverged,jacConverged=FALSE))
        } else {
            return(list(loadings=B,Phi=Phi,uniquenesses=omega,theta=theta,
                Acov.theta=NULL,s=s,sigma=sig,chisq=NULL,df=NULL,
                chisq.null=NULL,df.null=NULL,orthogonal=FALSE,
                estConverged=z$stats$estConverged,
                rotationConverged=z$stats$rotationConverged,jacConverged=FALSE))
        }
    }

    # Asymptotic covariance matrix of theta. Don't forget the n.obs!
    # Leave t- and p-values to calling program.

    ACov.theta <- J %*% Upsilon %*% t(J)
    se.theta   <- sqrt(diag(ACov.theta)/n.obs)

    # and its jacobian Delta = d sigma/d theta'.
    # This is easy analytics.

    Delta <- jacobian.sigma.tsfa(theta, p, orthogonal=z$stats$orthogonal)

    # Compute the asymptotic covariance matrix of s - sigma.

    M <- NCOL(y)  #  PG added
    Mstar      <- (M*(M+1))/2
    IM2        <- diag(1,Mstar)
    Acov.s.sig <- (IM2 - Delta %*% J) %*% Upsilon %*% t(IM2 - Delta %*% J)

    # Rank of this Acov = degrees of freedom.

    df <- Mstar - M*p - M + (p*(p-1))/2

    # Generalized inverse of this Acov.

    AsmsInv <- MPinverse(Acov.s.sig, rank=df)

    # Compute the chi-square.
    # Even though s and sig are vectors, this should work.
    # Maybe consider replacing n.obs by n.obs-1, as usual in SEM,
    # or the Bartlett-correction
    #   n.obs - 1 - (2M + 5)/6
    # sometimes advocated for CFA, or the original Bartlett correction
    #   n.obs - 1 - (2M + 5)/6 - 2p/3
    # advocated for EFA and used, e.g., in SPSS FACTOR.
    # We could make this an option. The latter Bartlett correction may not
    # work well when sequential chi-square tests are computed or
    # for comparing fit according to fit indexes, because it depends on
    # the number of factors, which is not constant across models.

    chisq <- n.obs * t(s-sig) %*% AsmsInv %*% (s-sig)

    ####### I think it's better to move the null model stuff to the top!!!

    # chi-square for the null model.
    # The null model is the independence model: 0 factors, only a
    # diagonal error covariance matrix.

    # For ML estimation, its estimate is the diagonal
    # matrix equal to diag(S). Hence, S - Sigma = S with its diagonal
    # set to zero.

    if (est=="factanal") {
        diag(SmplCov) <- 0
        s.null        <- c(vecs(SmplCov))
        chisq.null    <- n.obs * t(s.null) %*% AsmsInv %*% s.null
    } else {
        warning(paste("Estimation method ",est,
                   " not implemented yet for null model.",sep=""))
        s.null        <- rep(NA,Mstar)
        chisq.null    <- NULL
    }
    df.null <- Mstar - M

    # Make a distinction between orthogonal and oblique.

Acov.theta <- ACov.theta   # PG added

    if (z$stats$orthogonal) {
        return(list(loadings=B,uniquenesses=omega,theta=theta,
                Acov.theta=Acov.theta,s=s,sigma=sig,chisq=chisq,df=df,
                chisq.null=chisq.null,df.null=df.null,orthogonal=TRUE,
                estConverged=TRUE, rotationConverged=TRUE,jacConverged=TRUE))
    } else {
        return(list(loadings=B,Phi=Phi,uniquenesses=omega,theta=theta,
                Acov.theta=Acov.theta,s=s,sigma=sig,chisq=chisq,df=df,
                chisq.null=chisq.null,df.null=df.null,orthogonal=FALSE,
                estConverged=TRUE, rotationConverged=TRUE,jacConverged=TRUE))
    }
}

#-------------------------------------------------------------------------------

jacobian.theta.eiv <- function(cov, nobs=NULL, loadings, uniquenesses, Phi,
                           identity=1:ncol(loadings)) {

    # Jacobian d theta/d s', where theta(-hat) is an implicit function of s.
    # theta = vector of parameters in the eiv parameterization of a factor
    #         analysis model
    # s = vecs(S) = vector of nonduplicated sample (co)variances.

    # From the implicit function theorem:
    # d theta/d s' = - (d^2 L/d theta d theta')^{-1} (d^2 L/d theta d s')

    # Htt = d^2 L/d theta d theta'
    Htt <- hessian.FA.eiv(cov, nobs, loadings, uniquenesses, Phi, identity)

    # Hts = d^2 L/d theta d s'
    Hts <- cross.hessian.FA.eiv(nobs, loadings, uniquenesses, Phi, identity)

    # Note that nobs drops out, so is not really necessary.
    return(-solve(Htt,Hts))
}

################  from cov03#########
# Some general utility functions.

const.pi <- function() {
    # Safer to have pi as a function, so that it can be used to denote something
    # else as a variable.
    #
    # Just typing the constant to 20 decimals is faster.
    # From Wikipedia, but the first 16 decimals coincide with what R gives
    # (which is only 16 decimals).
    return(3.14159265358979323846) # 26433 83279 50288 41971 69399 37510
}

tr <- function(A) {

  # trace of a matrix (generally believed to be square)
  # I still find it strange that R does not provide such a function.

  return(sum(diag(A)))
}

duplic.times <- function(A, B=NULL) {

    # This function computes D_M' (A %x% B), where
    # A is an M x p matrix and B is an M x q matrix, and
    # D_M is the duplication matrix of order M.

    if (is.null(B)) B <- A
    A  <- as.matrix(A)
    B  <- as.matrix(B)
    M  <- nrow(A)
    if (nrow(B) != M) stop("numbers of rows of A and B must be equal.")
    if (ncol(B) == 1) { # swap(A,B)
        C <- B
        B <- A
        A <- C
    }

    C <- NULL
    for (j in 1:M) {
        # It appears that we do not need to use drop=FALSE.
        C <- rbind(C, A[j,] %x% B[j,])
        if (j < M) {
            for (i in (j+1):M) {
                # It appears that we do not need to use drop=FALSE.
                C <- rbind(C, (A[i,] %x% B[j,] + A[j,] %x% B[i,]))
            }
        }
    }

    return(C)
}

duplicMPtimes <- function(A, B=NULL) {

    # This function computes D_M^+ (A %x% B), where
    # A is an M x p matrix and B is an M x q matrix, and
    # D_M^+ is the Moore-Penrose inverse of the duplication matrix
    # of order M.
    #
    # If A = B and p = q = 1, then this simply computes a column vector
    # (matrix with 1 column) with all nonduplicated elements A_i A_j.
    # If B = NULL, B = A. Furthermore, if the number of columns
    # of either A or B is 1, then D_M^+ (A %x% B) = D_M^+ (B %x% A),
    # which is used below.

    if (is.null(B)) B <- A
    A  <- as.matrix(A)
    B  <- as.matrix(B)
    M  <- nrow(A)
    if (nrow(B) != M) stop("numbers of rows of A and B must be equal.")
    if (ncol(B) == 1) { # swap(A,B)
        C <- B
        B <- A
        A <- C
    }

    C <- NULL
    if (ncol(A) == 1) {
        for (j in 1:M) {
            for (i in j:M) {
                # (j,j) is automatically handled well by this formula.
                # It appears that we do not need to use drop=FALSE.
                C <- rbind(C, (A[i,1] * B[j,] + A[j,1] * B[i,])/2)
            }
        }
    } else {
        for (j in 1:M) {
            for (i in j:M) {
                # (j,j) is automatically handled well by this formula.
                # It appears that we do not need to use drop=FALSE.
                C <- rbind(C, (A[i,] %x% B[j,] + A[j,] %x% B[i,])/2)
            }
        }
    }

    return(C)
}

dup.kron.dup <- function(A, B=NULL) {

    # This function computes D_M' (A %x% B) D_k, where
    # A and B are M x k matrices, D_M and D_k are duplication matrices.
    #
    # If B = NULL, B = A.

    if (is.null(B)) B <- A
    A  <- as.matrix(A)
    B  <- as.matrix(B)
    M  <- nrow(A)
    k  <- ncol(A)
    if (nrow(B) != M) stop("numbers of rows of A and B must be equal.")
    if (ncol(B) != k) stop("numbers of columns of A and B must be equal.")

    C <- NULL
    for (j in 1:k) {
        # It is now more efficient to single out (j,j).

        C <- cbind(C, duplic.times(A[,j],B[,j]))

        if (j < k) {
            for (i in (j+1):k) {
                C <- cbind(C, duplic.times(A[,j],B[,i])
                              + duplic.times(A[,i],B[,j]) )
            }
        }
    }

    return(C)
}

dMP.kron.dMP <- function(A, B=NULL) {

    # This function computes D_M^+ (A %x% B) (D_M^+)', where
    # A and B are M x M matrices and D_M^+ is the Moore-Penrose inverse
    # of the duplication matrix of order M.
    #
    # If B = NULL, B = A.

    if (is.null(B)) B <- A
    A  <- as.matrix(A)
    B  <- as.matrix(B)
    M  <- nrow(A)
    if (nrow(B) != M) stop("numbers of rows of A and B must be equal.")
    if (ncol(A) != M) stop("A is not square.")
    if (ncol(B) != M) stop("B is not square.")

    C <- NULL
    for (j in 1:M) {
        # It is now more efficient to single out (j,j).

        C <- cbind(C, duplicMPtimes(A[,j],B[,j]))

        if (j < M) {
            for (i in (j+1):M) {
                C <- cbind(C, (duplicMPtimes(A[,j],B[,i])
                              + duplicMPtimes(A[,i],B[,j]))/2 )
            }
        }
    }

    return(C)
}

dMP.kron.d <- function(A, B=NULL) {

    # This function computes D_M^+ (A %x% B) D_k, where
    # A and B are M x k matrices, D_M^+ is the Moore-Penrose inverse
    # of the duplication matrix of order M, and D_k is the duplication
    # matrix of order k.
    #
    # If B = NULL, B = A.

    if (is.null(B)) B <- A
    A  <- as.matrix(A)
    B  <- as.matrix(B)
    M  <- nrow(A)
    k  <- ncol(A)
    if (nrow(B) != M) stop("numbers of rows of A and B must be equal.")
    if (ncol(B) != k) stop("numbers of columns of A and B must be equal.")

    C <- NULL
    for (j in 1:k) {
        # It is now more efficient to single out (j,j).

        C <- cbind(C, duplicMPtimes(A[,j],B[,j]))

        if (j < k) {
            for (i in (j+1):k) {
                C <- cbind(C, duplicMPtimes(A[,j],B[,i])
                              + duplicMPtimes(A[,i],B[,j]) )
            }
        }
    }

    return(C)
}

vecs <- function(A) {

  # "symmetric" vec of a matrix: only the nonduplicated elements
  # This does not check for symmetry, it assumes it and uses it.

  vv <- NULL
  n  <- ncol(A)
  for (j in 1:n) {
    vv <- c(vv, A[j:n,j])
  }

  return(matrix(vv,ncol=1))
}

unvecs <- function(v) {

  # Compute symmetric matrix A from its nonduplicated elements
  # v = "symmetric" vec of a matrix A: only the nonduplicated elements

  v <- as.matrix(v) # Column vector = matrix with 1 column
  p <- nrow(v)      # Number of elements

  # Determine n = number of rows (and columns) of A.
  # p = n(n+1)/2  <=>  n = sqrt(2p + 0.25) - 0.5

  n = floor(sqrt(2 * p + 0.25)) # To be sure it's exactly the desired integer

  # Now process the vector

  m <- 1
  A <- matrix(nrow=n, ncol=n)
  for (j in 1:n) {
    A[j:n,j] <- v[m:(m+n-j),1,drop=FALSE]
    A[j,j:n] <- t(A[j:n,j,drop=FALSE])
    m <- m + n - j + 1
  }

  return(A)
}

strictlower <- function(A) {

  # Return the below-diagonal elements of a (typically square) matrix
  # in a vector.

  p <- nrow(A)
  q <- ncol(A)
  if (q >= p) q <- p-1
  if (q < 1) return(NULL)

  v <- NULL
  for (j in 1:q) {
      v <- c(v, A[(j+1):p,j])
  }

  return(v)
}

corr.to.corrmat <- function(v) {

  # Inverse of strictlower, but more restricted: It is assumed that
  # the original matrix was a square correlation matrix, so 1's are
  # added on the diagonal, and the elements above the diagonal are
  # copied from the ones below the diagonal.

  v <- as.matrix(v) # Column vector = matrix with 1 column
  p <- nrow(v)      # Number of elements

  # Determine n = number of rows (and columns) of A.
  # p = n(n-1)/2  <=>  n = sqrt(2p + 0.25) + 0.5

  n = ceiling(sqrt(2 * p + 0.25)) # To be sure it's exactly the desired integer

  # Now process the vector

  m <- 1
  A <- diag(1,n) # start with identity
  for (j in 1:(n-1)) {
    A[(j+1):n,j] <- v[m:(m+n-j-1),1,drop=FALSE]
    A[j,(j+1):n] <- t(A[(j+1):n,j,drop=FALSE])
    m <- m + n - j
  }

  return(A)
}

MPinverse <- function(X, rank, tol) {

    # Compute the Moore-Penrose generalized inverse of X.
    # Adapted from ginv, which only allows specifying a tolerance
    # (tol) for determining when a singular value should be treated
    # as zero. However, in our case, we already know the exact rank
    # that the matrix must have, so this function uses the rank
    # given by the user. If no rank is given, ginv is called.

    if(is.null(rank)) {
        # No explicit rank given, use ginv. This should not be necessary
        # for our usage, but you never know what happens with your
        # routines later.
        if(require("MASS")) return(ginv(X,tol)) else stop("package MASS not found.")
    }

    if (length(dim(X)) > 2 || !is.numeric(X))
        stop("'X' must be a numeric matrix")
    if (!is.matrix(X))
        X <- as.matrix(X)
    Xsvd <- svd(X)

    # Some experimentation, conventional definitions of the svd,
    # and the source code for the `ginv' function suggest that the
    # singular values are sorted in descending order, although this
    # is not explicitly mentioned in the help or (equivalently)
    # the reference manual.
    # The definition of Positive is the only difference with ginv.
    # It dependes critically on this assumed ordering of the
    # singular values. Hence the following check.

    increasing <- Xsvd$d[1:(length(Xsvd$d)-1)] < Xsvd$d[2:length(Xsvd$d)]
    if(any(increasing)) {
        stop(paste("Singular values are not ordered descendingly. \n",
             "Please notify the authors of the tsfa package."))
    }

    # Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0) # ginv
    Positive <- seq(Xsvd$d) <= rank                # our change

    if (all(Positive))
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive))
        array(0, dim(X)[2:1])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
        t(Xsvd$u[, Positive, drop = FALSE]))
}

################  from cov04#########
#----------------------------------------------------------------------
# Some functions to squeeze a multivariate mean in the right form
# for the sandwich package.

mvmean <- function(y) {

    # Multivariate mean. To trick the sandwich package.

    mm <- colMeans(y)

    # Residuals; from tsfa.R

    ee <- sweep(y,2,mm, "-")

    z <- list(coefficients=mm, data=y)
    class(z) <- "mvmean"

    return(z)
}

coef.mvmean <- function(object) {

    # Extract coefficients from an "mvmean" object.

    return(object$coefficients)
}

fitted.mvmean <- function(object) {

    # Compute fitted values from an "mvmean" object.

    mu    <- coefficients(object)
    y     <- object$data  # Is there no generic accessor method for this?
    n.obs <- nrow(y)
    ny    <- ncol(y)

    return(t(matrix(mu,nrow=ny,ncol=n.obs)))
}

residuals.mvmean <- function(object) {

    # Extract residuals from an "mvmean" object.

    return(object$data - fitted(object))
}

estfun.mvmean <- function(object) {

    # Estimating function of an "mvmean" object, for sandwich.
    # Note that the estfun functions implemented in sandwich,
    # for usage with the "lm" object class, do not work for
    # multivariate dependent data.

    return(resid(object))
}

#----------------------------------------------------------------------
# Two functions for conversion from parameter matrices to parameter
# vector and vice versa.

mat.to.theta <- function(loadings, uniquenesses, Phi=NULL) {

    # Collect all estimated parameters from the parameter matrices
    # B, Phi, and omega, from a (TS)FA model in one parameter vector
    # theta.

    if (is.null(Phi)) {
        # Orthogonal solution, no Phi.
        return(c(loadings, uniquenesses))
    } else {
        # Oblique solution, with Phi.
        return(c(loadings, vecs(Phi), uniquenesses))
    }
}

theta.to.mat <- function(theta,p,orthogonal=FALSE) {

    # Assemble the parameter matrices B, Phi, and omega, from a (TS)FA
    # model from the parameter vector theta that contains their elements.
    # p = number of factors
    # orthogonal = indicates whether an orthogonal solution is in theta,
    #              in which case there is no Phi, or an oblique solution,
    #              in which case there is a Phi as well.

    # I don't check consistency here.

    n <- length(theta)
    if (orthogonal) {
        # Orthogonal solution, no Phi.

        # Number of indicators. This'd better be an integer!
        M <- n/(p+1)

        B <- matrix(theta[1:(M*p)],M,p)
        omega <- theta[(M*p + 1):n]
        return(list(loadings=B,uniquenesses=omega))
    } else {
        # Oblique solution, with Phi.
        nphi <- (p*(p+1))/2

        # Number of indicators. Again, let's hope that this is an integer.
        M <- (n-nphi)/(p+1)

        B     <- matrix(theta[1:(M*p)],M,p)
        Phi   <- unvecs(theta[(M*p+1):(M*p+nphi)])
        omega <- theta[(M*p + nphi + 1):n]
        return(list(loadings=B,Phi=Phi,uniquenesses=omega))
    }
}

#----------------------------------------------------------------------
# sigma as a function of theta.

sigma.tsfa <- function(theta, p, orthogonal=FALSE) {

    # Compute sigma(theta), where
    # sigma = vecs(Sigma) is the vector containing the nonduplicated
    #   elements of the model-implied covariance matrix Sigma of a
    #   TSFA model, and theta is the vector containing the parameters.
    # Sigma = B Phi B' + diag(omega) and
    # theta contains the elements of B, Phi, and omega.

    # First, extract the parameter matrices from the vector theta.

    zz    <- theta.to.mat(theta, p, orthogonal)

    # Second, copy them to simpler symbols.

    B     <- zz$loadings
    omega <- zz$uniquenesses

    # Depending on orthogonal, return vecs(Sigma).

    if (!orthogonal) {
        return(vecs(B %*% zz$Phi %*% t(B) + diag(omega)))
    } else {
        return(vecs(tcrossprod(B) + diag(omega)))
    }
}

# The Jacobian of sigma as a function of theta.

jacobian.sigma.tsfa <- function(theta, p, orthogonal=FALSE) {

    # Compute the jacobian Delta = d sigma/d theta, where
    # sigma = vecs(Sigma) is the vector containing the nonduplicated
    #   elements of the model-implied covariance matrix Sigma of a
    #   TSFA model, and theta is the vector containing the parameters.
    # Sigma = B Phi B' + diag(omega) and
    # theta contains the elements of B, Phi, and omega.

    # The formula is an adaptation of formula (8.1) on p. 190 in
    # Wansbeek & Meijer (2000):
    #   Delta = D_M^+ (2 B Phi %x% I_M), (B %x% B)D_k, H_M),
    # where D_k is the duplication matrix, D_M^+ = (D_M'D_M)^{-1} D_M',
    # and H_M is the diagonalization matrix.

    # First, extract the parameter matrices from the vector theta.

    zz  <- theta.to.mat(theta, p, orthogonal)
    B   <- zz$loadings
    if (!orthogonal) {
        Phi <- zz$Phi
    } else {
        Phi <- diag(1, p)
    }
    # We don't need omega.

    # Identity matrix of order M.

    M  <- nrow(B)
    IM <- diag(1, M)

    # Derivatives with respect to vec(B).

    Delta <- duplicMPtimes(2 * B %*% Phi, IM)

    # Add derivatives with respect to vecs(Phi).

    if (!orthogonal) Delta <- cbind(Delta, dMP.kron.d(B))

    # Add derivatives with respect to omega.

    for(j in 1:M) {
        Delta <- cbind(Delta, duplicMPtimes(IM[,j]))
    }

    return(Delta)
}

hessian.sigma.tsfa <- function(theta, p, orthogonal=FALSE) {

    # Compute the 3-way array with elements = d^2 sigma_m/d theta d theta', where
    # sigma = vecs(Sigma) is the vector containing the nonduplicated
    #   elements of the model-implied covariance matrix Sigma of a
    #   TSFA model, and theta is the vector containing the parameters.
    # Sigma = B Phi B' + diag(omega) and
    # theta contains the elements of B, Phi, and omega.

    # Formulas:
    # d^2 sigma_jj / dB_jr dB_js   = 2 Phi_rs
    # d^2 sigma_jj / dB_jr dPhi_rs = 2 B_js
    # d^2 sigma_ij / dB_ir dB_js   =   Phi_rs   (i != j)
    # d^2 sigma_ij / dB_ir dPhi_rs =   B_js     (i != j)
    # d^2 sigma_ij / dB_jr dPhi_rs =   B_is     (i != j)

    # First, extract the parameter matrices from the vector theta.

    zz  <- theta.to.mat(theta, p, orthogonal)
    B   <- zz$loadings
    if (!orthogonal) {
        Phi <- zz$Phi
    } else {
        Phi <- diag(1, p)
    }
    # We don't need omega.

    # initialize array of 2nd derivatives

    M  <- nrow(B)         # number of indicators
    Mstar <- (M*(M+1))/2  # number of nonduplicated elements of covariance matrix
    npar  <- M * p + M
    if (!orthogonal) npar <- npar + (p*(p+1))/2
    D <- array(0,dim=c(npar,npar,Mstar))

    # loop over elements of sigma

    t <- 0
    for (j in 1:M) {
        t <- t + 1   # index of sigma_jj
        w <- M * p   # index of Phi_rs
        for (s in 1:p) {
            v <- (s-1) * M + j  # index of B_js
            for (r in s:p) {
                w <- w + 1          # index of Phi_rs
                u <- (r-1) * M + j  # index of B_jr

                # 2nd derivative of sigma_jj w.r.t. B_jr and B_js
                D[u,v,t] <- 2 * Phi[r,s]
                D[v,u,t] <- D[u,v,t]

                if (!orthogonal) {
                    # 2nd derivative of sigma_jj w.r.t. B_jr and Phi_rs
                    D[u,w,t] <- 2 * B[j,s]
                    D[w,u,t] <- D[u,w,t]

                    # 2nd derivative of sigma_jj w.r.t. B_js and Phi_rs
                    D[v,w,t] <- 2 * B[j,r]
                    D[w,v,t] <- D[v,w,t]
                }
            }
        }

        if (j < M) { # This annoying clause is necessary because R treats (M+1):M as
                     # a downward stepping loop instead of an empty one.
            for (i in (j+1):M) {
                t <- t + 1
                w <- M * p   # index of Phi_rs
                for (s in 1:p) {
                    for (r in s:p) {
                        w <- w + 1          # index of Phi_rs

                        u <- (r-1) * M + i  # index of B_ir
                        v <- (s-1) * M + j  # index of B_js

                        # 2nd derivative of sigma_ij w.r.t. B_ir and B_js
                        D[u,v,t] <- Phi[r,s]
                        D[v,u,t] <- D[u,v,t]

                        u <- (r-1) * M + j  # index of B_jr
                        v <- (s-1) * M + i  # index of B_is

                        # 2nd derivative of sigma_ij w.r.t. B_jr and B_is
                        D[u,v,t] <- Phi[r,s]
                        D[v,u,t] <- D[u,v,t]

                        if (!orthogonal) {
                            # 2nd derivative of sigma_ij w.r.t. B_jr and Phi_rs
                            D[u,w,t] <- B[i,s]
                            D[w,u,t] <- D[u,w,t]

                            # 2nd derivative of sigma_ij w.r.t. B_is and Phi_rs
                            D[v,w,t] <- B[j,r]
                            D[w,v,t] <- D[v,w,t]

                            u <- (r-1) * M + i  # index of B_ir
                            v <- (s-1) * M + j  # index of B_js

                            # 2nd derivative of sigma_ij w.r.t. B_ir and Phi_rs
                            D[u,w,t] <- B[j,s]
                            D[w,u,t] <- D[u,w,t]

                            # 2nd derivative of sigma_ij w.r.t. B_js and Phi_rs
                            D[v,w,t] <- B[i,r]
                            D[w,v,t] <- D[v,w,t]
                        }
                    }
                }
            }
        }
    }

    return(D)
}

#----------------------------------------------------------------------

### Estimation of the parameters as a function of the covariance matrix.
### Why do we use `p', where factanal uses `factors' ?

theta.tsfa <- function(s, p, est="factanal",
        estArgs=list(scores="none", control=list(opt=list(maxit=10000))),
        rotation=if(p==1) "none" else "quartimin", rotationArgs=NULL,
        GPFargs=list(Tmat=diag(p),normalize=TRUE, eps=1e-5, maxit=1000),
        BpermuteTarget=NULL, factorNames=paste("Factor", seq(p))) {

    # This is a function that returns all TSFA parameters (B, Phi, omega)
    # stacked in the vector theta, as a function of the vector s that
    # contains all nonduplicated elements of the sample covariance
    # matrix S. It is essential here that it is the covariance matrix
    # and not the correlation matrix! It is also important that all options
    # like rotation and normalize are the same as used in the main estimation.
    # eps may need to be chosen considerably tighter here, though,
    # and BpermuteTarget will typically be the original estimate.
    # For the typical application of this function (compute numerical
    # derivatives), factorNames will not be important.

    # Not a good choice to use Sigma here, suggests population cov.mat.
    # instead of sample cov.mat., but I prefer to retain most of the
    # estTSF code unchanged.

    Sigma <- unvecs(s)

    #------- From estTSF ---------------------------------------

    z <- estFAmodel(Sigma, p, n.obs=NA, est="factanal",
           estArgs=estArgs, rotation=rotation, rotationArgs=rotationArgs,
           GPFargs=GPFargs, BpermuteTarget=BpermuteTarget,
           factorNames=factorNames)

    #------- End from estTSF -----------------------------------

    ny    <- nrow(Sigma)
    npar  <- ny * p + ny              # loadings + uniquenesses
    if (!z$stats$orthogonal) {        # oblique rotation
        npar <- npar + (p*(p+1))/2    # factor covariances
    }

    if (!(z$stats$estConverged) || !(z$stats$rotationConverged)) {
       # No convergence: return missing and let the calling program
       # handle it from there.
       return(rep(NA,npar))
    }

    hatOmega <- diag(z$Omega)
    # Perhaps better to use loadings(z). Is there a corresponding function
    # for Omega?
    hatB     <- z$loadings

    if (z$stats$orthogonal) { # orthogonal rotation
        return(mat.to.theta(loadings=hatB,uniquenesses=hatOmega))
    } else { # oblique rotation
        return(mat.to.theta(loadings=hatB,uniquenesses=hatOmega,Phi=z$Phi))
    }
}

#----------------------------------------------------------------------

data.for.HAC.with.seasons <- function(y, nseason=1, diff.=TRUE) {

    # Compile data for computing a heteroskedasticity and autocorrelation
    # consistent (HAC) estimator of the sample covariance matrix.
    # This new data matrix consists of the original data matrix y
    # with additional columns containing cross-products
    # y_i y_j, j=1,...,ncol(y), i=j,...,ncol(y).
    # If diff. == TRUE then the first differences of y are used instead of y.
    #
    # nseason = the number of seasons. If <= 1, then no seasons, and the
    # means and covariances are assumed to be the same at each time point.
    # If > 1, then means and covariances are to be computed for each
    # season separately. The return matrix zz is a matrix consisting of
    # dummy variables, with as many columns as there are seasons.

    yy    <- if (diff.) diff(y) else y
    n.obs <- periods(yy)
    n.var <- ncol(yy)
    if (nseason < 1) {
        k <- 1
    } else {
        k <- nseason
    }

    # Add the products of the columns as additional columns to yy.

    for (j in 1:n.var) {
        for (i in j:n.var) {
            yy <- cbind(yy, yy[,i] * yy[,j])
        }
    }

    # Compute the matrix zz of season dummies.
    # Probably R's standard time series facilities or Paul's
    # packages have a function season(y) or something, so that
    # this whole construction is unnecessary and inefficient, but
    # it is more efficient in terms of my own programming time.

    zz <- matrix(0, nrow=n.obs, ncol=k)
    n  <- 1
    s  <- 0
    while (n <= n.obs) {
        # s = season number
        s <- s + 1
        if (s > k) s <- s - k
        zz[n,s] <- 1
        n <- n + 1
    }

    return(list(yy=yy,zz=zz))
}

HAC.s.from.HAC.m <- function(A, mu) {

    # A = estimate of asymptotic covariance matrix of the sample mean of
    # m = (y', vecs(yy')')', and mu is the sample mean of y.
    #
    # Based on this, this function computes a HAC estimator of the
    # covariance matrix of s = vecs(S), using
    # the delta method and the explicit function s of the means m
    # of y: s_{ij} = sum(y_i y_j)/N - m_i m_j = m_k - m_i m_j
    # for some k.
    #
    # The order of the s_{ij} in the vector s is the same as the order
    # of the corresponding elements m_k in m, so we start with simply
    # copying the lower-right block of A.

    ny      <- length(mu)
    mu      <- matrix(mu, ncol=1)
    nystar  <- (ny*(ny+1))/2
    Upsilon <- A[(ny+1):(ny+nystar), (ny+1):(ny+nystar)]

    # Now make adjustments for the subtraction of the means.
    # Let D = ds/dm' = (D_1, I), where D_1 = ds/d mu'. Then
    # Upsilon = DAD' = A_{22} + D_1 A_{12} + A_{21} D_1' + D_1 A_{11} D_1'.
    # D_1 = - D_M^+ (I %x% mu + mu %x% I) = -2 D_M^+ (mu %x% I),
    # where D_M^+ is the Moore-Penrose inverse of the duplication matrix
    # of order M = ny.
    #
    # It follows that D_1 A_{12} = (A_{21} D_1')' = -2 D_M^+ (mu %x% A_{12})
    # and D_1 A_{11} D_1' = 4 D_M^+ (mu mu' %x% A_{11}) (D_M^+)'

    U2      <- 2 * duplicMPtimes(mu, A[1:ny, (ny+1):(ny+nystar)])
    Upsilon <- Upsilon - U2 - t(U2)
    Upsilon <- Upsilon + 4 * dMP.kron.dMP(tcrossprod(mu), A[1:ny, 1:ny])

    return(Upsilon)
}

#-------------------------------------------------------------------------------
# Some functions for ML using only the covariance matrix.
# Therefore, use "mwl" designator (for Maximum Wishart Likelihood), to distinguish
# from multivariate normal loglikelihood that includes the means as well.

loglike.mwl <- function(S, Sigma, nobs=100) {

    # Loglikelihood. Take this literally, so include the constants and do not
    # reverse. (Hence, this function must be maximized.)
    # The `100' default for nobs is only meant to give useful answers if the
    # number of observations is not given (or not relevant in case of artificial
    # tests).

    M <- nrow(S)
    L <- -0.5 * nobs * (M * log(2 * const.pi()) + log(det(Sigma)) + tr(solve(Sigma,S)))

    return(L)
}

gradient.mwl <- function(S, Sigma, nobs=100) {

    # Derivative of loglike.mwl w.r.t. sig = vecs(Sigma), i.e., the nonduplicated
    # elements of Sigma.
    # The `100' default for nobs is only meant to give useful answers if the
    # number of observations is not given (or not relevant in case of artificial
    # tests).

    iSig <- solve(Sigma)
    SES  <- iSig %*% (S - Sigma) %*% iSig

    # If SES is not symmetric, it's forced to be. (It should be, though.)

    g <- 0.5 * nobs * vecs(SES + t(SES) - diag(diag(SES)))

    # This is a matrix with 1 column.

    return(g)
}

hessian.mwl <- function(S, Sigma, nobs=100) {

    # Hessian matrix of loglike.mwl w.r.t. sig = vecs(Sigma),
    # i.e., the nonduplicated elements of Sigma.
    # Thus, it is the matrix d^2 L/d sig d sig'.
    # The `100' default for nobs is only meant to give useful answers if the
    # number of observations is not given (or not relevant in case of artificial
    # tests).

    iSig <- solve(Sigma)
    SES  <- iSig %*% (S - Sigma) %*% iSig
    HH   <- -nobs * dup.kron.dup(SES, iSig) - 0.5 * nobs * dup.kron.dup(iSig)

    return(HH)
}

info.mwl <- function(Sigma, nobs=100) {

    # Information matrix of loglike.mwl w.r.t. sig = vecs(Sigma),
    # i.e., the nonduplicated elements of Sigma.
    # Thus, it is the matrix -E d^2 L/d sig d sig'.
    # The `100' default for nobs is only meant to give useful answers if the
    # number of observations is not given (or not relevant in case of artificial
    # tests).

    iSig <- solve(Sigma)
    HH   <- 0.5 * nobs * dup.kron.dup(iSig)

    return(HH)
}

cross.hessian.mwl <- function(Sigma, nobs=100) {

    # cross-Hessian matrix of loglike.mwl w.r.t. sig = vecs(Sigma) and s = vecs(S).
    # Thus, it is the matrix d^2 L/d sig d s'.
    # The `100' default for nobs is only meant to give useful answers if the
    # number of observations is not given (or not relevant in case of artificial
    # tests).

    iSig <- solve(Sigma)
    HH   <- 0.5 * nobs * dup.kron.dup(iSig)

    return(HH)
}

#-------------------------------------------------------------------------------
# Loglikelihood and derivatives for generic FA model


loglike.FA <- function(cov, nobs=NULL, loadings, uniquenesses, Phi=NULL) {

    # Loglikelihood of a factor analysis model, using only the covariance matrix.
    # (If the intercepts are unrestricted, the means do not contribute to the
    # loglikelihood anyway.)
    #
    # cov = sample covariance matrix

    M <- nrow(cov)

    # Put the matrices in a parameter vector. (This is a bit inefficient,
    # because the reverse action is done in the sigma function, but it's
    # more convenient to use a tested function than write another one, however
    # similar.)

    theta <- mat.to.theta(loadings, uniquenesses, Phi)

    # Compute sig = vecs(Sigma).

    p          <- ncol(loadings)
    orthogonal <- FALSE
    if (is.null(Phi)) orthogonal <- TRUE
    sig <- sigma.tsfa(theta, p, orthogonal)

    # Return loglikelihood.

    return(loglike.mwl(cov, unvecs(sig), nobs))
}

gradient.FA <- function(cov, nobs=NULL, loadings, uniquenesses, Phi=NULL) {

    # Derivative of loglike.FA w.r.t. the elements of B, Phi, omega.

    # Put the matrices in a parameter vector. (This is a bit inefficient,
    # because the reverse action is done in the jacobian function, but it's
    # more convenient to use a tested function than write another one, however
    # similar.)

    theta <- mat.to.theta(loadings, uniquenesses, Phi)

    # Jacobian of sig = vecs(Sigma) w.r.t. theta

    p          <- ncol(loadings)
    orthogonal <- FALSE
    if (is.null(Phi)) orthogonal <- TRUE
    Jsig <- jacobian.sigma.tsfa(theta, p, orthogonal)

    # Compute sig = vecs(Sigma).

    sig <- sigma.tsfa(theta, p, orthogonal)

    # Derivatives of L w.r.t. sig.

    g <- gradient.mwl(cov, unvecs(sig), nobs)

    # Return the derivatives as a matrix with 1 column.

    return(crossprod(Jsig,g))
}

hessian.FA <- function(cov, nobs=NULL, loadings, uniquenesses, Phi=NULL) {

    # Hessian matrix of loglike.FA w.r.t. the elements of B, Phi, omega.

    # Formula:
    # d^2 L/d theta d theta' = J'(d^2 L/d sig d sig')J
    #                          + sum_{ij} (dL/d sig_ij)(d^2 sig_ij/d theta d theta')

    # Put the matrices in a parameter vector. (This is a bit inefficient,
    # because the reverse action is done in the jacobian function, but it's
    # more convenient to use a tested function than write another one, however
    # similar.)

    theta <- mat.to.theta(loadings, uniquenesses, Phi)
    p     <- ncol(loadings)
    orthogonal <- FALSE
    if (is.null(Phi)) orthogonal <- TRUE

    # Compute sig = vecs(Sigma) and its Jacobian w.r.t. theta.

    sig  <- sigma.tsfa(theta, p, orthogonal)
    Jsig <- jacobian.sigma.tsfa(theta, p, orthogonal)

    # 3-way array of 2nd derivatives of sigma w.r.t. theta
    D <- hessian.sigma.tsfa(theta, p, orthogonal)

    # Hessian matrix of L w.r.t. sig.

    Sigma <- unvecs(sig)
    Hss   <- hessian.mwl(cov, Sigma, nobs)

    # First term.

    H <- t(Jsig) %*% Hss %*% Jsig

    # Gradient of L w.r.t. sigma

    g <- gradient.mwl(cov, Sigma, nobs)

    # Add term sum_m g(m) d^2 sig(m)/d theta d theta'
    for (m in 1:length(sig)) {
        H <- H + g[m] * D[,,m]
    }

    # Return the Hessian matrix w.r.t. the FA parameter vector.

    return(H)
}

info.FA <- function(nobs=NULL, loadings, uniquenesses, Phi=NULL) {

    # Information matrix of loglike.FA w.r.t. the elements of B, Phi, omega.

    # Put the matrices in a parameter vector. (This is a bit inefficient,
    # because the reverse action is done in the jacobian function, but it's
    # more convenient to use a tested function than write another one, however
    # similar.)

    theta <- mat.to.theta(loadings, uniquenesses, Phi)

    # Jacobian of sig = vecs(Sigma) w.r.t. theta

    p          <- ncol(loadings)
    orthogonal <- FALSE
    if (is.null(Phi)) orthogonal <- TRUE
    Jsig <- jacobian.sigma.tsfa(theta, p, orthogonal)

    # Compute sig = vecs(Sigma).

    sig <- sigma.tsfa(theta, p, orthogonal)

    # Information matrix of L w.r.t. sig.

    H <- info.mwl(unvecs(sig), nobs)

    # Return the information matrix w.r.t. the FA parameter vector.

    return(t(Jsig) %*% H %*% Jsig)
}

cross.hessian.FA <- function(nobs=NULL, loadings, uniquenesses, Phi=NULL) {

    # cross-Hessian matrix of loglike.FA w.r.t. theta and s = vecs(S).
    # Thus, it is the matrix d^2 L/d theta d s' = J' d^2 L/d sig d s',
    # where J = d sig/d theta' and sig = vecs(Sigma).
    # The `100' default for nobs is only meant to give useful answers if the
    # number of observations is not given (or not relevant in case of artificial
    # tests).

    # Put the matrices in a parameter vector. (This is a bit inefficient,
    # because the reverse action is done in the jacobian function, but it's
    # more convenient to use a tested function than write another one, however
    # similar.)

    theta <- mat.to.theta(loadings, uniquenesses, Phi)

    # Jacobian of sig = vecs(Sigma) w.r.t. theta

    p          <- ncol(loadings)
    orthogonal <- FALSE
    if (is.null(Phi)) orthogonal <- TRUE
    Jsig <- jacobian.sigma.tsfa(theta, p, orthogonal)

    # Compute sig = vecs(Sigma).

    sig <- sigma.tsfa(theta, p, orthogonal)

    # cross-Hessian matrix of L w.r.t. sig and s'.

    Sigma <- unvecs(sig)
    H.sig.s <- cross.hessian.mwl(Sigma, nobs)

    # Return result

    return(crossprod(Jsig,H.sig.s))
}

#-------------------------------------------------------------------------------
# Functions for eiv parameterization

free.params.eiv <- function(M, p, identity=1:p) {

    # Create a vector with the indices of the free elements of the
    # eiv parameterization of a factor analysis model.
    # M = number of indicators
    # p = number of factors
    # identity = vector of indices of identity submatrix of loadings matrix

    npar <- M * p + M + (p*(p+1))/2
    vv   <- NULL
    for (j in 1:p) {
        # Parameter numbers of j-th colunm of B
        colj <- ((j-1)*M+1) : ((j-1) * M + M)
        # Exclude identity rows and add to index vector vv
        vv   <- c(vv,colj[-identity])
    }
    # Add remaining parameters
    vv <- c(vv,(M*p+1):npar)

    return(vv)
}

gradient.FA.eiv <- function(cov, nobs=NULL, loadings, uniquenesses, Phi,
                           identity=1:ncol(loadings)) {

    # gradient vector of loglike.FA w.r.t. the elements of B, Phi, omega
    # in the eiv parameterization, i.e., with Phi free and B having an I_p submatrix.

    # This is a subvector of the gradient.FA vector, with the elements
    # corresponding to the fixed elements of B removed.

    g <- gradient.FA(cov, nobs, loadings, uniquenesses, Phi)

    # Compute vector of indices of free parameters
    M  <- nrow(loadings)
    p  <- ncol(loadings)
    vv <- free.params.eiv(M, p, identity)

    # Return the gradient w.r.t. the eiv parameter vector.

    return(g[vv])
}

hessian.FA.eiv <- function(cov, nobs=NULL, loadings, uniquenesses, Phi,
                           identity=1:ncol(loadings)) {

    # Hessian matrix of loglike.FA w.r.t. the elements of B, Phi, omega
    # in the eiv parameterization, i.e., with Phi free and B having an I_p submatrix.

    # This is a submatrix of the hessian.FA matrix, with the rows and columns
    # corresponding to the fixed elements of B removed.

    H <- hessian.FA(cov, nobs, loadings, uniquenesses, Phi)

    # Compute vector of indices of free parameters
    M  <- nrow(loadings)
    p  <- ncol(loadings)
    vv <- free.params.eiv(M, p, identity)

    # Return the Hessian matrix w.r.t. the eiv parameter vector.

    return(H[vv,vv])
}

cross.hessian.FA.eiv <- function(nobs=NULL, loadings, uniquenesses, Phi,
                           identity=1:ncol(loadings)) {

    # cross-Hessian matrix of loglike.FA w.r.t. the elements of B, Phi, omega
    # in the eiv parameterization, i.e., with Phi free and B having an I_p submatrix,
    # in the rows, and w.r.t. s' = vecs(S)' in the columns.

    # This is a submatrix of the cross.hessian.FA matrix, with the rows
    # corresponding to the fixed elements of B removed.

    H <- cross.hessian.FA(nobs, loadings, uniquenesses, Phi)

    # Compute vector of indices of free parameters
    M  <- nrow(loadings)
    p  <- ncol(loadings)
    vv <- free.params.eiv(M, p, identity)

    # Return the Hessian matrix w.r.t. the eiv parameter vector and s.

    return(H[vv,])
}




#  initial scratchings for a dfs based loadings estimate

estGDdfa <- function(y, p, FAdiff=FALSE,  DFAdiff=FALSE,
                FA.est="factanal", 
		FA.estArgs=list(scores="none", control=list(opt=list(maxit=10000))),
                DFA.order=1, 
       DFA.MaxLik.algorithm="optim",
       DFA.MaxLik.algorithm.args=list(method="BFGS", upper=Inf, lower=-Inf,
                    hessian=TRUE, control=list(maxit=1000)),
		rotation=if(p==1) "none" else "quartimin", 
		rotationArgs=NULL, 
		GPFargs=list(Tmat=diag(p),normalize=TRUE, eps=1e-5, maxit=1000), 
		BpermuteTarget=NULL,
                factorNames=paste("Factor", seq(p))) {
  
  # 1/ est TSFA model as starting point
  # 2/ optimize as DFA
  # 3/ rotate

  # SS specifications, ( AR(1) )
  # state = z(t) = (xi_t, x_{t+1})' +1 or -1 ?
  #  H= [B | 0] or [0|B] , G=0 unless intercept
  
  #       AR | 0 
  # F=    ------
  #       I  | 0
  
  
  #  Q= svd(cov(Xi))
  #  R = svd(cov(y))  z0, ...rootP0=

   Step1 <- estTSFmodel(y, p, diff.=FAdiff, 
        FAest=FA.est, 
	FAestArgs=FA.estArgs,
        rotation="none")

   dfaModel <- SS( H=cbind(loadings(mStep1), matrix(0, p, p*AR.order)))

   Step2 <-  estMaxLik(dfaModel,
       TSdata(output=if(diff.DFA) diff(y) else y,  input=Xdata),
       algorithm=DFA.MaxLik.algorithm,
       algorithm.args=DFA.MaxLik.algorithm.args)
   
     B <- TSmodel(Step2)$H[,1:nfact, drop=FALSE]
   #  xihat <- smoother(Step2,
   #     TSdata(output=YY,input=Xdata))$smooth$state[,1:nfact, drop=FALSE]

   Step3 <-RotateFactorModel(B, Phi, Sigma,  uniquenesses,
        rotation=rotation, 
	rotationArgs=rotationArgs, 
	GPFargs=GPFargs,
	BpermuteTarget=BpermuteTarget,
        factorNames=factorNames)
   
   Step3$FAmodel  <- Step1
   Step3$DFAmodel <- Step2
   Step3
   }

# See estFAmodel call  possibly put this in estFAmodel
RotateFactorModel <- function(loadings, Phi, Sigma,  uniquenesses,
    	rotation= if(NCOL(B)==1) "none" else "quartimin", 
    	rotationArgs=NULL,
    	GPFargs=list(Tmat=diag(NCOL(B)), normalize=TRUE, eps=1e-5, maxit=1000),
    	BpermuteTarget=NULL,
    	factorNames=paste("Factor", seq(NCOL(B))),
    	indicatorNames=NULL ) {
   stds <- sqrt(diag(Sigma))
   if (rotation == "none") {
      loadings.std <- loadings	      
      Phi  <- NULL
      rotationConverged <- TRUE
      orthogonal <- TRUE
      }
   else {     
      rotB <- do.call(rotation, c(list(loadings), GPFargs, rotationArgs))
      loadings.std <-  rotB$loadings
      rotationConverged <- rotB$convergence
      orthogonal <- rotB$orthogonal
  
      # Make sure columns are ordered properly and have the correct signs.
      Phi <- rotB$Phi
      if (! is.null(BpermuteTarget)) {
	 z  <- permusign(diag(stds) %*% loadings.std, BpermuteTarget, Phi=Phi)
	 loadings.std <- diag(1/stds) %*% z$loadings
	 Phi  <- z$Phi
	 }
      }
   loadings <- diag(stds) %*% loadings.std
   dimnames(loadings) <- list(indicatorNames, factorNames)

   ### Compute Bartlett predictor 
   Sigma.std <- if(is.null(Phi)) loadings.std %*% t(loadings.std) + diag(uniquenesses) else
		     loadings.std %*% Phi %*% t(loadings.std) + diag(uniquenesses)
   SinvB <- solve(Sigma.std, loadings.std) 
   LB.std   <- solve(crossprod(loadings.std, SinvB), t(SinvB))

   LB <- LB.std %*% diag(1/stds)
   dimnames(LB) <- list(factorNames, indicatorNames)
   list(loadings=loadings,
        Omega=diag(stds * uniquenesses * stds), 
        Phi=Phi, LB=LB, LB.std=LB.std,  
        rotationStats=list(rotationConverged=rotationConverged,
	orthogonal=orthogonal, uniquenesses=uniquenesses, call=match.call()))
}

## This file contains functions to run dMLCA with patient-level covariates
## The code here are used in simulation.
## Sample size on each site are the same.
## Code for dMLCA with no patient-level covariates and different sample size on each site are 
## documented in a separate file, which is used in MIS-C subphenotyping

## Compile C functions:
## STEP 1: Delete the current 'Func.All.o' and 'Func.All.so' file, if any.
## STEP 2: In terminal, under the path of the file 'Func.All.c', run 'R CMD SHLIB Func.All.c' (no quote).
## STEP 3: In R, run
dyn.load('Func.All.so')

class_prior <- function(X, beta, alpha, K, n){
  temp1 <- X %*% beta + alpha[rep(seq_len(K), each = n), ]
  temp1 <- pmin(temp1, 7)
  temp1 <- pmax(temp1, -7)
  temp1 <- exp(temp1)
  temp1/as.vector(temp1 %*% rep(1,C))
}
## distributed latent class regression, with patient-level covariates
## Y: binary, each row is Y_ki = Y[(k-1)*n + i，] length q
## X: each row is X_ki = Y[(k-1)*n + i，] length p
## K: number of sites
## C: number of latent classes
## n: number of obs on each site (assume equal on different sites)
## p.initial: q*C
## alpha.initial: K*C, last column is 0
## beta.initial: p*C, last column is 0
DistLCR_NR <- function(Y, X, K, C, n, p.initial, alpha.initial, beta.initial, maxit=1000, tol=1e-6){
  
  p.last <- p.initial
  alpha.last <- alpha.initial
  beta.last <- beta.initial
  
  if(is.null(dim(X))) numX <- 1 else numX <- dim(X)[2]
  if(is.null(dim(beta.last))) beta.last <- matrix(beta.last,1)
  if(is.null(dim(alpha.last))) alpha.last <- t(alpha.last)
  
  conv <- FALSE
  iter1 <- 0
  diff_all <- c()
  for(iter in 1:maxit){
    iter1 <- iter1 + 1
    
    ## calculate prior class probability lambda
    lambda_y <- class_prior(X, beta.last, alpha.last, K, n)
    
    ## p_cj (c=1,...,C; j=1,...,q) (assume C=3; q=5)
    ## pdf.C requires the input probability matrix (P(Y|Z)) to be:
    ## 1-p_11, p_11; 1-p_21, p_21; 1-p_31, p_31; 1-p_12, p_12;...,1-p_32, p_32;...; 1-p_35, p_35; 
    p.last_vec <- as.vector(t(p.last))
    p.last.0.1 <- cbind(1-p.last_vec,p.last_vec)
    
    ## w_kic
    w <- postClass.C(lambda_y, p.last.0.1, Y+1, C, K)
    
    ## update p
    ## for Y with 1 dim (q=1)
    if(is.null(dim(Y))){
      numerator_p <- t(w * Y) %*% rep(1, K*n)
      denominator_p <- as.vector(t(w) %*% rep(1, K*n))
      p.current <- t(numerator_p/as.vector(denominator_p))
    } else
      ## for Binary y, the second column is P(Y=1)
      p.current_vec <- matrix(probHat.C(w, Y+1, C, K),ncol=2, byrow=TRUE)[,2]
    p.current <- matrix(p.current_vec, ncol=C, byrow=TRUE)
    
    ## update alpha, beta
    ## 1,2-order derivative of alpha_kc
    D.alpha <- d2lldalpha2.C(w, lambda_y, K=K)
    D.alpha.order.1 <- D.alpha$grad #site first, then class
    D.alpha.order.1 <- as.matrix(D.alpha.order.1) ##row: site
    D.alpha.order.2 <- D.alpha$hess
    ## 1,2-order order derivative of beta_c
    X_mat <- as.matrix(X)
    D.beta <- dLL2dBeta.C(w,lambda_y, X_mat, K)
    D.beta.order.1 <- D.beta$grad
    D.beta.order.2 <- D.beta$hess
    ## Gradient of beta and alpha
    D.alpha.order.1_mat <- matrix(D.alpha.order.1,K,C-1)
    G.alpha.beta <- c(D.beta.order.1, as.vector(t(D.alpha.order.1_mat))) ##beta: predictor then class, alpha: class than site
    ## Second order dirivative
    D.alpha.beta.order.2 <- d2lldab2.C(w, lambda_y, X_mat, K)
    
    ## Construct the Hessian matrix
    len.beta <- (C-1)*numX
    H.dim <- len.beta + (C-1)*K
    H.alpha.beta <- matrix(0, H.dim, H.dim)
    H.alpha.beta[1:len.beta, 1:len.beta] <- D.beta.order.2
    H.alpha.beta[(len.beta+1):H.dim, (len.beta+1):H.dim] <- D.alpha.order.2
    H.alpha.beta[1:len.beta, (len.beta+1):H.dim] <- D.alpha.beta.order.2
    H.alpha.beta[(len.beta+1):H.dim, 1:len.beta] <- t(D.alpha.beta.order.2)
    ## inverse of the Hessian
    H.inv <- ginv(H.alpha.beta)
    beta.alpha.last <- c(beta.last[,1:(C-1)], t(alpha.last[,1:(C-1)]))
    
    ## Newton-Raphson Step
    beta.alpha.current <- beta.alpha.last - H.inv %*% G.alpha.beta
    beta.current <- cbind(matrix(beta.alpha.current[1:len.beta], numX),rep(0, numX))
    alpha.current <- cbind(matrix(beta.alpha.current[(len.beta+1):H.dim], K, C-1, byrow=TRUE), rep(0, K))
    
    diff <- sqrt(sum((p.last-p.current)^2)+sum((beta.alpha.last-beta.alpha.current)^2))
    diff_all <- c(diff_all, diff)
    
    beta.last <- beta.current
    alpha.last <- alpha.current
    p.last <- p.current
    if(diff <= tol) {conv = TRUE; break}
  }
  ## Calculate the optimal loglikelihood
  lambda_y <- class_prior(X, beta.last, alpha.last, K, n)
  ## p_cj (c=1,...,C; j=1,...,q) (assume C=3; q=5)
  ## pdf.C requires the input probability matrix (P(Y|Z)) to be:
  ## 1-p_11, p_11; 1-p_21, p_21; 1-p_31, p_31; 1-p_12, p_12;...,1-p_32, p_32;...; 1-p_35, p_35; 
  p.last_vec <- as.vector(t(p.last))
  p.last.0.1 <- cbind(1-p.last_vec,p.last_vec)
  ## All .C function need Y=1,2,...
  pdf.y <- pdf.C(p.last.0.1, Y+1, C, K)
  ## loglikelihood(Y)
  loglike_ki <- log(as.vector((pdf.y * lambda_y) %*% rep(1,C)))
  loglike <- sum(loglike_ki)
  list(phat = p.last, alphahat = alpha.last, betahat = beta.last,
       lambda_y = lambda_y, loglikelihood = loglike,
       niter = iter1, conv = conv,diff=diff_all)
}

## main function for distributed latent class regression with patient-level covariates
## initial = "random" or "manual" (given in p.initial, alpha.initial, beta.initial)
## rep: how many times random initial
DistLCR <- function(Y, X, K, C, n, p.initial=NULL, alpha.initial=NULL, beta.initial=NULL, 
                    maxit=1000, tol=1e-3, initial = "random", rep = 5){
  if(is.null(dim(X))) numX <- 1 else numX <- dim(X)[2]
  if(is.null(dim(Y))) q <- 1 else q <- dim(Y)[2]
  if(initial == "random"){
    llik <- c()
    loglik <- -Inf
    for (s in 1:rep){
      #cat("s=",s,'\n')
      converge <- FALSE
      count.try <- 0
      while(!converge & count.try <= 5){
        count.try <- count.try+1
        alpha.initial <- cbind(matrix(rnorm(K*(C-1)),K,C-1), rep(0,K))
        beta.initial <- cbind(matrix(rnorm((C-1)*numX),numX,C-1),rep(0,numX))
        p.initial <- matrix(runif(q*C, 0.1, 0.9),q,C)
        
        fit <- -1; no.try <- 0
        while(is.numeric(fit) & no.try <= 5){
          no.try <- no.try + 1
          fit <- tryCatch(DistLCR_NR(Y,X,K,C,n, p.initial, alpha.initial, beta.initial, tol=tol, maxit=maxit), 
                          error=function(o) -1)
          if(is.numeric(fit)){
            alpha.initial <- cbind(matrix(rnorm(K*(C-1)),K,C-1), rep(0,K))
            beta.initial <- cbind(matrix(rnorm((C-1)*numX),numX,C-1),rep(0,numX))
            p.initial <- matrix(runif(C, 0.1, 0.9),q,C)
          }}
        
        if(!is.numeric(fit)) converge <- fit$conv
      }
      if(!converge) next
      llik <- c(llik, fit$loglikelihood)
      if(fit$loglikelihood > loglik) {loglik <- fit$loglikelihood
      alpha.initial.best <- alpha.initial; beta.initial.best <- beta.initial;
      p.initial.best <- p.initial; fit.best <- fit}
    }
    alpha.initial <- alpha.initial.best
    beta.initial <- beta.initial.best
    p.initial <- p.initial.best
  }
  if(initial == "manual"){
    if (method == "NR"){
    fit <- DistLCR_NR(Y,X,K,C,n, p.initial, alpha.initial, beta.initial, tol=tol, maxit=maxit)}
    fit.best <- fit
    llik <- fit$loglikelihood
  }
  list(fit = fit.best, llik=llik, alpha.initial=alpha.initial, beta.initial=beta.initial, p.initial=p.initial)
}
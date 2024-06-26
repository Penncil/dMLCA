## Main functions:
## DistLCR: code for simulation; dMLCA with patient-level covariates; equal sample size on sites
## DistLCA: code for MIS-C subphenotyping; dMLCA without patient-level covariates; different sample size among sites

## Compile C functions:
## STEP 1: Delete the current 'Func.All.o' and 'Func.All.so' file, if any.
## STEP 2: In terminal, under the path of the file 'Func.All.c', run 'R CMD SHLIB Func.All.c' (no quote).
## STEP 3: In R, run
dyn.load('Func.All.so')

library(MASS) 
library(gtools)

class_prior <- function(X, beta, alpha, K, n){
  temp1 <- X %*% beta + alpha[rep(seq_len(K), each = n), ]
  temp1 <- pmin(temp1, 7)
  temp1 <- pmax(temp1, -7)
  temp1 <- exp(temp1)
  temp1/as.vector(temp1 %*% rep(1,C))
}

pdf.C <- function(p.last, y, C, K){
  dim.p <- dim(p.last)
  dim.y <- dim(y)
  q <- dim.y[2]
  nK <- dim(y)[1]
  numChoices <- apply(y,2,max)
  lik_vec <- .C("ylik_Ksites",
                as.double(t(p.last)),
                as.integer(t(y)), 
                as.integer(K), 
                as.integer(nK/K), 
                as.integer(q),
                as.integer(numChoices), 
                as.integer(C), 
                lik = double(nK*C))
  matrix(lik_vec$lik, ncol=C, byrow=TRUE)
}
## calculate w
postClass.C <- function(lambda, p.last, y, C, K) {
  dim.p <- dim(p.last)
  dim.y <- dim(y)
  q <- dim.y[2]
  nK <- dim.y[1]
  numChoices <- apply(y,2,max)
  w_vec <- .C("postclass_Ksites",
              as.double(t(lambda)),
              as.double(t(p.last)),
              as.integer(t(y)),
              as.integer(K),
              as.integer(q),
              as.integer(nK/K),
              as.integer(numChoices),
              as.integer(C),
              posterior = double(nK*C)
  )
  matrix(w_vec$posterior, ncol=C, byrow=TRUE)
}
## update p
probHat.C <- function(w, y, C, K) {
  dim.y <- dim(y)
  q <- dim.y[2]
  nK <- dim.y[1]
  numChoices <- apply(y,2,max)
  p_vec <- .C("probhat_Ksites",
              as.integer(t(y)),
              as.integer(K),
              as.double(t(w)),
              as.integer(q),
              as.integer(nK/K),
              as.integer(numChoices),
              as.integer(C),
              ph = double(sum(numChoices)*C)
  )
  p_vec$ph
}
## 1,2-order derivatives of beta
## x is the covariate matrix, no intercept.
dLL2dBeta.C <- function(w,lambda, x, K) {
  ## x must be a matrix.
  C <- dim(lambda)[2]
  numx <- dim(x)[2]
  dbeta_vec <-  .C("d2lldbeta2_Ksites",
                   as.integer(K),
                   as.double(t(w)),
                   as.double(t(lambda)),
                   as.double(t(x)),
                   as.integer(dim(x)[1]/K),
                   as.integer(C),
                   as.integer(numx),
                   grad = double((C-1)*numx),
                   hess = double(((C-1)*numx)^2)                
  )
  list(grad=dbeta_vec$grad,hess=matrix(dbeta_vec$hess,ncol=((C-1)*numx),byrow=TRUE))
}
## 1,2-order derivatives of alpha
d2lldalpha2.C <- function(w, lambda, K){
  C <- dim(lambda)[2]
  dalpha_vec <- .C("d2lldalpha2_Ksites",
                   as.integer(K),
                   as.double(t(w)),
                   as.double(t(lambda)),
                   as.integer(dim(w)[1]/K),
                   as.integer(C), 
                   grad = double((C-1)*K),
                   hess = double(((C-1)*K)^2))
  list(grad=dalpha_vec$grad, hess=matrix(dalpha_vec$hess, ncol=(C-1)*K, byrow=TRUE))
}
## 2-order derivatives of alpha and beta
d2lldab2.C <- function(w, lambda, x, K){
  C <- dim(lambda)[2]
  numx <- dim(x)[2]
  d2ab_vec <- .C("d2lldab2_Ksites",
                 as.integer(K),
                 as.double(t(w)),
                 as.double(t(lambda)),
                 as.double(t(x)),
                 as.integer(dim(x)[1]/K),
                 as.integer(C),
                 as.integer(numx),
                 hess = double( ((C-1)*numx)*(C-1)*K))
  hess <- matrix(d2ab_vec$hess, ncol=(C-1)*K, byrow=TRUE)
  hess
}

lambdaHat.C <- function(w, C, K) {
  dim.w <- dim(w)
  nK <- dim.w[1]
  l_vec <- .C("lambdahat_Ksites",
              as.integer(K),
              as.double(t(w)),
              as.integer(nK/K),
              as.integer(C),
              lh = double(K*C)
  )
  l_vec$lh
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
   p.last_vec <- as.vector(t(p.last))
  p.last.0.1 <- cbind(1-p.last_vec,p.last_vec)
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

########################################################
## categorical data, different sample site across sites
pdf_cat.C <- function(vp, y, C, K, n_k){
  
  dim.y <- dim(y)
  q <- dim.y[2]
  
  numChoices <- apply(y,2,max)
  lik_vec <- .C("ylik_Ksites_cat",
                as.double(vp$vecprobs),
                as.integer(t(y)), 
                as.integer(K), 
                as.integer(n_k), 
                as.integer(q),
                as.integer(vp$numChoices), 
                as.integer(C), 
                lik = double(sum(n_k)*C))
  matrix(lik_vec$lik, ncol=C, byrow=TRUE)
}

postClass_cat.C <- function(lambda, vp, y, C, K, n_k) {
  
  dim.y <- dim(y)
  q <- dim.y[2]
  
  numChoices <- apply(y,2,max)
  w_vec <- .C("postclass_Ksites_cat",
              as.double(t(lambda)),
              as.double(vp$vecprobs),
              as.integer(t(y)),
              as.integer(K),
              as.integer(q),
              as.integer(n_k),
              as.integer(vp$numChoices),
              as.integer(C),
              posterior = double(sum(n_k)*C)
  )
  matrix(w_vec$posterior, ncol=C, byrow=TRUE)
}

probHat_cat.C <- function(w, y, C, K,n_k) {
  dim.y <- dim(y)
  q <- dim.y[2]
  nK <- dim.y[1]
  numChoices <- apply(y,2,max)
  p_vec <- .C("probhat_Ksites_cat",
              as.integer(t(y)),
              as.integer(K),
              as.double(t(w)),
              as.integer(q),
              as.integer(n_k),
              as.integer(numChoices),
              as.integer(C),
              ph = double(sum(numChoices)*C)
  )
  p_vec$ph
}

lambdaHat_cat.C <- function(w, C, K, n_k) {
  
  l_vec <- .C("lambdahat_Ksites_cat",
              as.integer(K),
              as.double(t(w)),
              as.integer(n_k),
              as.integer(C),
              lh = double(K*C)
  )
  l_vec$lh
}

vectorize <- function(probs) {
  classes <- nrow(probs[[1]])
  vecprobs <- unlist(lapply(probs,t))
  numChoices <- sapply(probs,ncol)
  return(list(vecprobs=vecprobs,numChoices=numChoices,classes=classes))
}

unvectorize <- function(vp) {
  probs <- list()
  idx <- c(0,cumsum(vp$numChoices*vp$classes))
  for (i in 1:length(vp$numChoices)){
    probs[[i]] <- matrix(vp$vecprobs[(idx[i]+1):idx[i+1]],nrow=vp$classes,byrow=TRUE)
  }
  return(probs)
}

DistLCA <- function(formula, data, K, C=2, n_k, nrep=1, maxit=1000, tol=1e-3, verbose=TRUE){
  
  mframe <- model.frame(formula, data, na.action = NULL)
 
  Y <- model.response(mframe)
  q <- dim(Y)[2]
  numChoices <- apply(Y,2,max)
  loglike0 <- -Inf
  
  for(rep in 1:nrep){
    ## generate initial values 
    probs <- list()
    ## Generate random initial p
    for (j in 1:q) {
      probs[[j]] <- matrix(runif(C * numChoices[j]),
                           nrow = C, ncol = numChoices[j])
      probs[[j]] <- probs[[j]]/rowSums(probs[[j]])
    }
    probs.init <- probs
    
    p.initial <- vectorize(probs)$vecprobs
    ## Generate random initial lambda
    lambda.initial <- t(rmultinom(K, 100, rep(1/C,C)))/100
    
    p.last <- p.initial
    lambda.last <- lambda.initial
    vp <- list()
    vp$vecprobs <- p.last
    vp$numChoices <- numChoices
    vp$classes <- C
    conv <- FALSE
    iter1 <- 0
    diff_all <- c()
    for(iter in 1:maxit){
      iter1 <- iter1 + 1
      
      lambda_y <- lambda.last[rep(1:K, n_k),]
      w <- postClass_cat.C(lambda_y, vp, Y, C, K, n_k)
      ## update p
      p.current <- probHat_cat.C(w, Y, C, K, n_k)
      
      ## update lambda
      lambda.current <- matrix(lambdaHat_cat.C(w, C, K, n_k), ncol=C, byrow=TRUE)
      
      diff <- sqrt(sum((p.last-p.current)^2)+sum((lambda.last-lambda.current)^2))
      diff_all <- c(diff_all, diff)
      
      lambda.last <- lambda.current
      p.last <- p.current
      vp$vecprobs <- p.last
      if(diff <= tol) {conv = TRUE; break}
    }
    ## Calculate the optimal loglikelihood
    pdf.y <- pdf_cat.C(vp, Y, C, K, n_k)
    ## loglikelihood(Y)
    lambda_y <- lambda.last[rep(1:K, n_k),]
    loglike_ki <- log(as.vector((pdf.y * lambda_y) %*% rep(1,C)))
    loglike <- sum(loglike_ki)
    if (nrep > 1 & verbose) {
      cat("Model ", rep, ": llik = ", loglike,
          " ... best llik = ", max(loglike,loglike0), "\n", sep = "")
      flush.console()
    }
    if(loglike > loglike0){
      p.opt <- p.last
      lambda.opt <- lambda.last
      conv.opt <- conv
      iter1.opt <- iter1
      diff.opt <- diff_all
      loglike0 <- loglike}
  }
  
  
  vp$vecprobs <- p.opt
  phat<- unvectorize(vp)
  names(phat) <- colnames(Y)
  
  for (j in 1:q) {
    rownames(phat[[j]]) <- paste("class ", 1:C,
                                 ": ", sep = "")
    if (is.factor(data[, match(colnames(Y), colnames(data))[j]])) {
      lev <- levels(data[, match(colnames(Y), colnames(data))[j]])
      colnames(phat[[j]]) <- lev
    }
    else {
      colnames(phat[[j]]) <- paste("Pr(", 1:ncol(phat[[j]]),
                                   ")", sep = "")
    }}
  w <- postClass_cat.C(lambda.opt[rep(1:K, n_k),], vp, Y, C, K, n_k)
  predclass <- apply(w, 1, function(x) which.max(x))
  
  list(probs = phat, lambdahat = lambda.opt, loglikelihood = loglike0,
       predclass = predclass,posterior=w,
       niter = iter1.opt, conv = conv.opt, diff=diff.opt)
}
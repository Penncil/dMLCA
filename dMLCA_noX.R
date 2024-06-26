# This file contains functions to run dMLCA without patient-level covariates
## The code here are used in MIS-C subphenotyping.
## Sample size on each site can be different.

## Compile C functions:
## STEP 1: Delete the current 'Func.All.o' and 'Func.All.so' file, if any.
## STEP 2: In terminal, under the path of the file 'Func.All.c', run 'R CMD SHLIB Func.All.c' (no quote).
## STEP 3: In R, run
dyn.load('Func.All.so')
## load required R functions
source("functions.R")

## categorical data, different sample site across sites

###### use vp as input is the best way, fix later#####
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

DistLCA_Cat <- function(formula, data, K, C=2, n_k, nrep=1, maxit=1000, tol=1e-3, verbose=TRUE){
  
  mframe <- model.frame(formula, data, na.action = NULL)
  ## model.frame {stats}: return a data.frame with the variables needed to 
  ## use formula and any ... arguments.
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
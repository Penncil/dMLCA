## example to run dMLCA 
source("dMLCA.R")

DataGenerator <- function(K, C, n ,q, alpha, beta, wright.F=c(0.01,0.02,0.04), p_y = NULL){
 
  N <- K*n
  ## X
  X1 <- rbinom(N , 1, 0.5)
  X2 <- rnorm(N, 0, 1)
  X3 <- rnorm(N, 0, 3)
  X <- matrix(cbind(X1, X2, X3), N)
  ## lambda
  lambda <- class_prior(X, beta, alpha, K, n)
  ## Y
  p.ancestral <- runif(q, min=0.1, max=0.9)
  if (is.null(p_y)){
    temp1 <- ((1-wright.F)/wright.F) %*% matrix(p.ancestral,1)
    temp2 <- ((1-wright.F)/wright.F) %*% matrix(1-p.ancestral,1)
    p_y <- matrix(rbeta(q*C, temp1, temp2), q, C, byrow=TRUE)  
  }
  class_c <- apply(lambda, 1, function(x) which(rmultinom(1,1,x)==1)) ## faster ???
  Y <- matrix(rbinom(q*N, 1, p_y[,class_c]), N, q, byrow=TRUE)

  ## Output in the form of a (n*K) by q matrix
  ## columns are in the order: n obs on site 1, n obs on site 2, ..., n obs on site K.
  list(Y=Y, X=X, p_y=p_y,
       class_y=class_c, lambda_y=lambda)
}

## with covariates
## Step 1: generate simulated data 
## Step 2: run dMLCA
## Set true value of the parameters
C <- 3
K <- 5
q <- 5
n <- 500

beta <- matrix(c(0.5, 1, 0, 1, 0, 0, 0.5, -1, 0), 3, C, byrow=TRUE)
lambda_0 <- matrix(c(15, 20, 65, 34, 33, 33, 50, 25, 25, 65, 20, 15, 10, 15, 75)/100, K, C, byrow=TRUE)
alpha <- log(lambda_0/lambda_0[,C])
p0 <- matrix(c(0.1, 0.5,0.9), q, C)

Data <- DataGenerator(K, C, n ,q, alpha, beta, wright.F=c(0.01,0.02,0.04), p_y = p0)
Y <- Data$Y
X <- Data$X
p0 <- Data$p_y
result <- DistLCR(Y, X, K, C, n, maxit=1000, tol=1e-3, initial = "random", rep = 5)

## without covariates
beta <- matrix(0, 3, C, byrow=TRUE)

Data <- DataGenerator(K, C, n ,q, alpha, beta, wright.F=c(0.01,0.02,0.04), p_y = p0)
Y <- data.frame(Data$Y)
for (i in 1:q) {Y[,i]=factor(Y[,i])}

colnames(Y) <- c('A','B','C','D','E')
fmla <- as.formula(paste('cbind(',paste(colnames(Y),collapse = ','),')~1'))

result <- DistLCA(formula=fmla, data=Y, K=K, C=C, n_k=rep(n,K), nrep=5, maxit=10000, tol=1e-6)


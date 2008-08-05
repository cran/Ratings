

# parse the passed seeds
# 1] if a scalar is passed, it is used by Mersennse twister
# 2] if a list of length two is passed, a parallel-friendly stream is
#    created using L'Ecuyer
"form.seeds" <-
  function(seed) {
    if(length(seed)==1) {
      if(is.na(seed)) seed <- 12345
      seed <- as.integer(seed)
      if(seed < 0) {
        cat("Error: Mersenne seed negative.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)                       
      }
      seeds <- list(0, rep(seed,6), 0)
    }
    if(length(seed)==2) {
      if(!is.list(seed)) {
        cat("Error: List must be passed to use L'Ecuyer.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)          
      }
      lec.seed <- seed[[1]]
      lec.substream <- as.integer(seed[[2]])
      if(is.na(lec.seed[1])) lec.seed <- rep(12345, 6)
      if(length(lec.seed) != 6) {
        cat("Error: L'Ecuyer seed not of length six.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)          
      }
      if(!all(lec.seed >= 0))  {
        cat("Error: At least one L'Ecuyer seed negative.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)          
      }
      if( max(lec.seed[1:3]) >= 4294967087){
        cat("Error: At least one of first three L'Ecuyer seeds\n")
        cat("  greater than or equal to 4294967087\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)          
      }
      if( all(lec.seed[1:3]) == 0 ){
        cat("Error: first three L'Ecuyer seeds == 0\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)          
      }
      if( max(lec.seed[4:6]) >= 4294944443){
        cat("Error: At least one of last three L'Ecuyer seeds\n")
        cat("  greater than or equal to 4294944443\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)          
      }         
      if( all(lec.seed[4:6]) == 0 ){
        cat("Error: last three L'Ecuyer seeds == 0\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)          
      }
      if(lec.substream < 1) {
        cat("Error: L'Ecuyer substream number not positive.\n")
        stop("Please respecify and call ", calling.function(), " again.",
             call.=FALSE)               
      }
      seeds <- list(1, lec.seed, lec.substream) 
    }
    if(length(seed)>2) {
      cat("Error: Seed passed as length greater than two.\n")
      stop("Please respecify and call ", calling.function(), " again.",
           call.=FALSE)        
    }
    return(seeds)
  }



## model for ordinal rating data
##
## i,...,I raters
## j,...,J things rated
## k,...,K rating categories
##
## Y is an I by J matrix with elements in {1,2,...,K, NA}
##      NA denotes missing data assumed to be MAR
##
## y^*_{ij} = alpha_i + beta_i * theta_j + epsilon_{ij}
## epsilon_{ij} ~iid N(0,1)
##
##
##          {1   iff  y^*_{ij}  \in (-infty, 0]
##          {2   iff  y^*_{ij}  \in (0, \gamma_2]
## y_{ij} = ...
##          {K   iff  y^*_{ij}  \in (\gamma_{K-1}, \infty)
##
##
## beta.constraint  inequality contraint on beta
##                       NULL:           no constraint
##                       negative value: beta constrained < 0
##                       positive value: beta constrained > 0 
## theta.neg.index  the index (start with 1) of the negative constrained theta
## theta.pos.index  the index (start with 1) of the positive constrained theta
## vinva            the prior precision of alpha[i]
## vinvb            the prior precision of beta[i]
## ma               the prior mean of alpha[i]
## mb               the prior mean of beta[i]
##
ordrating <- function(Y, beta.constraint=NULL,
                      theta.neg.index=NULL, theta.pos.index=NULL,
                      vinva=0.2, vinvb=0.2, ma=0, mb=1,
                      theta.start=NULL, gamma.start=NULL,
                      burnin=1000, mcmc=10000, thin=1,
                      tune=1, verbose=0, seed=NA){

  ## error checking for Y
  if (min(Y, na.rm=TRUE) < 1){
    cat("Error: min(Y) < 1.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)          
  }
  unique.Y <- unique(as.vector(as.integer(Y)))
  unique.Y <- unique.Y[!is.na(unique.Y)]
  if (length(unique.Y) < max(Y, na.rm=TRUE)){
    warning("Warning: Not all categories are used in Y\n")
  }

  ## error checking for vinva and vinvb
  if (vinva <= 0){
    cat("Error: vinva <= 0.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)          
  }
  if (vinvb <= 0){
    cat("Error: vinvb <= 0.\n")
    stop("Please respecify and call ", calling.function(), " again.",
         call.=FALSE)          
  }

  
  ## get key dimensions
  K <- max(Y, na.rm=TRUE)
  I <- nrow(Y)
  J <- ncol(Y)

  ## set up constraint on beta
  if (is.null(beta.constraint)){
    beta.constraint <- 0
  }
  else if(beta.constraint < 0){
    beta.constraint <- -1
  }
  else if(beta.constraint > 0){
    beta.constraint <- 1
  }

  ## set up constraints on theta
  if (is.null(theta.neg.index)){
    theta.neg.index <- -999
  }
  if (is.null(theta.pos.index)){
    theta.pos.index <- -999
  }
  
  
  
  ## get names
  if (is.null(rownames(Y))){
    rownames(Y) <- 1:nrow(Y)
  }
  if(is.null(colnames(Y))){
    colnames(Y) <- 1:ncol(Y)
  }
  rater.names <- rownames(Y)
  rated.names <- colnames(Y)

  
  ## seeds
  seeds <- form.seeds(seed) 
  lecuyer <- seeds[[1]]
  seed.array <- seeds[[2]]
  lecuyer.stream <- seeds[[3]]
  
  
  ## starting values for theta
  if (is.null(theta.start)){
    theta.start <- apply(Y, 2, mean, na.rm=TRUE)
    theta.start <- (theta.start - mean(theta.start)) / sd(theta.start)
    if (beta.constraint != 0){
      theta.start <- theta.start * beta.constraint
    }

    if (theta.pos.index > 0 & theta.neg.index > 0){
      if (theta.start[theta.pos.index] > theta.start[theta.neg.index]){
        theta.start[theta.pos.index] <- 3
        theta.start[theta.neg.index] <- -3
      }
      else{
        theta.start <- theta.start * -1
        theta.start[theta.pos.index] <- 3
        theta.start[theta.neg.index] <- -3
      }
    }
  }
  else{
    if (length(theta.start) != J){
      stop("theta.start not of length J\n")
    }
  }


  
  ## convert Y into a sparse format where the missing values are
  ## not explicitly represented
  Y.sparse <- NULL
  for (i in 1:I){
    for (j in 1:J){
      if (!is.na(Y[i,j])){
        Y.sparse <- rbind(Y.sparse, c(i-1, j-1, Y[i,j]))
      }
    }
  }

  rm(Y)


  
  ## starting values for alpha and beta
  alpha.start <- rep(0, I)
  beta.start <- rep(beta.constraint, I)

  ## starting values for gamma
  if (is.null(gamma.start)){
    gamma.start <- rep(NA, K+1)
    gamma.start[1] <- -300
    gamma.start[K+1] <- 300
    gamma.start[2] <- 0
    if (K > 2){
      gamma.start[3:K] <- seq(from=1, to=K, length.out=(K-2)) 
    }
  }

  ## setup storage matrix for the sample
  npar <- 2*I + J + (K-1)
  store.mat <- matrix(-999, mcmc/thin, npar)

  posterior <- .C("ordratingpost",
                  samdata = as.double(store.mat),
                  samrow = as.integer(nrow(store.mat)),
                  samcol = as.integer(ncol(store.mat)),
                  Ysparse = as.integer(Y.sparse),
                  Ysparserow = as.integer(nrow(Y.sparse)),
                  Ysparsecol = as.integer(ncol(Y.sparse)),
                  I = as.integer(I),
                  J = as.integer(J),
                  K = as.integer(K),
                  betaconstraint = as.integer(beta.constraint),
                  thetanegindex = as.integer(theta.neg.index-1),
                  thetaposindex = as.integer(theta.pos.index-1),
                  vinva = as.double(vinva),
                  vinvb = as.double(vinvb),
                  ma = as.double(ma),
                  mb = as.double(mb),
                  burnin = as.integer(burnin),
                  mcmc = as.integer(mcmc),
                  thin= as.integer(thin),
                  tune = as.double(tune),
                  verbose = as.integer(verbose),
                  alpha = as.double(alpha.start),
                  beta = as.double(beta.start),
                  theta = as.double(theta.start),
                  gamma = as.double(gamma.start),
                  lecuyer=as.integer(lecuyer),
                  seedarray=as.integer(seed.array),
                  lecuyerstream=as.integer(lecuyer.stream),
                  PACKAGE="Ratings"
                  )


  store.mat <- matrix(posterior$samdata, posterior$samrow, posterior$samcol,
                      byrow=FALSE)
  output <- mcmc(data=store.mat, start=1, end=mcmc, thin=thin)
  
  parnames <- c(paste("alpha", rater.names, sep=""),
                paste("beta", rater.names, sep=""),
                paste("theta", rated.names, sep=""),
                paste("gamma", 1:(K-1), sep="")
                )

  varnames(output) <- parnames

  return(output)
  
} ## end ordrating

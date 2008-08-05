## posterior predictive probabilities from ordrating model


## the posterior predictive part
## returns a number of number of products by number of rating categories
## matrix with typical element tau_{pc} 
## out is an output object from ordrating
tauCalculate <- function(out, ndraws=500){
  
  
  theta <- out[,grep("^theta", colnames(out))]
  theta.median <- apply(theta, 2, median)
  alpha <- out[,grep("^alpha", colnames(out))]
  beta <- out[,grep("^beta", colnames(out))]
  gamma <- out[,grep("^gamma", colnames(out))]
    
  ncat <- ncol(gamma) + 1
  npaper <- ncol(theta)
  nrater <- ncol(alpha)
  M <- nrow(out)
  
  starprobmat <- matrix(0, npaper, ncat)

  samp.inds <- sample(1:M, size=ndraws, replace=FALSE)

  count <- 1
  for (i in samp.inds){
    ## linear predictor
    eta <- as.numeric(alpha[i,]) +  outer(beta[i,], theta[i,])

    ## fill in starprob
    starprobmat[,1] <- starprobmat[,1] + apply(pnorm(-eta), 2, sum)
    for (k in 2:(ncat-1)){
      starprobmat[,k] <- starprobmat[,k] +
        apply(pnorm(gamma[i,k] - eta) - pnorm(gamma[i,k-1] - eta), 2, sum)
    }
    starprobmat[,ncat] <- starprobmat[,ncat] +
      apply(1 - pnorm(gamma[i,ncat-1] - eta), 2, sum) 

    if (count %% 100 == 0){
      cat("mcmc iteration", count, "of", ndraws, "\n")
    }

    count <- count + 1
  }

  starprobmat <- starprobmat / (nrater * ndraws)

  rownames(starprobmat) <- colnames(theta)
  colnames(starprobmat) <- paste("ProbY", 1:ncat, sep="")
  
  starprobmat <- starprobmat[order(theta.median),]

  return(starprobmat)
  
}




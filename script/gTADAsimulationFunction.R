
gTADAsimulator <- function(N, muAll, 
                       alpha0, geneSet = NULL, ## Gene set parameters
                       beta0 = NULL, ## Mutation-rate parameter
                       gamma.mean.dn, beta.dn,
                       gamma.mean.CC, beta.CC,
                       rho1, nu1, rho0, nu0, tradeoff=FALSE,
                       adjustMutationRate = FALSE) {
      m <- dim(muAll)[1] # number of genes
      K <- dim(muAll)[2] # number of mutational categories
      exp0 <- rep(alpha0[1], m)
#      if (!is.null(geneSet)){
 #         for (i in 1:dim(geneSet)[2])
  #             exp0 <- exp0 + alpha0[1+i]*geneSet[, i]
   #   if (adjustMutationRate)
    #      exp0 <- exp0 + beta0*log(apply(muAll, 1, sum) + 1, base = 10)
     # }

      pp0 <- exp(exp0)/(1 + exp(exp0))
      z0 <- sapply(pp0, function(x)
          sample(0:1, 1, prob = c(1 - x, x)))
      

      z = z0
        gamma.dn <- array(1, dim=c(m,K))
        gamma.CC <- array(1, dim=c(m,K))
        q <- array(0, dim=c(m,K))
        x <- array(0, dim=c(m,3*K))
        k <- sum(z==1)
      print(k)
      print(gamma.mean.dn)
      print(m) 
        for (j in 1:K) {
                mu <- muAll[, j]
                    # sample de novo
                    gamma.dn[z==1, j] <- rgamma(k, gamma.mean.dn[j]*beta.dn[j], beta.dn[j])
                    col <- 3*(j-1)+1
                    x[,col] <- rpois(m, 2 * mu  * gamma.dn[,j] * N$dn[j])

                    # sample case-control
                    gamma.CC[z==1, j] <- rgamma(k, gamma.mean.CC[j]*beta.CC[j], beta.CC[j])
                    q[z==0, j] <- rgamma(m-k, rho0[j], nu0[j])
                    if (tradeoff==FALSE) {
                              q[z==1, j] <- rgamma(k, rho1[j], nu1[j])
                    } else {
                              q[z==1, j] <- mu[z==1] / (delta[j] * gamma.CC[z==1, j])
                    }
                    x[,col+1] <- rpois(m, q[,j] * gamma.CC[,j] *N$ca[j])
                    x[,col+2] <- rpois(m, q[,j] * N$cn[j])

        }

        sample.info <- cbind(z, gamma.dn, gamma.CC, q, x, muAll)

        return (list(sample=x, sample.info=sample.info))
}

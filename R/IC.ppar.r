IC.ppar <- function(object)
{
#computes unconditional loglik, AIC, BIC, and cAIC
#object of class ppar
#k...penalty parameter

probmat <- pmat(object)                         #matrix with solving probabilities
loglik <- sum(rowSums(log(probmat),na.rm=TRUE))            #unconditional log-likelihood

npar <- dim(object$W)[2]                        #number of item parameters
N <- dim(object$X)[1]                           #number of persons
AIC <- -2*loglik + 2*npar
BIC <- -2*loglik + log(N)*npar
cAIC <- -2*loglik + log(N)*npar + npar

result <- list(j.loglik=loglik,AIC=AIC,BIC=BIC,cAIC=cAIC)
class(result) <- "ICr"
result
}

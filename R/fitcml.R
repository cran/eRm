"fitcml" <-
function(mt_ind,nr,x_mt,rtot,W)
{ 

#cml function for call in nlm
cml <- function(eta)
{
beta <- as.vector(W%*%eta)
#beta <- c(0, beta1)
beta_list <- as.list(split(beta,mt_ind))    #beta as list
parlist <- lapply(beta_list,exp)           #initial epsilon as list

g_iter <- NULL                                                  #computation of the first gamma derivatives
K <- length(parlist)
for (t in 1:(K-1)) {                                                    #building up J1,...,Jt,...,Js

  if (t==1) {                                                                 #first iteration step
    gterm <- c(1,parlist[[t]])                                  #0th element included
  }else
  {
    gterm <- g_iter                                                     #gamma previous iteration with 0th el
    g_iter <- NULL
  }

  parvek <- c(1,parlist[[t+1]])                                 #eps vector in current iteration with 0th el

  h <- length(parvek)                                                   #dimensions for matrix
  mt <- length(gterm)
  #rtot1 <- rtot+1
  rtot1 <- h+mt-1                              #number of possible raw scores (0 included)

  gtermvek <- rep(c(gterm,rep(0,h)),h)  #building up matrix for gamma term
  gtermvek <- gtermvek[-((length(gtermvek)-h+1):length(gtermvek))]      #eliminating last h 0's
  gmat <- matrix(gtermvek,nrow=rtot1,ncol=h)
  emat <- matrix(rep(parvek,rep(rtot1,h)),ncol=h,nrow=rtot1)              #building up matrix for eps term
  gmat_new <- gmat*emat                                               #merge matrices
  g_iter <- rowSums(gmat_new)                                   #gamma functions in current iteration are rowsums
}

Lg <- list(gamma=g_iter[2:(rtot+1)])                                    #final gamma vector stored in gamma (without gamma0)

#----------------- log-likelihood -----------------------
L1 <- nr%*%log(Lg$gamma)    
L2 <- x_mt%*%beta
L1-L2
#----------------- end likelihood -----------------------
}
############ end gamma functions ##########################

eta <- c(rep(0,dim(W)[2]))                                  #starting values for eta parameters
fit <- nlm(cml,eta,hessian=TRUE)        
#fit1 <- optim(eta,cml,method="BFGS",hessian=T)             #optim and nlm statement
#list(fit=fit)
}


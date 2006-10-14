`fitcml` <-
function (mt_ind,nrlist,x_mt,rtot,W,ngroups,gind,x_mtlist)
{

#cml function for call in nlm
cml <- function(eta)
{

beta <- as.vector(W%*%eta)

#beta_list <- as.list(split(beta[1:(length(beta)/ngroups)],mt_ind[1:(length(mt_ind)/ngroups)]))                         #beta as list
#beta_list <- as.list(split(beta[(length(beta)/ngroups+1):length(beta)],mt_ind[(length(beta)/ngroups+1):length(beta)]))
#mt_ind1 <- mt_ind[1:(length(mt_ind)/ngroups)]

#beta_list <- as.list(split(beta,mt_ind))

Lg <- tapply(beta,gind, function(betag) {         #computation of the gamma-functions for each subgroup
                                                                      
      beta_list <- as.list(split(betag,mt_ind[1:(length(betag))]))  #list of virtual item parameters per item
      parlist <- lapply(beta_list,exp)                                #initial epsilon as list

      g_iter <- NULL                                                  #computation of the gamma functions
      K <- length(parlist)

      for (t in 1:(K-1)) {                                            #building up J1,...,Jt,...,Js

        if (t==1) {                                                   #first iteration step
          gterm <- c(1,parlist[[t]])                                  #0th element included
        }else
        {
         gterm <- g_iter                                              #gamma previous iteration with 0th el
         g_iter <- NULL
        }

        parvek <- c(1,parlist[[t+1]])                                 #eps vector in current iteration with 0th el
        h <- length(parvek)                                           #dimensions for matrix
        mt <- length(gterm)
        rtot1 <- h+mt-1                                               #number of possible raw scores (0 included)

        gtermvek <- rep(c(gterm,rep(0,h)),h)                          #building up matrix for gamma term
        gtermvek <- gtermvek[-((length(gtermvek)-h+1):length(gtermvek))]      #eliminating last h 0's
        gmat <- matrix(gtermvek,nrow=rtot1,ncol=h)
        emat <- matrix(rep(parvek,rep(rtot1,h)),ncol=h,nrow=rtot1)            #building up matrix for eps term
        gmat_new <- gmat*emat                                                 #merge matrices
        g_iter <- rowSums(gmat_new)                                           #gamma functions in current iteration are rowsums
      }

     list(g_iter[2:(rtot+1)])                                                 #final gamma vector stored in gamma (without gamma0)
     })
#----------------- end gamma functions ------------------

#----------------- log-likelihood -----------------------


L1 <- sum(mapply("%*%",nrlist,lapply(Lg,log)))                                            #summing up L1-terms (group-wise)
#L1 <- nr%*%log(Lg$gamma)      #nr aufsplitten
beta.glist <- split(beta,gind)
L2 <- sum(mapply("%*%",x_mtlist,beta.glist))                                              #summing up L2-terms (group-wise)
#L2 <- x_mt%*%beta
L1-L2                                                                                     #actual likelihood value
#----------------- end likelihood -----------------------
}


eta <- rep(0,dim(W)[2])                                   #starting values for eta parameters
#eta <- jitter(eta)

fit <- nlm(cml,eta,hessian=TRUE)
#fit <- optim(eta,cml,method="BFGS",hessian=TRUE) 
}


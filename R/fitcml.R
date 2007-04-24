`fitcml` <-
function (mt_ind,nrlist,x_mt,rtot,W,ngroups,gind,x_mtlist,NAstruc,g_NA,st.err)
{

#cml function for call in nlm
cml <- function(eta)
{

beta <- as.vector(W%*%eta)

beta.list <- split(beta,gind)      

#---------------- NEW -----------------------
#if (length(g_NA) > 1) {
#  beta.vec1 <- unlist(mapply(function(bl,gn) {
#                         rep(bl,gn)
#                       },beta.list,g_NA,SIMPLIFY=TRUE))
#  bl.ind <- rep(1:sum(g_NA), each = length(beta.list[[1]]))
#  beta.list1 <- split(beta.vec1,bl.ind)
#} else {
beta.list1 <- beta.list
#}
#---------------- end NEW -------------------

betaNA <- mapply(function(x,y) {rbind(x,y)},beta.list1,NAstruc,SIMPLIFY=FALSE)         #beta and NAstructure as list (over Groups)


Lg <- lapply(betaNA, function(betaNAmat) {        #gamma functions for each Group x NAgroup combination 

         beta.vec <- betaNAmat[1,]                #get parameter vector beta
         
         Lg.NA <- apply(matrix(betaNAmat[-1,],ncol=length(beta.vec)),1, function(NAvec) {                 #likelihood for each NAgroup within Groups                                          
            
            beta_list <- as.list(split(beta.vec[NAvec==1],mt_ind[1:(length(beta.vec[NAvec==1]))]))        #list of virtual item-category parameters per item
                      
            parlist <- lapply(beta_list,exp)                                #initial epsilon as list
      
            #gamma functions
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
      
           Lg.NA <- as.vector(g_iter[2:(rtot+1)])                                                 #final gamma vector stored in gamma (without gamma0)
           return(Lg.NA)
           }) })          
#----------------- end gamma functions ------------------

#----------------- log-likelihood -----------------------
                               
#=========to be deleted
#L1t <- (mapply(function(x,z) {
#                   x[!is.na(z)]%*%na.exclude(z)
#                   },nrlist,lapply(Lg,log)))          #sum up L1-terms (group-wise)
#L2t <- (mapply("%*%",x_mtlist,beta.list1))            #sum up L2-terms (group-wise)
#print(L1t-L2t)
#==========end delete


L1 <- sum(mapply(function(x,z) {
                   x[!is.na(z)]%*%na.exclude(z)
                   },nrlist,lapply(Lg,log)))          #sum up L1-terms (group-wise)

L2 <- sum(mapply("%*%",x_mtlist,beta.list1))           #sum up L2-terms (group-wise)

L1-L2                                               #actual likelihood value
#print(L1-L2)                                              
#----------------- end likelihood -----------------------
}


eta <- rep(0,dim(W)[2])                                   #starting values for eta parameters

options(warn=-1)                                          #turn off warnings for NA/Inf
fit <- nlm(cml,eta,hessian=st.err)                        #NLM optimizer

#options(warn=0)

#fit <- optim(eta,cml,method="BFGS",hessian=TRUE) 
}


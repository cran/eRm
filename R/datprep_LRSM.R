"datprep_LRSM" <-
function(X,W,mpoints,Groups)
{
  #TFrow <- (rowSums(X)==0)                       #el. persons with 0 rawscore
  #X <- X[!TFrow,]

  ngroups <- max(Groups)                          #number of groups
  N <- dim(X)[1]                                  #number of persons
  K <- dim(X)[2]/mpoints                          #number of items
  hmax <- max(X)                                  #highest category
  mt_vek <- rep(hmax,K)                           #number of categories - 1 for each item                  
  mt_vek_0 <- mt_vek+1                            #number of categories for each item
  
  X01_0 <- matrix(rep(0,(N*sum(mt_vek_0)*mpoints)),nrow=N) #empty 0/1 matrix  
  K1 <- dim(X)[2]
  cummt0 <- c(0,cumsum(rep(mt_vek_0,mpoints))[1:(K1-1)])+1     #index vector for 0th category
  indmatp <- apply(X,1,function(xi) {xi+cummt0})  #preparing index matrix for 1 responses
  imp1 <- as.vector(indmatp)
  imp2 <- rep(1:N,rep(K1,N))
  indmat <- cbind(imp2,imp1)                      #final index matrix for 1 responses
  X01_0[indmat] <- 1                              #0/1 matrix with 0th category
  X01 <- X01_0[,-cummt0]
  
  #automatized generation of the design matrix W
  if (length(W)==1) {                             #generating design matrix
    e_it <- gl(K,hmax)                            #factor for item parameters
    c_cat <- gl(hmax,1,K*hmax)                    #factor for category par
    Xm <- model.matrix(~e_it+c_cat)[,-1]          #design matrix with 0/1 contrasts (without intercept)
    catvek <- 1:hmax                              #preparing the item design vectors
    e_itnew <- catvek*Xm[,1:(K-1)]                  
    Xm[,1:(K-1)] <- e_itnew
    W11 <- Xm                                     #first part (same as RSM) without virtual items
    ZW <- dim(W11)[1]
    
    W1 <- NULL
    for (i in 1:(mpoints*ngroups)) W1 <- rbind(W1,W11)  #first part with virtual items
    
    if (mpoints > 1) {                            #more than 1 measurement points
      if (ngroups > 1) {                          #more than 1 group/more mpoints
        t_mp1 <- rep(1:mpoints,rep(ZW*ngroups,mpoints))
        t_mp <- factor(t_mp1)
        g_ng1 <- rep(rep(1:ngroups,rep(ZW,ngroups)),mpoints)
        g_ng <- factor(g_ng1)
        W2 <- model.matrix(~t_mp+g_ng+t_mp*g_ng)[,-1]     #full design (main effects g and mp, interactions)
      } else {                                    #1 group/more mpoints
        t_mp <- gl(mpoints,ZW)             #factor for measurement points
        W2 <- model.matrix(~t_mp)[,-1] }
    } else if (ngroups > 1) {                     #1 mpoint/more groups
        g_ng <- gl(ngroups,ZW)
        W2 <- model.matrix(~g_ng)[,-1] 
    } else if (ngroups == 1) W2 <- NULL           #1 mpoint/1 group
        
  W2_cat <- W2*catvek                             #imposing item categories
  W <- cbind(W1,W2_cat)                           #design matrix completed
  }
  
  list(X=X,X01=X01,mt_vek=mt_vek,W=W)
}


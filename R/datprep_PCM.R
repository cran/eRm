`datprep_PCM` <-
function(X,W)
{
#... X: data matrix with response categories to be converted into 0/1 matrix

  #TFrow <- (rowSums(X)==0)  #el. persons with 0/K rawscore
  #X <- X[!TFrow,]

  #converting into 0/1 matrix
  N <- dim(X)[1]                                  #number of persons
  mt_vek <- apply(X,2,max)                        #number of categories - 1 for each item
  mt_vek_0 <- mt_vek+1                            #number of categories for each item
  X01_0 <- matrix(rep(0,(N*sum(mt_vek_0))),nrow=N) #empty 0/1 matrix
  K <- length(mt_vek)                             #number of items
  cummt0 <- c(0,cumsum(mt_vek_0)[1:(K-1)])+1      #index vector for 0th category
  indmatp <- apply(X,1,function(xi) {xi+cummt0})  #preparing index matrix for 1 responses
  imp1 <- as.vector(indmatp)
  imp2 <- rep(1:N,rep(K,N))
  indmat <- cbind(imp2,imp1)                      #final index matrix for 1 responses
  X01_0[indmat] <- 1                              #0/1 matrix with 0th category
  X01 <- X01_0[,-cummt0]                          #delete 0-category answers --> final 0/1 pattern matrix (dim N*sum(mt_vek))

  #automatized generation of the design matrix W
  if (length(W)==1) {
    W1 <- diag(1,(sum(mt_vek)-1))                   #build up design matrix
    w1 <- rep(0,(sum(mt_vek)-1))                    #first item category parameter set to 0
    W <- rbind(w1,W1)                               #PCM design matrix 
  }

  list(X=X,X01=X01,mt_vek=mt_vek,W=W)
#Output: X01      ... 0/1 response matrix of dimension N*rtot
#        mt_vek   ... vector of length K with number of categories - 1 (for each item)
#        W        ... design matrix of dimension sum(mt_vek)*sum(mt_vek)
}


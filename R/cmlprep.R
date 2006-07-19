"cmlprep" <-
function(X01,mt_vek,mpoints,Groups)
{
  
  if (length(Groups)==1) { 
     x_mt <- colSums(X01)                     #item category raw scores as vector
     ngroups <- 1
  } else                                      #if groups defined
  {
    Xord <- X01[order(Groups),]               #sorting by group specification
    Xlistvek <- split(Xord,sort(Groups))      #list of vectors
    Xlistmat <- lapply(Xlistvek,function(x) {x <- matrix(x,ncol=dim(X01)[2],byrow=TRUE)})       #list of splitted data matrices
    x_mtmat <- sapply(Xlistmat,colSums)
    x_mt <- as.vector(x_mtmat)                #item (category) raw scores for each group
    ngroups <- max(Groups)                    #number of groups
    
    #old code with Groups as matrix
    #x_mtmat <- apply(Groups,1,function(g) {
    #                            x_mtrow <- colSums(X01[g[1]:g[2],])})
    #x_mt <- as.vector(x_mtmat)
    #ngroups <- dim(Groups)[1]
  }
  
  end1 <- length(mt_vek)*mpoints*ngroups
  mt_ind <- rep(1:end1,rep(mt_vek,mpoints*ngroups))     #category index vector (for converting x_mt into list)
  x_tmt <- split(x_mt,mt_ind)                       #list for likelihood
  rtot <- sum(mt_vek)*mpoints

  ics <-  rep(sequence(mt_vek),mpoints)             #item category scores for each item as vector
  rv <- apply(X01,1,function(x) ics%*%x)            #person raw scores
  rvtab <- table(rv)                                #raw score frequencies
  dnamevek <- as.numeric(unlist(dimnames(rvtab)))   #different raw scores for 0 fill up
  nr <- rep (0,rtot+1)                              #setting 0 raw score frequencies
  nr[dnamevek+1] <- rvtab                           #vector with person raw scores from 1:rtot (with 0 fill up)
  nr <- nr[-1]

list(x_mt=x_mt,mt_ind=mt_ind,x_tmt=x_tmt,rtot=rtot,nr=nr)
}


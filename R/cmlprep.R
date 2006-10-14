`cmlprep` <-
function(X01,mt_vek,mpoints,Groups,W)
{

  if (length(Groups)==1) {
     x_mt <- colSums(X01)                     #item category raw scores as vector
     x_mtlist <- list(x_mt) 
     ngroups <- 1
  } else {                                     #if groups defined
    ngroups <- max(Groups)                    #number of groups
    x_mtlist <- by(X01,Groups,colSums)        #item-category person raw scores for each group (as list)
    x_mt <- as.vector(unlist(x_mtlist))       #as vector: g1|g2|...
  }

  end1 <- length(mt_vek)*mpoints*ngroups
  mt_ind <- rep(1:end1,rep(mt_vek,mpoints*ngroups)) #category index vector (for converting x_mt into list)
  x_tmt <- split(x_mt,mt_ind)                       #list for likelihood: item-wise * ngroups
  rtot <- sum(mt_vek)*mpoints

  ics <-  rep(sequence(mt_vek),mpoints)             #item category scores for each item as vector
  rv <- apply(X01,1,function(x) ics%*%x)            #person raw scores

  if (ngroups > 1) {                                      #groups
    w.last <-  W[,dim(W)[2]]
    seglen <- sum(mt_vek)                               #length of beta vektor (called segment)
    gind <- rep(rep(1:ngroups,rep(seglen,ngroups)),mpoints)                 #index vector for group extraction
  } else {
    gind <- rep(1,dim(W)[1])
  }
  
  rvlist <- split(rv,Groups)
  nrlist <- lapply(rvlist,function(rvel) {                                    #matrix with raw score frequencies for each group (transposed)
                            rvtab <- table(rvel)                              #raw score frequencies
                            dnamevek <- as.numeric(unlist(dimnames(rvtab)))   #different raw scores for 0 fill up
                            nr <- rep (0,rtot+1)                              #setting 0 raw score frequencies
                            nr[dnamevek+1] <- rvtab                           #vector with person raw scores from 1:rtot (with 0 fill up)
                            nr <- nr[-1]
                            return(nr)
                          })
list(x_mt=x_mt,mt_ind=mt_ind,x_tmt=x_tmt,rtot=rtot,nrlist=nrlist,gind=gind,x_mtlist=x_mtlist)
}


pifit.internal <- function(object)
{
#object of class ppar
#function is called in itemfit.ppar and personfit.ppar


  X <- object$X
  mt_vek <- apply(X,2,max,na.rm=TRUE)             #number of categories - 1 for each item
  mt_ind <- rep(1:length(mt_vek),mt_vek)
  mt_seq <- sequence(mt_vek)
  gmemb <- object$gmemb

  pmat <- pmat(object)                          #matrix with model probabilites

  #-----------------matrix with expected response patterns--------------
  Emat.cat <- t(apply(pmat,1,function(x) x*mt_seq))
  if ((object$model == "RM") || (object$model == "LLTM")) { 
    Emat <- Emat.cat
  } else {
    E.list <- tapply(1:length(mt_ind),mt_ind, function(ind) {rowSums(Emat.cat[,ind],na.rm=TRUE)})
    Emat <- matrix(unlist(E.list),ncol=dim(X)[2],dimnames=list(rownames(pmat),colnames(X)))
  } 
  #------------------------variance term for standardized residuals------
  pmat.l0 <- tapply(1:length(mt_ind),mt_ind, function(ind) {
                            vec0 <- 1-rowSums(as.matrix(pmat[,ind]),na.rm=TRUE)     #prob for 0th category
                            cbind(vec0,pmat[,ind])
                            })
  pmat0 <- matrix(unlist(pmat.l0),nrow=length(gmemb))   #prob matrix 0th category included
  mt_vek0 <- integer(0)                                 #add 0th category to all indices
  for (i in mt_vek) mt_vek0 <- c(mt_vek0, 0:i)
  mt_ind0 <- rep(1:length(mt_vek),mt_vek+1)
  colnames(Emat) <- NULL
  Emat0 <- t(apply(Emat[,mt_ind0],1,function(x) {mt_vek0 - x}))
  Vmat.cat <- (Emat0)^2*pmat0
  V.list <- tapply(1:length(mt_ind0),mt_ind0, function(ind) {rowSums(Vmat.cat[,ind],na.rm=TRUE)})
  Vmat <- matrix(unlist(V.list),ncol=dim(X)[2],dimnames=list(rownames(pmat),colnames(X)))

  result <- list(Emat=Emat,Vmat=Vmat)

}
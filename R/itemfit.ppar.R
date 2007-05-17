`itemfit.ppar` <-
function(object)
# computes Chi-square based itemfit statistics (Smith, p.77ff)
# for object of class "ppar" (from person.parameter)
# Note: They have to be computed for each NA group separately since for the 
# pmat-computation the person parameter are needed (which are different for the subgroups) 
{
  X <- object$X
  mt_vek <- apply(X,2,max,na.rm=TRUE)             #number of categories - 1 for each item
  mt_ind <- rep(1:length(mt_vek),mt_vek)
  
  if (any(is.na(X))) {
    dichX <- ifelse(is.na(X),1,0)
    strdata <- apply(dichX,1,function(x) {paste(x,collapse="")})
    gmemb <- as.vector(data.matrix(data.frame(strdata)))
  } else {
    gmemb <- rep(1,dim(X)[1])
  }

  rp <- rowSums(X,na.rm=TRUE)
  maxrp <- sum(mt_vek)
  TFrow <- ((rp==maxrp) | (rp==0))              #don't regard persons with 0/full raw score
  X01 <- object$X01[!TFrow,] 
  gmemb <- gmemb[!TFrow]
  
  X01.l <- by(X01,gmemb,function(x) x)   #create list of X01 due to NAstructure
  Pmat.l <- pmat(object)                        #list of prob matrices due to Nastructure
  
  st.res <- mapply(function(x,p) {              #standardized residuals
                   res <- (x-p)/sqrt(p*(1-p))   #residual matrix
                   return(res)
                  }, X01.l,Pmat.l,SIMPLIFY=FALSE)

  sq.res <- lapply(st.res,function(x) {x^2})        #squared standardized residuals
  
  #item outfit
  ifit <- lapply(sq.res,colSums)
  
  i.df <- lapply(X01.l,function(x) {dim(x)[1]})     #degrees of freedom (eventually -1 since one parameter is fixed??)
  #i.df <- lapply(X01.l,function(x) {1}) 

  result <- list(i.fit=ifit,i.df=i.df)
  class(result) <- "ifit"
  result
}


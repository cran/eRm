`personfit.ppar` <-
function(object)
# computes Chi-square based itemfit statistics (Smith, p.77ff)
# for object of class "ppar" (from person.parameter)
{

  X <- object$X
  if (any(is.na(X))) {
    dichX <- ifelse(is.na(X),1,0)
    strdata <- apply(dichX,1,function(x) {paste(x,collapse="")})
    gmemb <- as.vector(data.matrix(data.frame(strdata)))
  } else {
    gmemb <- rep(1,dim(X)[1])
  }

  X01.l <- by(object$X01,gmemb,function(x) x)   #create list of X01 due to NAstructure
  Pmat.l <- pmat(object)                        #list of prob matrices due to Nastructure

  st.res <- mapply(function(x,p) {              #standardized residuals
                   res <- (x-p)/sqrt(p*(1-p))   #residual matrix
                   return(res)
                  }, X01.l,Pmat.l,SIMPLIFY=FALSE)

  sq.res <- lapply(st.res,function(x) {x^2})        #squared standardized residuals

  #personfit
  p.fit <- lapply(sq.res,rowSums,na.rm=TRUE)

  #person infit
  
  sum.p.df <- lapply(sq.res,colSums)
  p.df <- lapply(sum.p.df,function(y) length(na.omit(y)))

  result <- list(p.fit=p.fit,p.df=p.df)
  class(result) <- "pfit"
  result
}


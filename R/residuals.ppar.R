`residuals.ppar` <-
function(object,...)
# computes standardized and squared residuals
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
  
  st.res <- lapply(st.res,function(st) {
                    colnames(st)<- paste("Item",(1:dim(st)[2]))
                    return(st)})

  sq.res <- lapply(sq.res,function(sq) {
                    colnames(sq)<- paste("Item",(1:dim(sq)[2]))
                    return(sq)})
                    
  result <- list(st.res=st.res,sq.res=sq.res)
  class(result) <- "resid"
  result
}


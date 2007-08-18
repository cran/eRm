`personfit.ppar` <-
function(object)
# computes Chi-square based itemfit statistics (Smith, p.77ff)
# for object of class "ppar" (from person.parameter)
{

  if (length(object$pers.ex)==0) {
    X <- object$X
  } else { 
    X <- object$X[-object$pers.ex,]
  }
  
  #rp <- rowSums(X,na.rm=TRUE)
  #mt_vek <- apply(X,2,max,na.rm=TRUE)
  #maxrp <- sum(mt_vek)
  #TFrow <- ((rp==maxrp) | (rp==0))              #exclude full and 0 responses
     
  VE <- pifit.internal(object)                  #compute expectation and variance term
  Emat <- VE$Emat
  Vmat <- VE$Vmat
  
  st.res <- (X-Emat)/sqrt(Vmat)
  #st.res <- (X[!TFrow,]-Emat)/sqrt(Vmat)
  
  sq.res <- st.res^2                            #squared standardized residuals
  pfit <- rowSums(sq.res,na.rm=TRUE)            
  
  pdf <- apply(X,1,function(x) {length(na.exclude(x))})
  
  #pdf <- apply(X[!TFrow,],1,function(x) {length(na.exclude(x))})   #degress of freedom (#of persons per item)
  
  result <- list(p.fit=pfit,p.df=pdf,st.res=st.res)
  class(result) <- "pfit"
  result
}


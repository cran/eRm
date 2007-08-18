`itemfit.ppar` <-
function(object)
# computes Chi-square based itemfit statistics 
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
  ifit <- colSums(sq.res,na.rm=TRUE)            
                                      
  idf <- apply(X,2,function(x) {length(na.exclude(x))})
  #idf <- apply(X[!TFrow,],2,function(x) {length(na.exclude(x))})   #degress of freedom (#of persons per item)
  
  result <- list(i.fit=ifit,i.df=idf,st.res=st.res)
  class(result) <- "ifit"
  result
}


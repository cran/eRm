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
  Cmat <- VE$Cmat

  st.res <- (X-Emat)/sqrt(Vmat)
  #st.res <- (X[!TFrow,]-Emat)/sqrt(Vmat)

  sq.res <- st.res^2                            #squared standardized residuals
  pfit <- rowSums(sq.res,na.rm=TRUE)

  pdf <- apply(X,1,function(x) {length(na.exclude(x))})

  #pdf <- apply(X[!TFrow,],1,function(x) {length(na.exclude(x))})   #degress of freedom (#of persons per item)

  p.outfitMSQ <- pfit/pdf

  qsq.outfitMSQ <- (rowSums(Cmat/Vmat^2, na.rm=TRUE)/pdf^2) - 1/pdf
  q.outfitMSQ <- sqrt(qsq.outfitMSQ)

  psumVmat<-rowSums(Vmat)
  p.infitMSQ <- rowSums(sq.res*Vmat, na.rm = TRUE)/psumVmat

  qsq.infitMSQ <- rowSums(Cmat-Vmat^2, na.rm=TRUE)/psumVmat^2
  q.infitMSQ <- sqrt(qsq.infitMSQ)

  p.outfitZ <- ((p.outfitMSQ)^(1/3)-1)*(3/q.outfitMSQ)+(q.outfitMSQ/3)
  p.infitZ <- ((p.infitMSQ)^(1/3)-1)*(3/q.infitMSQ)+(q.infitMSQ/3)

  result <- list(p.fit = pfit, p.df = pdf, st.res = st.res, p.outfitMSQ = p.outfitMSQ,
                 p.infitMSQ = p.infitMSQ,
                 p.outfitZ = p.outfitZ, p.infitZ = p.infitZ)
  class(result) <- "pfit"
  result
}


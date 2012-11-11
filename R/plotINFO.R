plotINFO <- function(ermobject,type="both",theta=seq(-5,5,0.01),...)
  {
    if(type=="both") par(mfrow=c(2,1))
    if(type=="item"|type=="both")
      {
        iinfo <- item_info(ermobject,theta)
        info <- lapply(iinfo,function(x)x$i.info)
        pltinfo <-matrix(unlist(info),ncol=dim(ermobject$X)[2]) 
        matplot(theta,pltinfo,type="l",main="Item Information",ylab="Information",xlab="Latent Trait")
        itmnames <- paste("Item",1:dim(ermobject$X)[2])
        legend("topright",itmnames, pch = NULL, lty= 1:5, col=1:6)
      }
    if(type=="test"|type=="both")
      {
        tinfo <- test_info(ermobject,theta)
        plot(theta,tinfo,type="l",main="Test Information",xlab="Latent Trait",ylab="Scale Information")
      }
    par(mfrow=c(1,1))
  }

`plotICC.Rm` <-
function(object, item.subset = "all", empirical = FALSE, xlim = c(-4,4), ylim = c(0,1), 
         xlab = "Latent Dimension", ylab = "Probability to Solve",...)
# produces ICC plots
# object of class Rm
{
  
  X <- object$X  

  if (empirical) {                                       #empirical ICC for Rasch model only
    th.est <- person.parameter(object,se=FALSE,splineInt=TRUE)
    thetapar <- th.est$thetapar
    if (length(thetapar)==1) {                           #Too complicated with NA'groups (for each NAgroup separate plots...)
      emp.plot <- TRUE
      thetapar.u <- unique(round(unlist(thetapar),5))
    } else {
      emp.plot <- FALSE
      warning("Empirical ICCs are not produced for different NA groups!\n")
    }} 
  if (!empirical) emp.plot <- FALSE
    
  if ((object$model != "RM") && (empirical)){
    warning("Empirical ICCs can only be plotted for a dichotomous Rasch model!\n")
    emp.plot <- FALSE
  }         
  
  theta <- seq(xlim[1],xlim[2],by=0.1)                          #x-axis
  p.list <- plist.internal(object,theta)                        #matrix of probabilities
  th.ord <- order(theta)
 
  if (any(item.subset=="all")) {
    textlab <- colnames(object$X)
    ivec <- 1:length(p.list)
  } else {
      if (is.character(item.subset)) {                         #item names specified
      ivectemp <- t(as.matrix(1:length(p.list)))
      colnames(ivectemp) <- colnames(object$X)
      ivec <- ivectemp[,item.subset]
      textlab <- item.subset
      textlab[ivec] <- textlab
      it.legend <- item.subset
    } else {                                                    #numeric vector specified
      textlab <- colnames(object$X)[item.subset]
      textlab[item.subset] <- textlab
      ivec <- item.subset
    }
  }
  
  if (object$model=="RM") {                                     #Rasch model
    p.list <- lapply(p.list,function(x) {x[,-1]})               #Delete 0-probabilites
    p.mat <- matrix(unlist(p.list),ncol=length(p.list))         #matrix with solving probabilities
    text.ylab <- p.mat[(1:length(theta))[theta==median(theta)],]
  }
  
    if (object$model != "RM"){ 
      for (i in ivec) {                                 #runs over items
         yp <- as.matrix(p.list[[i]])
         yy <- yp[th.ord,]
         get(getOption("device"))()
         matplot(sort(theta),yy,type="l",lty=1,col=1:(dim(yp)[2]),
                 main=paste("ICC plot for item ",textlab[i]),xlim=xlim,
                 ylim=ylim,xlab=xlab,ylab=ylab)
         legend(xlim[1],0.5,paste(c("Category"),0:(dim(yp)[2]-1)), col=1:(dim(yp)[2]),lty=1,lwd=1,...)
      }
    } else {
      if (any(item.subset=="all")) par(mfrow=c(2,2))
      for (i in ivec) {                                 #runs over items
         if (any(item.subset=="all")) {
           if ((i %% 4) == 0) {
              get(getOption("device"))()
              par(mfrow=c(2,2))
           }} else {get(getOption("device"))()}  
         yp <- as.matrix(p.list[[i]])
         yy <- yp[th.ord,]
         matplot(sort(theta),yy,type="l",lty=1,col=1:(dim(yp)[2]),
                 main=paste("ICC plot for item ",textlab[i]),xlim=xlim,
                 ylim=ylim,xlab=xlab,ylab=ylab,...)
         if (emp.plot) {
           freq.table <- as.matrix(table(rowSums(X),X[,i]))   
           rel.freq <- freq.table[,2]/rowSums(freq.table) 
           idx <- as.numeric(rownames(freq.table))
           lines(th.est$pred.list[[1]]$y[idx+1],rel.freq,type="l",...)
         }    
      }
    }
}


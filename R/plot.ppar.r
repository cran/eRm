plot.ppar <- function(x,xlab="Person Raw Scores",ylab="Person Parameters (Theta)",main=NULL,...)
# plot of the person raw scores against the person parameters
# x...object of class "ppar" (resulting from person.parameter.eRm)
{
  pl <- x$pred.list                              #list with spline interpolations

  if (is.null(pl)) stop("Spline interpolation required in person.parameter.eRm!")

  X <- x$X
  if (length(x$pers.ex) > 0) {
    X <- X[-x$pers.ex,]
    #gmemb <- x$gmemb[-x$pers.ex]
  }
  gmemb <- x$gmemb
  X.list <- split(as.data.frame(X),as.factor(gmemb))
  
  if (length(pl) > 1) {
    for (i in 1:length(pl)) main.text <- paste("Person Parameter Plot of Group",i)
  } else {
    main.text <- "Plot of the Person Parameters"
  }
  
  if (!is.null(main)) main.text <- main
  
  for (i in 1:length(pl)) {
    get(getOption("device"))()
    plot(rowSums(X.list[[i]],na.rm=TRUE),x$thetapar[[i]],xlim=c(min(pl[[i]]$x),max(pl[[i]]$x)),
         ylim=c(min(pl[[i]]$y),max(pl[[i]]$y)),xlab=xlab,ylab=ylab,
         main=main.text,...)
    lines(pl[[i]]$x,pl[[i]]$y)
  }
}



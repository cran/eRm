"plot.Rm" <-
function(x,y,...){                                           #graphical model check
if (length(x$G)!=1) {
   betapar1 <- x$W %*% x$likall123[[2]][[1]]$estimate
   betapar2 <- x$W %*% x$likall123[[3]][[1]]$estimate
   P <- max(abs(betapar1),abs(betapar2))
   plot(betapar1,betapar2,main="Graphical Model Test",xlab="Beta Group1",ylab="Beta Group2",ylim=c(-P,P),xlim=c(-P,P),type="n")
   text(betapar1,betapar2,(1:length(betapar1)))
   abline(0,1)
} else
  stop("Insufficient number of persons! No graphical model test can be produced.")
invisible(NULL)
}


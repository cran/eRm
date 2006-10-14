`plot.Rm` <-
function(x,y,...){                                           #graphical model check
if (length(x$LR)>1) {
   betapar1 <- x$W %*% x$etaparG1                            #beta parameters for group 1
   bp1sum0 <- scale(betapar1,scale=FALSE)                    #sum 0 standardization
   betapar2 <- x$W %*% x$etaparG2
   bp2sum0 <- scale(betapar2,scale=FALSE)
   
   P <- max(abs(bp1sum0),abs(bp2sum0))
   plot(bp1sum0,bp2sum0,main="Graphical Model Test",xlab="Beta Group1",ylab="Beta Group2",ylim=c(-P,P),xlim=c(-P,P),type="n")
   text(bp1sum0,bp2sum0,(1:length(bp1sum0)))
   abline(0,1)
} else
  stop("Insufficient number of subjects! No graphical model test can be produced.")
invisible(NULL)
}


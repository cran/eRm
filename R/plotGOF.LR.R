`plotGOF.LR` <-
function(x,beta.subset="all", xlab="Beta Group 1",ylab="Beta Group 2",
         ylim=c(-3,3),xlim=c(-3,3),type="p",...)
{                                           
# graphical model check
# beta.subset...plot only a subset of beta-parameters; either "all" or an index vector
# x...object of class LR (from LRtest)
# ...additional graphic parameters

if (length(x$likgroup) > 2) warning("Only the parameters for the first two subgroups are plotted!") 

nparg1 <- length(x$betalist[[1]])
nparg2 <- length(x$betalist[[2]])
if (nparg1 != nparg2) stop("Unequal number of parameters in the subgroups! Plot cannot be produced, choose another split in LRtest!")

beta1 <- x$betalist[[1]]
beta2 <- x$betalist[[2]]

if (is.character(beta.subset)) {
  if (beta.subset=="all") {
    beta.subset <- 1:length(beta1)
    textlab <- names(beta1)
  } else {
    textlab <- beta.subset
  }
} else {
  textlab <- names(beta1)[beta.subset]
}
yshift <- (ylim[2]-ylim[1])/30

plot(beta1[beta.subset],beta2[beta.subset],main="Graphical Model Check",xlab=xlab,
ylab=ylab,ylim=ylim,xlim=xlim,type=type,...)
text(beta1[beta.subset],beta2[beta.subset]+yshift,labels=textlab,...)
abline(0,1)

}


`plotGOF.LR` <-
function(x,beta.subset="all", xlab="Beta Group 1",ylab="Beta Group 2",tlab="item",
         ylim=c(-3,3),xlim=c(-3,3),type="p",pos="4",...)
{
# graphical model check
# beta.subset...plot only a subset of beta-parameters; either "all" or an index vector
# x...object of class LR (from LRtest)
# tlab ... labelling: "item" abbreviated beta parameter name, "number" number from beta par list, "none"
# pos (where the textlabel appears)
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
    #textlab <- names(beta1)
    switch(EXPR=tlab,
      item=textlab <- substr(names(beta1),6,100),  #remove "beta " from names
      number=textlab <- 1:length(beta1),
      identify=labs <- substr(names(beta1),6,100)
    )
  } else {
    textlab <- beta.subset
  }
} else {
  #textlab <- names(beta1)[beta.subset]
    switch(EXPR=tlab,
      item=textlab <- substr(names(beta1)[beta.subset],6,100),  #remove "beta " from names
      number=textlab <- beta.subset,
      identify=labs <- substr(names(beta1)[beta.subset],6,100)
    )
}


#yshift <- (ylim[2]-ylim[1])/30
yshift<-0

plot(beta1[beta.subset],beta2[beta.subset],main="Graphical Model Check",xlab=xlab,
ylab=ylab,ylim=ylim,xlim=xlim,type=type,...)
abline(0,1)
if(exists("textlab")) {
      text(beta1[beta.subset],beta2[beta.subset]+yshift,labels=textlab,pos=pos,...)
}
if(exists("labs")) {
      options(locatorBell = FALSE)
      xycoords<-cbind(beta1,beta2)
      nothing<-identify(xycoords,labels = labs,atpen=TRUE,offset=1)
}

}


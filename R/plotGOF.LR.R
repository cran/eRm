`plotGOF.LR` <-
function(x,beta.subset="all", main="Graphical Model Check", xlab=NULL,ylab=NULL,tlab="item",
         ylim=c(-3,3),xlim=c(-3,3),type="p",pos="4", conf=NULL, ctrline=NULL,...)
{
# graphical model check
# beta.subset...plot only a subset of beta-parameters; either "all" or an index vector
# x...object of class LR (from LRtest)
# tlab ... labelling: "item" abbreviated beta parameter name, "number" number from beta par list,
#            "identify" interactive, "none"
# pos ... (where the textlabel appears)
# conf ... confidence ellipses: NULL or
#              list(gamma=0.95, col="red", ia=TRUE, lty="dashed", which=all items in beta.subset)
# ctrline ... control lines (confidence bands): NULL or list(gamma=0.95,lty="solid", col="blue")
# ...     additional graphic parameters

if (length(x$likgroup) > 2) warning("Only the parameters for the first two subgroups are plotted!")


if (is.null(xlab)) xlab<-paste("Beta for Group: ",x$spl.gr[1],sep="")
if (is.null(ylab)) ylab<-paste("Beta for Group: ",x$spl.gr[2],sep="")

nparg1 <- length(x$betalist[[1]])
nparg2 <- length(x$betalist[[2]])
if (nparg1 != nparg2) stop("Unequal number of parameters in the subgroups! Plot cannot be produced, choose another split in LRtest!")

beta1 <- x$betalist[[1]] * -1  # -1 to obtain difficulty parameters
beta2 <- x$betalist[[2]] * -1



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
  ##beta.subset<-sort(beta.subset)

    switch(EXPR=tlab,
      item=textlab <- substr(names(beta1)[beta.subset],6,100),  #remove "beta " from names
      number=textlab <- beta.subset,
      identify=labs <- substr(names(beta1)[beta.subset],6,100)
    )
}


#yshift <- (ylim[2]-ylim[1])/30
yshift<-0

plot(beta1[beta.subset],beta2[beta.subset],main=main,xlab=xlab,
       ylab=ylab,ylim=ylim,xlim=xlim,type=type,...)
abline(0,1)
if(exists("textlab")) {
      text(beta1[beta.subset],beta2[beta.subset]+yshift,labels=textlab,pos=pos,...)
}
if(exists("labs")) {
      options(locatorBell = FALSE)
      xycoords <- cbind(beta1[beta.subset], beta2[beta.subset])
      nothing<-identify(xycoords,labels = labs,atpen=TRUE,offset=1)
}

# se's needed for ellipses and control lines

if(is.list(conf) || is.list(ctrline)){

   if(any(is.na(unlist(x$selist)))) {
      warning("Confidence ellipses or control lines cannot be plotted.\n  LR object without standard errors. Use option 'se=TRUE' in LRtest()")
      conf <- ctrline <- NULL
   } else {
      s1 <- x$selist[[1]]
      s2 <- x$selist[[2]]
      v1 <- s1^2
      v2 <- s2^2
      suspicious.se<-any(cbind(s1,s2)[beta.subset]>10)
      if(suspicious.se){
         warning("Suspicious size of standard error(s).\n  Check model specification, split criterion, data.")
      }
   }

   #if(any(abs(cbind(beta1,beta2)[beta.subset])>8)){
   #   warning("Suspicious size of parameter estimate(s).\n  Check model specification, split criterion, data.")
   #   if(is.null(conf)) conf$ia <- FALSE
}


# confidence ellipses

if(is.list(conf)){


    # (interactive) plot of confidence ellipses

    ## function ellipse() from package car
    ellipse <-
    function (center, shape, radius, center.pch = 19, center.cex = 1.5,
        segments = 51, add = TRUE, xlab = "", ylab = "", las = par("las"),
        col = palette()[2], lwd = 2, lty = 1, ...)
    {
        if (!(is.vector(center) && 2 == length(center)))
            stop("center must be a vector of length 2")
        if (!(is.matrix(shape) && all(2 == dim(shape))))
            stop("shape must be a 2 by 2 matrix")
        angles <- (0:segments) * 2 * pi/segments
        unit.circle <- cbind(cos(angles), sin(angles))
        ellipse <- t(center + radius * t(unit.circle %*% chol(shape)))
        if (add)
            lines(ellipse, col = col, lwd = lwd, lty = lty, ...)
        else plot(ellipse, xlab = xlab, ylab = ylab, type = "l",
            col = col, lwd = lwd, lty = lty, las = las, ...)
        if (center.pch)
            points(center[1], center[2], pch = center.pch, cex = center.cex,
                col = col)
    }

    # select items for which ellipses are drawn  ## rh 2011-05-31
    if(is.null(conf$which)) conf$which<-beta.subset#seq_along(beta.subset)
    ##conf$which <- sort(conf$which)
    if(!all(conf$which %in% beta.subset))
        stop("Incorrect item number(s) for which ellipses are to be drawn")
    if(is.null(conf$col)) {
        conf$c <- rep("red",length.out=length(beta1))
    } else if (!is.null(conf$which)){
##        conf$c <- rep(NA,length.out=length(beta.subset))
        conf$c <- rep(NA,length.out=length(conf$which))
        if (length(conf$c)!=length(conf$which))
           stop("which and col must have the same length in specification of conf")
        else
           conf$c[conf$which]<-conf$col
    }
    conf$col <- conf$c

    if(is.null(conf$gamma)) conf$gamma <- 0.95
    if(is.null(conf$lty)) conf$lty <- "dotted"
    if(is.null(conf$ia)) conf$ia <- FALSE

    z <- qnorm((conf$gamma+1)/2)

    ci1u <- beta1 + z*s1
    ci1l <- beta1 - z*s1
    ci2u <- beta2 + z*s2
    ci2l <- beta2 - z*s2



    if(conf$ia) {


         identifyEll <- function(x, y, ci1u, ci1l, ci2u,ci2l, v1, v2, conf, n=length(x), ...)
         ## source: example from help("identify")
         ## a function to use identify to select points, and overplot the
         ## points with a cofidence ellipse as they are selected
         {
             xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
             sel <- rep(FALSE, length(x)); res <- integer(0)
             while(sum(sel) < n) {
                 ans <- identify(x[!sel], y[!sel], n=1, plot=FALSE, ...)
                 if(!length(ans)) break
                 ans <- which(!sel)[ans]
                 i <- ans
            lines(rep(x[i],2),c(ci2u[i],ci2l[i]),col=conf$col[1], lty=conf$lty)
            lines(c(ci1u[i],ci1l[i]), rep(y[i],2),col=conf$col[1],lty=conf$lty)
            ellipse(center=c(x[i],y[i]),matrix(c(v1[i],0,0,v2[i]),2),z,segments=200,center.cex=0.5,lwd=1, col=conf$col[1])
                 #points(x[ans], y[ans], pch = pch)
                 sel[ans] <- TRUE
                 res <- c(res, ans)
             }
             #res
         }
         identifyEll(beta1[beta.subset],beta2[beta.subset],
                             ci1u[beta.subset], ci1l[beta.subset], ci2u[beta.subset], ci2l[beta.subset],
                             v1[beta.subset], v2[beta.subset], conf)
    } else {

         # non-interactive: plot of all ellipses at once

         x<-beta1
         y<-beta2
         for (i in beta.subset) {
            if(i %in% conf$which){
              lines(rep(x[i],2),c(ci2u[i],ci2l[i]),col=conf$col[i], lty=conf$lty)
              lines(c(ci1u[i],ci1l[i]), rep(y[i],2),col=conf$col[i],lty=conf$lty)
              ellipse(center=c(x[i],y[i]),matrix(c(v1[i],0,0,v2[i]),2),z,segments=200,center.cex=0.5,lwd=1, col=conf$col[i])
            }
         }
    }
}


# 95% control lines (Wright)

if(is.list(ctrline)){

    if(is.null(ctrline$gamma)) ctrline$gamma <- 0.95
    if(is.null(ctrline$col)) ctrline$col <- "blue"
    if(is.null(ctrline$lty)) ctrline$lty <- "solid"

    z <- qnorm((ctrline$gamma+1)/2)

    d<-(beta1+beta2)/2
    se.d<-sqrt(v1+v2)
    d<-sort(d)
    se.d<-se.d[order(d)]
    upperx<-d-z*se.d/2
    uppery<-d+z*se.d/2
    lines(upperx,uppery, col=ctrline$col, lty=ctrline$lty)
    lines(uppery,upperx, col=ctrline$col, lty=ctrline$lty)


}

}

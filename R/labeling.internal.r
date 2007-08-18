labeling.internal <- function(model,X,X01,W,etapar,betapar,mpoints,ngroups)
{
#labeling for W, eta, beta.

if (is.null(colnames(W))) {                             #eta labels
    names(etapar) <- paste("eta",1:dim(W)[2])
    colnames(W) <- names(etapar)
  } else {
    names(etapar) <- colnames(W)
  }

if (mpoints == 1) {                                     #no mpoints labels
  if ((model=="RM") || (model=="LLTM")) {               #no category labels
    betanames <- paste("beta",colnames(X))
  } else {
    indmt <- apply(X,2,max,na.rm=TRUE)
    catnames <- sequence(indmt)
    itnames <- rep(colnames(X),indmt)
    betanames <- paste("beta",paste(itnames,catnames,sep=".c"))
  }

} else {                                                                 #repeated measurement models
  indmt0 <- apply(X,2,max,na.rm=TRUE)
  indmt <- rep(apply(X,2,max,na.rm=TRUE),ngroups)
  catnames <- sequence(indmt)                                            #category names
  if (substr(colnames(X)[1],1,2)=="I1") {                                #if item names specified by user
    itemind <- rep(paste("I",1:(dim(X)[2]/mpoints),sep=""),mpoints)
    itnames <- rep(itemind,indmt0)
  } else {
    itnames <- rep(colnames(X),indmt0)
  }
  if (model == "LLTM") {
    icnames <- itnames
  } else {
    icnames <- paste(itnames,catnames,sep=".c")
  }
  t.lab <- paste("t",rep(1:mpoints,each=length(icnames)/mpoints),sep="") #time labels
  if (ngroups > 1) {
    g.lab <- rep(paste("g",rep(1:ngroups,each=length(icnames)/mpoints/ngroups),sep=""),mpoints)
    betanames <- paste(icnames,t.lab,g.lab)
  } else {
    betanames <- paste(icnames,t.lab)
  }
}


if (is.null(rownames(W))) {                      #no labels provided
    rownames(W) <- betanames
    names(betapar) <- betanames
  } else {
    names(betapar) <- rownames(W)
 }

list(W=W,etapar=etapar,betapar=betapar)
}
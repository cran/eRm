`datcheck` <-
function(X,W,mpoints,groupvec,model)
{
  if (is.data.frame(X))  {X <- as.matrix(X)}                  #X as data frame allowed

  if (is.null(colnames(X))) {                                 #determine item names
    if (mpoints > 1) {
      mpind <- paste("t",rep(1:mpoints,each=(dim(X)[2]/mpoints),1),sep="") #time points
      itemind <- paste("I",1:(dim(X)[2]/mpoints),sep="")
      colnames(X) <- paste(itemind,mpind)
    } else {
      colnames(X) <- paste("I",1:dim(X)[2],sep="")                         #item labels
  }}
  if (is.null(rownames(X))) rownames(X) <- paste("P",1:dim(X)[1],sep="")   #person labels

#----------------------- check groupvec --------------------------

  if ((length(groupvec) > 1) && (length(groupvec) != dim(X)[1])) {
    stop("Wrong specification of groupvec!")}

  if (min(groupvec)!=1) {
    stop("Group specification must start with 1!")}

  if (length(unique(groupvec))!=(max(groupvec))) {
    stop("Group vector is incorrectly specified (perhaps a category is missing)!")} # rh 2011-03-03

  if ((max(groupvec) > 1) && (mpoints==1)) {
    stop("Model not identifiable! Group contrasts can only be imposed for repeated measurement designs.") }

#  if ((length(groupvec) > 1) && any(is.na(X))) {
#    stop("Model with repeated measures, group specification and NAs cannot be computed!") }

#----------------------- check X --------------------------------
allna.vec <- apply(X,2,function(y) {all(is.na(y))})                 #eliminate items with all NA's
if (any(allna.vec)) {stop("There are items with full NA responses which must be deleted!")}

allna.vec <- apply(X,1,function(y) {all(is.na(y))})                 #eliminate items with all NA's
if (any(allna.vec)) {stop("There are persons with full NA responses which must be deleted!")}

allna.vec <- apply(X,1,function(y) {sum(is.na(y))})
if (any(allna.vec == (dim(X)[2]-1))) {stop("Subjects with only 1 valid response must be removed!")}

ri.min <- apply(X,2,min,na.rm=TRUE)                                 #if no 0 responses
if (any(ri.min > 0)) {
  cat("Warning message: The following items have no 0-responses: \n")
  cat(colnames(X)[ri.min>0],sep=", ")
  cat("\n")
  cat("Responses are shifted such that lowest category is 0. \n")
  cat("\n")
}
X <- t(apply(X,1,function(y) {y-ri.min}))                           #shift down to 0

ri <- apply(X,2,sum,na.rm=TRUE)                                     #item raw scores
n.NA <- colSums(apply(X,2,is.na))                                   #number of NA's per column
maxri <- (dim(X)[1]*(apply(X,2,max,na.rm=TRUE)))-n.NA               #maximum item raw scores with NA
TFcol <- ((ri==maxri) | (ri==0))
X.n <- X[,!TFcol]                                                   #new matrix with excluded items
item.ex <- (1:dim(X)[2])[TFcol]                                     #excluded items
if (length(item.ex) > 0) {
  if (mpoints == 1) {
    cat("Warning message: The following items were excluded due to complete 0/full responses: \n")
    cat(colnames(X)[item.ex],sep=", ")
    cat("\n")
  } else {
    cat("The following items show complete 0/full responses: \n")
    cat(colnames(X)[item.ex],sep=", ")
    cat("\n")
    stop("Estimation cannot be performed! Delete the correponding items for the other measurement points as well! \n")
}}

if ((model=="PCM") || (model=="LPCM")) {                         #check if there are missing categories for PCM (for RSM doesn't matter)
  tablist <- apply(X,2,function(x) list(as.vector(table(x))))
  tablen <- sapply(tablist,function(x) length(x[[1]]))
  xmax <- apply(X,2,max)+1
  indwrong <- which(tablen != xmax)
  if (length(indwrong) > 0) {
    cat("The following items do not have responses on each category: \n")
    cat(colnames(X)[indwrong],sep=", ")
    cat("\n")
    cat("Warning message: Estimation may not be feasible. Please check data matrix! \n")
    cat("\n")
  }
}


#-------------------------- ill conditioned for RM and LLTM --------------
if ((model=="RM") || (model=="LLTM")) {
  if (length(table(X.n)) != 2) stop("Dichotomous data matrix required!")
  k.t <- dim(X.n)[2]/mpoints                                    #check for each mpoint separately
  t.ind <- rep(1:mpoints,1,each=k.t)
  X.nlv <- split(t(X.n),t.ind)                                  #split X due to mpoints
  cn.lv <- split(colnames(X.n),t.ind)
  X.nl <- lapply(X.nlv,matrix,ncol=k.t,byrow=TRUE)
  for (i in 1:length(X.nl)) colnames(X.nl[[i]]) <- cn.lv[[i]]

  for (l in 1:mpoints) {                                       #check within mpoint
    X.nll <- X.nl[[l]]
    k <- ncol(X.nll)
    adj <- matrix(0,ncol=k,nrow=k)
    for (i in 1:k) for(j in 1:k) {
        adj[i,j]<- 1*any(X.nll[,i]> X.nll[,j],na.rm=TRUE)
    }
    cd <- component.dist(adj, connected = "strong")
    cm <- cd$membership
    cmp <- max(cm)
    if(cmp>1) {
         cmtab <- table(cm)
         maxcm.n <- as.numeric(names(cmtab)[cmtab!=max(cmtab)])
         suspcol <- (1:length(cm))[tapply(cm,1:length(cm),function(x) any(maxcm.n==x))]
         n.suspcol <- colnames(X.nll)[suspcol]
         cat("Suspicious items:",n.suspcol,"\n")
         stop("Estimation stopped due to ill-conditioned data matrix X!")
    }
}}
#----------------------- end ill-conditioned check -------------------------------

list(X=X.n,groupvec=groupvec)
}


`datcheck` <-
function(X,W,mpoints,groupvec,model)
{
  if (is.data.frame(X))  {X <- as.matrix(X)}                  #X as data frame allowed
    
  if (is.null(colnames(X))) {
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
    stop("Group vector is specified wrongly!")}
  
  if ((max(groupvec) > 1) && (mpoints==1)) {
    stop("Model not identifiable!") }
  
  if ((length(groupvec) > 1) && any(is.na(X))) {
    stop("Model with repeated measures, group specification and NAs cannot be computed!") }
  
#----------------------- check X --------------------------------
allna.vec <- apply(X,2,function(y) {all(is.na(y))})                 #eliminate items with all NA's
if (any(allna.vec)) {stop("There are items with full NA responses which must be deleted!")}

allna.vec <- apply(X,1,function(y) {all(is.na(y))})                 #eliminate items with all NA's
if (any(allna.vec)) {stop("There are persons with full NA responses which must be deleted!")}

ri.min <- apply(X,2,min,na.rm=TRUE)                                 #if no 0 responses
X <- t(apply(X,1,function(y) {y-ri.min}))

ri <- apply(X,2,sum,na.rm=TRUE)                                     #item raw scores
n.NA <- colSums(apply(X,2,is.na))                                   #number of NA's per column
maxri <- (dim(X)[1]*(apply(X,2,max,na.rm=TRUE)))-n.NA               #maximum item raw scores with NA
TFcol <- ((ri==maxri) | (ri==0))  
X.n <- X[,!TFcol]
item.ex <- (1:dim(X)[2])[TFcol]                                     #excluded items
if (length(item.ex) > 0) {
  cat("Warning message: The following items were excluded due to complete 0 or correct responses: \n")
  cat(colnames(X)[item.ex],sep=", ")
  cat("\n") 
  }   

#-------------------------- ill conditioned for RM --------------
if (model=="RM") {
  k <- ncol(X)
  adj <- matrix(0,nc=k,nr=k)
  for (i in 1:k) for(j in 1:k) {
      adj[i,j]<- 1*any(X[,i]> X[,j],na.rm=TRUE)
  }
  cd <- component.dist(adj, connected = "strong")
  cm <- cd$membership
  cmp <- max(cm)
  if(cmp>1) {
       cmtab <- table(cm)
       maxcm.n <- as.numeric(names(cmtab)[cmtab!=max(cmtab)])
       suspcol <- (1:length(cm))[tapply(cm,1:length(cm),function(x) any(maxcm.n==x))]
       cat("Suspicious colums in X:",suspcol,"\n")
       stop("Estimation stopped due to ill-conditioned data matrix X!")
  } 
}
#----------------------- end ill-conditioned check -------------------------------   

#check more measurement points for 0/full raw score
  
list(X=X.n,groupvec=groupvec)
}


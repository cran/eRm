"likLR" <-
function (X,W,mpoints,Groups,model)
{

W1 <- W

if ((model=="LLTM") | (model=="LRSM") | (model=="LPCM")) {              #for these models no LR-test can be computed
  G <- 1
  XlistX <- list(X)
} else
{
  npersons <- dim(X)[1]
  G1 <- rep(1,as.integer(npersons/2))
  G2 <- rep(2,(npersons-as.integer(npersons/2)))
  G <- c(G1,G2)
  Xlistvek <- split(X,G)                                   #list of subgroup matrices
  XlistX2 <- lapply(Xlistvek,function(x) {x <- matrix(x,ncol=dim(X)[2],byrow=TRUE)})  #list of subgroup matrices

  maxcheck <- sapply(XlistX2,function (x) {apply(x,2,max)})                           #check if LR-test can be applied
  dimmX <- nrow(maxcheck)*ncol(maxcheck)
  if (sum(maxcheck==t(maxcheck)[1,])< dimmX) {
    warning("Subgroup split for LR-test not possible")
    XlistX <- list(X)                                                                 #data matrix only
    G <- 1
  } else {
    XlistX <- c(list(X),XlistX2)}                                                     #full data & subgroup matrices
}

#parameter estimation
likall <- lapply(XlistX, function(xl,W=W1) {

                         if (model=="RM") Xprep <- datprep_RM(xl,W)
                         else if (model=="LLTM") Xprep <- datprep_LLTM(xl,W,mpoints,Groups)
                         else if (model=="RSM") Xprep <- datprep_RSM(xl,W)
                         else if (model=="PCM") Xprep <- datprep_PCM(xl,W)
                         else if (model=="LRSM") Xprep <- datprep_LRSM(xl,W,mpoints,Groups)
                         else if (model=="LPCM") Xprep <- datprep_LPCM(xl,W,mpoints,Groups)

                         Lprep <- cmlprep(Xprep$X01,Xprep$mt_vek,mpoints,Groups)
                         parest <- fitcml(Lprep$mt_ind,Lprep$nr,Lprep$x_mt,Lprep$rtot,Xprep$W)
                         list(parest,Xprep$W)
                         })

#likelihood ratio test
if (length(G) == 1) {
  LR <- NA
}else {
  LR <- lr(likall)
}

W1 <- likall[[1]][[2]]
rownames(W1) <- NULL
colnames(W1) <- NULL
                         
list(W=W1,likall=likall,LR=LR,G=G)                          #returns design matrix and results
}


`likLR` <-
function (X,W,mpoints,Groups,model)
{

W1 <- W

#LR-test
if ((model=="LLTM") | (model=="LRSM") | (model=="LPCM")) {              #for these models no LR-test can be computed
  G <- 1
  XlistX <- list(X)
} else
{ 
  XlistX2 <- list(NULL)                                    
  rv <- rowSums(X)                                         #person raw scores
  XlistX2[[1]] <- X[(rv > median(rv))==FALSE,]                            #group low r
  XlistX2[[2]] <- X[(rv > median(rv))==TRUE,]                             #group high r
  
  maxcheck <- sapply(XlistX2,function (x) {apply(x,2,max)})  #check if LR-test can be applied, i.e. if all items have at least 1 response
  dimmX <- nrow(maxcheck)*ncol(maxcheck)
  if (sum(maxcheck==t(maxcheck)[1,])< dimmX) {
    warning("Not enough subjects to perform LR-test!")
    XlistX <- list(X)                                                                 #data matrix only
    G <- 1
  } else {
    XlistX <- c(list(X),XlistX2)
    G <- 2
    }                                                     #full data & subgroup matrices
}

#parameter estimation
likall <- lapply(XlistX, function(xl,W=W1) {

                         if (model=="RM") Xprep <- datprep_RM(xl,W)
                         else if (model=="LLTM") Xprep <- datprep_LLTM(xl,W,mpoints,Groups)
                         else if (model=="RSM") Xprep <- datprep_RSM(xl,W)
                         else if (model=="PCM") Xprep <- datprep_PCM(xl,W)
                         else if (model=="LRSM") Xprep <- datprep_LRSM(xl,W,mpoints,Groups)
                         else if (model=="LPCM")  Xprep <- datprep_LPCM(xl,W,mpoints,Groups)

                         Lprep <- cmlprep(Xprep$X01,Xprep$mt_vek,mpoints,Groups,Xprep$W)                   
                         parest <- fitcml(Lprep$mt_ind,Lprep$nrlist,Lprep$x_mt,Lprep$rtot,Xprep$W,max(Groups),gind=Lprep$gind,x_mtlist=Lprep$x_mtlist)      
                         list(parest,Xprep$W)
                         })

#likelihood ratio test
if (G == 1) {
  LR <- NA
}else {
  LR <- lr(likall)
}

W1 <- likall[[1]][[2]]
rownames(W1) <- NULL
#evtl. noch if abfrage, falls namen spezifiziert wurden; evtl. eta-parameter
colnames(W1) <- NULL
                         
list(W=W1,likall=likall,LR=LR)                          #returns design matrix and results
}


`likLR` <-
function (X,W,mpoints,Groups,model,st.err,sum0)
{

#if (any(is.na(X))) {
 # pcX <- prelim_cat_pm(as.matrix(X))
 # rrow <- as.numeric(rownames(pcX$r))           #frequencies for different NA patterns
 # ogmemb <- rep(1:length(rrow),rrow)
 # gmemb <- ogmemb[pcX$ro]
#} else {
#  gmemb <- rep(1,dim(X)[1])
#}

if (any(is.na(X))) {
  dichX <- ifelse(is.na(X),1,0)
  strdata <- apply(dichX,1,function(x) {paste(x,collapse="")})
  gmemb <- as.vector(data.matrix(data.frame(strdata)))
} else {
  gmemb <- rep(1,dim(X)[1])
}

#parameter estimation
if (model=="RM") { Xprep <- datprep_RM(X,W,sum0)
} else if (model=="LLTM") { Xprep <- datprep_LLTM(X,W,mpoints,Groups,sum0)
} else if (model=="RSM") { Xprep <- datprep_RSM(X,W,sum0)
} else if (model=="PCM") { Xprep <- datprep_PCM(X,W,sum0)
} else if (model=="LRSM") { Xprep <- datprep_LRSM(X,W,mpoints,Groups,sum0)
} else if (model=="LPCM")  {Xprep <- datprep_LPCM(X,W,mpoints,Groups,sum0)
}

Lprep <- cmlprep(Xprep$X01,Xprep$mt_vek,mpoints,Groups,Xprep$W,gmemb)                   
parest <- fitcml(Lprep$mt_ind,Lprep$nrlist,Lprep$x_mt,Lprep$rtot,Xprep$W,
                 max(Groups),gind=Lprep$gind,x_mtlist=Lprep$x_mtlist,
                 Lprep$NAstruc,g_NA=Lprep$g_NA,st.err)      

W1 <- Xprep$W
rownames(W1) <- NULL
colnames(W1) <- NULL
options(warn=0)
                         
list(W=W1,parest=parest,X01=Xprep$X01)                          #returns design matrix and results
}


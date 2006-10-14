`LRSM` <-
function(X,W,mpoints=1,groupvec=1)
{
#library("gtools")
model <- "LRSM"
if (missing(W)) W <- NA
else W <- as.matrix(W)

XWcheck <- datcheck(X,W)                              #inital check of X and W
X <- XWcheck$X

lres <- likLR(X,W,mpoints,groupvec,model)
likall <- lres$likall[[1]]
LR <- lres$LR
                                
loglik <- -likall[[1]]$minimum                         #log-likelihood value
iter <- likall[[1]]$iterations                         #number of iterations
etapar <- likall[[1]]$estimate                         #eta estimates
se <- sqrt(diag(solve(likall[[1]]$hessian)))           #standard errors
betapar <- as.vector(lres$W%*% etapar)                 #beta estimates

result <- list(model=model,loglik=loglik,df=dim(lres$W)[2],iter=iter,etapar=etapar,se_eta=se,hessian=likall[[1]]$hessian,betapar=betapar,
               LR=LR,W=lres$W,mpoints=mpoints,ngroups=max(groupvec))
class(result) <- "eRm"                                 #classes: simple RM and extended RM
result
}


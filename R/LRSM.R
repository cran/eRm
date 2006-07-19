"LRSM" <-
function(X,W,mpoints=1,Groups=1)
{
#library("gtools")
model <- "LRSM"
if (missing(W)) W <- NA

XWcheck <- datcheck(X,W)                              #inital check of X and W
X <- XWcheck$X

lres <- likLR(X,W,mpoints,Groups,model)
likall <- lres$likall[[1]]
LR <- lres$LR
                                
loglik <- -likall[[1]]$minimum                         #log-likelihood value
iter <- likall[[1]]$iterations                         #number of iterations
etapar <- likall[[1]]$estimate                         #eta estimates
se <- sqrt(diag(solve(likall[[1]]$hessian)))           #standard errors
betapar <- as.vector(lres$W%*% etapar)                 #beta estimates

result <- list(loglik=loglik,iter=iter,etapar=etapar,se_eta=se,betapar=betapar,
               LR=LR,likall=likall,W=lres$W,mpoints=mpoints,ngroups=max(Groups))
class(result) <- c("Rm","eRm")                         #classes: simple RM and extended RM
result
}


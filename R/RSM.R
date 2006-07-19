"RSM" <-
function(X,W)
{
#...X: person*item scores matrix (starting from 0)

#-------------------main programm-------------------
#library("gtools")
Groups <- 1
mpoints <- 1
model <- "RSM"
if (missing(W)) W <- NA

XWcheck <- datcheck(X,W)                              #inital check of X and W
X <- XWcheck$X

lres <- likLR(X,W,mpoints,Groups,model)
likall123 <- lres$likall
likall <- lres$likall[[1]]
LR <- lres$LR
                                
loglik <- -likall[[1]]$minimum                         #log-likelihood value
iter <- likall[[1]]$iterations                         #number of iterations
etapar <- likall[[1]]$estimate                         #eta estimates
se <- sqrt(diag(solve(likall[[1]]$hessian)))           #standard errors
betapar <- as.vector(lres$W%*% etapar)                 #beta estimates

result <- list(loglik=loglik,iter=iter,etapar=etapar,se_eta=se,betapar=betapar,
               LR=lres$LR,likall=likall,W=lres$W,likall123=likall123,G=lres$G,
               mpoints=mpoints,ngroups=max(Groups))
class(result) <- c("Rm","eRm")                         #classes: simple RM and extended RM
result
}


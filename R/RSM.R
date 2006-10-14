`RSM` <-
function(X,W)
{
#...X: person*item scores matrix (starting from 0)

#-------------------main programm-------------------
groups <- 1
mpoints <- 1
model <- "RSM"
if (missing(W)) W <- NA
else W <- as.matrix(W)

XWcheck <- datcheck(X,W)                              #inital check of X and W
X <- XWcheck$X

lres <- likLR(X,W,mpoints,groups,model)
likall <- lres$likall[[1]]
LR <- lres$LR
                                
loglik <- -likall[[1]]$minimum                         #log-likelihood value
iter <- likall[[1]]$iterations                         #number of iterations
etapar <- likall[[1]]$estimate                         #eta estimates
se <- sqrt(diag(solve(likall[[1]]$hessian)))           #standard errors
betapar <- as.vector(lres$W%*% etapar)                 #beta estimates

if (length(LR)>1) {
  etaparG1 <- lres$likall[[2]][[1]]$estimate             #parameter vector group1 in LR-test (needed for plot method)
  etaparG2 <- lres$likall[[3]][[1]]$estimate             #parameter vector group2 
} else {
  etaparG1 <- NA
  etaparG2 <- NA
}


result <- list(model=model,loglik=loglik,df=dim(lres$W)[2],iter=iter,etapar=etapar,se_eta=se,hessian=likall[[1]]$hessian,betapar=betapar,
               LR=lres$LR,W=lres$W,etaparG1=etaparG1,etaparG2=etaparG2)

class(result) <- c("Rm","eRm")                         #classes: simple RM and extended RM
result
}


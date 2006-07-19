"RM" <-
function(X,W)
{
#...X: 0/1 person*item matrix

#-----------sourcing required files-----------------
#library("gtools")
#source("C:/erm/additional/datcheck.r")
#source("C:/erm/additional/datprep.r")
#source("C:/erm/additional/cmlprep.r")
#source("C:/erm/additional/cml.r")
#source("C:/erm/additional/lr.r")
#source("C:/erm/additional/likLR.r")

#----------------end sourcing-----------------------

#-------------------main programm-------------------
Groups <- 1
mpoints <- 1
model <- "RM"
if (missing(W)) W <- NA

XWcheck <- datcheck(X,W)                              #inital check of X and W
X <- XWcheck$X

lres <- likLR(X,W,mpoints,Groups,model)
likall123 <- lres$likall                              #values for subgroups for plot.eRm
likall <- lres$likall[[1]]                            #full groups for parameter estimation
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


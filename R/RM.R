`RM` <-
function(X, W, se = TRUE, sum0 = TRUE, etaStart)
{
#...X: 0/1 person*item matrix

#-------------------main programm-------------------
groupvec <- 1
mpoints <- 1
model <- "RM"

if (missing(W)) W <- NA
else W <- as.matrix(W)

if (missing(etaStart)) etaStart <- NA
else etaStart <- as.vector(etaStart)

XWcheck <- datcheck(X,W,mpoints,groupvec,model)                              #inital check of X and W
X <- XWcheck$X

lres <- likLR(X,W,mpoints,groupvec,model,st.err=se,sum0,etaStart)
parest <- lres$parest                             #full groups for parameter estimation
                                
loglik <- -parest$minimum                         #log-likelihood value
iter <- parest$iterations                         #number of iterations
etapar <- parest$estimate                         #eta estimates
if (se) {
  se <- sqrt(diag(solve(parest$hessian)))         #standard errors
} else {
  se <- rep(NA,length(etapar))
}
betapar <- as.vector(lres$W%*% etapar)            #beta estimates
X01 <- lres$X01 

npar <- dim(lres$W)[2]                            #number of parameters
N <- dim(X)[1]                                    #number of persons
AIC <- -2*loglik + 2*npar
BIC <- -2*loglik + log(N)*npar
cAIC <- -2*loglik + log(N)*npar + npar
IC <- list(AIC=AIC,BIC=BIC,cAIC=cAIC)

result <- list(X=X,X01=X01,model=model,loglik=loglik,IC=IC,npar=npar,iter=iter,
               etapar=etapar,se.eta=se,hessian=parest$hessian,betapar=betapar,W=lres$W)
class(result) <- c("dRm","Rm","eRm")                    #classes: dichotomous RM, RM (RM, PCM, RSM), and extended RM (all)
result
}


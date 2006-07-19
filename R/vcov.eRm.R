"vcov.eRm" <-
function(object,...) solve(object$likall[[1]]$hessian)      #VC-matrix of the parameter estimates


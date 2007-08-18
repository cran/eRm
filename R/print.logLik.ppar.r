`print.logLik.ppar` <-
function (x, digits = getOption("digits"),...)
{
    cat("'Unconditional (joint) log Lik.' ", format(c(x), digits = digits), " (df=",
        format(attr(x, "df")), ")\n", sep = "")
    invisible(x)
}
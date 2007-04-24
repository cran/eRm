`print.logLik.eRm` <-
function (x, digits = getOption("digits"),...)
{
    cat("'log Lik.' ", format(c(x), digits = digits), " (df=",
        format(attr(x, "df")), ")\n", sep = "")
    invisible(x)
}


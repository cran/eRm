# ==========================================================================
# use summary function to return percentage of misfitting person of object x
# ==========================================================================

summary.MisfittingPersons <- function(object, ...)
  # summary method for object of class "MisfittingPersons" 
  # x = object of class "MisfittingPersons"
{
  cat("\n")
  cat(" Percentage of Misfitting Persons:", "\n", round(object$PersonMisfit,4), "%", "\n")
  cat(" (i.e.,", object$count.misfit.Z, "out of",  object$total.persons, "persons have Chi-square based Z-values > 1.96)")
  cat("\n")
}


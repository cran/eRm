# ===============================================================
# use print function to return separation reliability of object x
# ===============================================================

print.MisfittingPersons <- function(x, ...)
  # print method for object of class "MisfittingPersons" 
  # x = object of class "MisfittingPersons"
{
  cat("\n")
  cat("Percentage of Misfitting Persons:", round(x$PersonMisfit,4), "%","\n")
  cat("\n")
}


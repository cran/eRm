PersonMisfit <- function(object) UseMethod("PersonMisfit")

PersonMisfit.ppar <- function (object) {    # takes person.parameter as object

  # extract person fit values   
  pfit <- personfit(object)
  person.misfit.Z <- pfit$p.infitZ
  person.misfit.Z[order(unlist(person.misfit.Z),decreasing=TRUE)]
  
  # count total number of persons and proportion of misfitting persons
  total.persons <- length(person.misfit.Z)
  count.misfit.Z<-sum(person.misfit.Z>1.96)
  PersonMisfit<-count.misfit.Z/total.persons*100
  
  # save result as a lsit
  result<-list(PersonMisfit=PersonMisfit, count.misfit.Z=count.misfit.Z, total.persons=total.persons)
  class(result) <- "MisfittingPersons"    # define the outcome of the function "PersonMisfit" as an object of class "MisfittingPersons" 
  return(result)
}

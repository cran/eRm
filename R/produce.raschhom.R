"produce.raschhom" <-
function(n.persons,n.items)
{
#produces rasch homogeneous data


person <- rnorm(n.persons)#draw person parameters from N(0,1)
item <- rnorm(n.items)#draw item parameters from N(0,1)

psolve <- matrix(NA,n.persons,n.items)#init matrix with solving probabilites

A <- outer(person,item,FUN="-")#without "for" loop
psolve <- exp(A)/(1+exp(A))

rpsolve <- matrix(runif(n.items*n.persons),n.persons,n.items)#draw random probabilities
X <- (rpsolve < psolve)*1#0/1-data matrix 

person.scores <- apply(X,1,sum)#person raw score
X <- X[order(person.scores),]#sorting X by person scores (for later raw score splitting)
person.scores <- sort(person.scores)#sorted person raw scores
 
item.scores <- apply(X,2,sum)#item raw score

list(X=X,person.scores=person.scores,item.scores=item.scores)
}


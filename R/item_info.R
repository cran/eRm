i_info<-function(hvec,itembeta,theta)
## calculates information (Samejima, 1969. Eq. 6-6 second line) for an item i as a function of theta
## Has internal function for getting the value and first and sewcond derivatives of the PCM for values of latent traits    
#@input: hvec...number of categories of item
#        itembeta...cumulative item parameters
#        theta ... supporting or sampling points on latent trait
#@output: a list with
#         $c.info...matrix of category information with categories in columns and values of theta (rows) (Samejima 1969, Eq. 6-6, first line)
#         $i.info...vector of item information at values of theta (Samejima 1969, Eq. 6-6, second line)
#@author: Thomas Rusch
#@date: 9.8.2023 
#
  {
   if(missing(theta)) theta <-seq(-5,5,0.01)
   p.ih<-function(hvec,itembeta,theta)
   #Calculates p.ih, first and second derviative of p.ih of a PCM for given item i
   #needs categories given (hvec) and the cumulative item parameters of the item (itembeta)
   #@output: a list with
   #         $p.ih...matrix of probabilities to fall into category h (colums) for given items as a function of theta (rows).
   #         $dp.ih ...the first derivative
   #         $d2p.ih ...the second derivative
  {
    beta <- c(0,itembeta) #eRm gives itempar with first fixed to zero
    numerator<-exp(outer(hvec,theta)+beta) #Numerator of the PCM, exp(h theta+b_ih) with each row being one h and columns being the theta values
    tmp<-hvec*numerator # h * exp(h theta+b_ih) (each row being a different h) and columns being the theta values
    tmp2<-hvec^2*numerator # h^2 * exp(h theta+b_ih) (each row being a different h) and columns being the theta values
    atheta <- apply(tmp,2,sum) #a(theta)=sum_l l exp(l*theta+beta_il) (so summed over categories) and number of elements being the thetas  
    btheta <- apply(numerator,2,sum) #b(theta)=sum_l exp(l*theta+beta_il) (so summed over categories)
    dtheta <- apply(tmp2,2,sum) #d(theta)=sum_l l^2 exp(l*theta+beta_il) (so summed over categories)
    pih<-t(numerator)/btheta   #now categories in column, thetas in rows
    dpih <- t(hvec*t(pih))- pih*(atheta/btheta) #first derivative: hPih - Pih a(theta)/b(theta)
    d2pih <- t(hvec^2*t(pih)) - 2*t(hvec*t(pih))*(atheta/btheta)+ 2*pih*(atheta^2/btheta^2)-pih*(dtheta/btheta) 
    return(list("pih"=pih,"dpih"=dpih,"d2pih"=d2pih))
  }
   tempus <- p.ih(hvec,itembeta,theta)
   c.info <- tempus$dpih^2/tempus$pih-tempus$d2pih #calculates category info (columns) for all theta(rows)
   i.info <-apply(c.info,1,sum) # calculates iteminfo for all theta(rows)
   return(list("c.info"=c.info,"i.info"=i.info))
 }

item_info <- function(ermobject,theta=seq(-5,5,0.01))
##Calculates information (Samejima, 1969) of all items as a function of the latent trait, theta
#        ermobject ... object of class eRm
#        theta ... supporting or sampling points on latent trait
#@output: a list where each element corresponds to an item and contains
#         $c.info...matrix of category information with categories in columns and values of theta (rows)
#         $i.info...vector of item information at values of theta
#@author: Thomas Rusch
#@date:13.6.2011
#
{
   vec.tmp <- get_item_cats(X=ermobject$X,nitems=dim(ermobject$X)[2],grp_n=dim(ermobject$X)[1])
   betapar <- ermobject$betapar
   veco <- unlist(lapply(vec.tmp,max))
   alloc.list<-vector("list",length(veco))
   hvec.list <- vector("list",length(veco))
   out.list <- vector("list",length(veco))
   for (i in 1:length(veco))
     {
       alloc.list[[i]] <- rep(i,veco[i])
       hvec.list[[i]] <- seq(0,veco[i])
     }
   uu<-unlist(alloc.list)
   itembeta.list <- split(betapar,uu)
   for (i in 1:length(itembeta.list))
     {
      out.list[[i]] <- i_info(hvec.list[[i]],itembeta.list[[i]],theta) #patch
    }
   return(out.list)
}



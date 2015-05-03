  #' Funkme_MetHast_Q 
#'
#' text...
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

formu2<-function(stateChain, Y, K=K){
      stateChain <-factor(stateChain, levels=c(1:K))
      # number in each state and sum of Ys
      xy<-as.data.frame(cbind(Y, stateChain))
      ny<-gapply(xy, FUN=function(x) length(x$Y), groups=stateChain)
      sumy<-gapply(xy, FUN=function(x) sum(x$Y), groups=stateChain)
      return(cbind(sumy, ny))}  

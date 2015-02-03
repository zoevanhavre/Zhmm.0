 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples #

 getq0NEW<-function(Q){
      require(markovchain)
 
      return(steadyStates(new("markovchain", states=as.character(c(1:dim(Q)[1])), transitionMatrix=Q)))
 }  
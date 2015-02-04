
 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

aminTheory<-function(amin, K){
	return((K-1)*(K-1+amin)*(1+( 1/( .5 -amin*(K-1) ) ))-amin*(K-1)+0.1)
	
}
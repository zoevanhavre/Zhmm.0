
 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

aMAXtheory<-function(amin, K, ref=c("Rousseau", "Gassiat")){
	if(ref=="Gassiat"){ simonSays<-K*((K-1)-1+1)-(K-1)*amin+0.01}
	if(ref=="Rousseau"){ simonSays<-(K-1)*(1+K-2+amin)*(1+(1/(.5-amin*(K-1)) ))-amin*(K-1)+0.01}

	return(simonSays)
	
}
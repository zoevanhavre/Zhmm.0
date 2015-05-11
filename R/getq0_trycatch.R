 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

 getq0_trycatch<-function(Q){
test1<-tryCatch(	getq0(Q), error= function(e) ALTERNATEq0(Q))
if(sum(test1>0)==length(test1)){
	return(test1)
}else{
	ALTERNATEq0(Q)
    }
}



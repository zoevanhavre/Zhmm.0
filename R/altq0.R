 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

ALTERNATEq0<-function(qnew){
	estq0<-round(qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew
			, 8)
			#return(estq0[ 1, (estq0 != 0)[1,] ])
			#return( estq0)
		okid<-c(1:dim(qnew)[1])[apply((estq0 == 1), 1, sum)==0]  
       	if(length(okid)==0){ okid<-1}                    
		return( estq0[okid[1], ])
}
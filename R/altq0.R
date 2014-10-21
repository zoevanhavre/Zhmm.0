 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

altq0<-function(qnew){ estq0<-round(qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew, 4);
#return(estq0[ 1, (estq0 != 0)[1,] ])
return( estq0)
}
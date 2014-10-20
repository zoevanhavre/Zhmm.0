 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))
 getq0<-function(Q){
      K<-dim(Q)[1]
      U<-matrix(rep(1/K, K*K), K,K)
      u<-rep(1/K, K)
      I<-diag(K)
      ipu<-I-Q+U
      u%*%solve(ipu)
      }  
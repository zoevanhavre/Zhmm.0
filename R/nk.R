 #' gibbsHMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

nk<-function(States, K=K){
      nkOrd<-matrix(nrow=dim(States)[1], ncol=K)
      for (i in 1:dim(States)[1]){
      rowS<-factor(States[i,], levels=1:K)
      nkOrd[i,]<-table(rowS)
      }
      return(nkOrd)
      }
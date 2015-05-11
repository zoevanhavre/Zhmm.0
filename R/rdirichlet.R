 #' gibbsHMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

 rdirichlet<-function(m=1,par){
          k=length(par);  mat=matrix(0,m,k)
          for (i in 1:m)  {sim=rgamma(k,shape=par,scale=1); mat[i,]=sim/sum(sim)}
#         mat<-sapply(mat, function(x) ifelse(x<0, 0, x ) )
          mat
          }
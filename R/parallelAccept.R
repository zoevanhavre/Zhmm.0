 #' UpdateStates
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))


parallelAccept<-function(w1, w2, a1, a2){
						w1[w1< 1e-200]<-1e-200             # truncate so super small values dont crash everyting
						w2[w2< 1e-200]<-1e-200
						T1<-dDirichlet(w2, a1, log=TRUE)
						T2<-dDirichlet(w1, a2, log=TRUE)
						B1<-dDirichlet(w1, a1, log=TRUE)
						B2<-dDirichlet(w2, a2, log=TRUE)
						MH<-min(1,	exp(T1+T2-B1-B2)) 
						Ax<-sample(c(1,0), 1, prob=c(MH,1-MH))
						return(Ax)}
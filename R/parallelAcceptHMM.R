 #' UpdateStates
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))


parallelAcceptHMM<-function(Qchain1, Qchain2, Alpha1, Alpha2){
						#w1[w1< 1e-200]<-1e-200             # truncate so super small values dont crash everyting
						#w2[w2< 1e-200]<-1e-200
Qchain1[Qchain1< 1e-200]<-1e-200             # truncate so super small values dont crash everyting
Qchain2[Qchain2< 1e-200]<-1e-200             # truncate so super small values dont crash everyting


K<-dim(Qchain1)[1]		

# Chain 2, prior 1
T1<-sum(sapply(c(1:K), function(x)dDirichlet(Qchain2[x,], Alpha1[x,], log=TRUE)))

#Chain 1, prior 2
T2<-sum(sapply(c(1:K), function(x)dDirichlet(Qchain1[x,], Alpha2[x,], log=TRUE)))

# Chain1 PRior 1
B1<-sum(sapply(c(1:K), function(x)dDirichlet(Qchain1[x,], Alpha1[x,], log=TRUE)))

# Chain 2, Prior 2
B2<-sum(sapply(c(1:K), function(x)dDirichlet(Qchain2[x,], Alpha2[x,], log=TRUE)))

# sum since logs??
	#T1<-dDirichlet(w2, a1, log=TRUE)
	#T2<-dDirichlet(w1, a2, log=TRUE)
	#B1<-dDirichlet(w1, a1, log=TRUE)
	#B2<-dDirichlet(w2, a2, log=TRUE)
	MH<-min(1,	exp(T1+T2-B1-B2)) 
	Ax<-sample(c(1,0), 1, prob=c(MH,1-MH))
	return(Ax)
					}


					
 #' UpdateStates
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))


parallelAcceptHMM_QV<-function(Qchain1, Qchain2, Alpha1, Alpha2, X1, X2){

q01<-ALTERNATEq0(Qchain1)
q02<-ALTERNATEq0(Qchain2)
q01[q01< 1e-200]<-1e-200  #?           # truncate so super small values dont crash everyting
q02[q02< 1e-200]<-1e-200 #? 

Qchain1[Qchain1< 1e-200]<-1e-200             # truncate so super small values dont crash everyting
Qchain2[Qchain2< 1e-200]<-1e-200             # truncate so super small values dont crash everyting


K<-dim(Qchain1)[1]		
X1<-as.numeric(table(factor(X1, levels=c(1:K))))
X2<-as.numeric(table(factor(X2, levels=c(1:K))))

X1[X1< 1e-200]<-1e-200  #?        
X2[X2< 1e-200]<-1e-200 	#?


# Chain 2, prior 1
T1<-dDirichlet(X2, q01, log=TRUE)+sum(sapply(c(1:K), function(x)dDirichlet(Qchain2[x,], Alpha1[x,], log=TRUE)))

#Chain 1, prior 2
T2<-dDirichlet(X1, q02, log=TRUE)+sum(sapply(c(1:K), function(x)dDirichlet(Qchain1[x,], Alpha2[x,], log=TRUE)))

topRatio<-T1+T2

# Chain1 PRior 1
B1<-dDirichlet(X1, q01, log=TRUE)+sum(sapply(c(1:K), function(x)dDirichlet(Qchain1[x,], Alpha1[x,], log=TRUE)))


# Chain 2, Prior 2
B2<-dDirichlet(X2, q02, log=TRUE)+sum(sapply(c(1:K), function(x)dDirichlet(Qchain2[x,], Alpha2[x,], log=TRUE)))

bottomRatio<-B1+B2
# sum since logs
	MH<-min(1,	exp(topRatio-bottomRatio)) 
	Ax<-sample(c(1,0), 1, prob=c(MH,1-MH))
	return(Ax)
					}


					
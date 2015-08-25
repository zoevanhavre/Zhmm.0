 #' Prior Parallel Tempering acceptance ratio
#'
#' acceptance ratio 
#' @param Qchain1, Qchain2, Alpha1, Alpha2,  X1all, X2all
#' @keywords dirichlet
#' @export
#' @examples #TBA

parallelAcceptHMM_Aug15<-function(Qchain1, Qchain2, Alpha1, Alpha2){
			  K<-dim(Alpha1)[1]
              Qchain1<-matrix(Qchain1, nrow=K, byrow=TRUE)
              Qchain2<- matrix(Qchain2, nrow=K, byrow=TRUE)
q01<-ALTERNATEq0(Qchain1)
q02<-ALTERNATEq0(Qchain2)
q01[q01< 1e-200]<-1e-200  #?           # truncate so super small values dont crash everyting
q02[q02< 1e-200]<-1e-200 #? 

Qchain1[Qchain1< 1e-200]<-1e-200             # truncate so super small values dont crash everyting
Qchain2[Qchain2< 1e-200]<-1e-200             # truncate so super small values dont crash everyting

K<-dim(Qchain1)[1]		

# Chain 2, prior 1
Q21<-sum(sapply(c(1:K), function(x)dDirichlet(Qchain2[x,], Alpha1[x,], log=TRUE)))

#Chain 1, prior 2
Q12<-sum(sapply(c(1:K), function(x)dDirichlet(Qchain1[x,], Alpha2[x,], log=TRUE)))

topRatio<-Q12+Q21

# Chain1 PRior 1
Q11<-sum(sapply(c(1:K), function(x)dDirichlet(Qchain1[x,], Alpha1[x,], log=TRUE)))

# Chain 2, Prior 2
Q22<-sum(sapply(c(1:K), function(x)dDirichlet(Qchain2[x,], Alpha2[x,], log=TRUE)))

bottomRatio<-Q22+Q11

# sum since logs
	MH<-min(1, exp(topRatio-bottomRatio)) 
	Ax<-sample(c(1,0), 1, prob=c(MH,1-MH))
	return(Ax)
					}



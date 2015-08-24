 #' Prior Parallel Tempering acceptance ratio
#'
#' acceptance ratio 
#' @param Qchain1, Qchain2, Alpha1, Alpha2,  X1all, X2all
#' @keywords dirichlet
#' @export
#' @examples #TBA

parallelAcceptHMM_Aug15<-function(Qchain1, Qchain2, Alpha1, Alpha2,  X1all, X2all){
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

C21<-q01[X2all[1]]
for( yo in 2:length(X2all)){C21*Qchain1[X2all[yo],X2all[yo]]}
Q21<-sum(sapply(c(1:K), function(x)dDirichlet(Qchain2[x,], Alpha1[x,], log=TRUE)))


#Chain 1, prior 2
C12<-q02[X1all[1]]
for( yo in 2:length(X1all)){C12*Qchain2[X1all[yo],X1all[yo]]}
Q12<-sum(sapply(c(1:K), function(x)dDirichlet(Qchain1[x,], Alpha2[x,], log=TRUE)))

topRatio<-log(C12)+Q12+log(C21)+Q21

# Chain1 PRior 1

C11<-q01[X1all[1]]
for( yo in 2:length(X1all)){C11*Qchain1[X1all[yo],X1all[yo]]}
Q11<-sum(sapply(c(1:K), function(x)dDirichlet(Qchain1[x,], Alpha1[x,], log=TRUE)))

# Chain 2, Prior 2
C22<-q02[X2all[1]]
for( yo in 2:length(X2all)){C22*Qchain2[X2all[yo],X2all[yo]]}
Q22<-sum(sapply(c(1:K), function(x)dDirichlet(Qchain2[x,], Alpha2[x,], log=TRUE)))

bottomRatio<-log(C22)+Q22+log(C11)+Q11

# sum since logs
	MH<-min(1, exp(topRatio-bottomRatio)) 
	Ax<-sample(c(1,0), 1, prob=c(MH,1-MH))
	return(Ax)
					}



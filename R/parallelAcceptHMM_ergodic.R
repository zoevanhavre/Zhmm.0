 #' UpdateStates
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))


parallelAcceptHMM_ergodic<-function(Qchain1, Qchain2, Alpha1, Alpha2){
q01<-ALTERNATEq0(Qchain1)
q02<-ALTERNATEq0(Qchain2)
q01[q01< 1e-200]<-1e-200  #?           # truncate so super small values dont crash everyting
q02[q02< 1e-200]<-1e-200 #? 

K<-dim(Qchain1)[1]		
# Chain 2, prior 1
T1<-dDirichlet(q02, Alpha1, log=TRUE)

#Chain 1, prior 2
T2<-dDirichlet(q01, Alpha2, log=TRUE)

topRatio<-T1+T2

# Chain1 PRior 1
B1<-dDirichlet(q01, Alpha1, log=TRUE)


# Chain 2, Prior 2
B2<-dDirichlet(q02, Alpha2, log=TRUE)

bottomRatio<-B1+B2
# sum since logs
	MH<-min(1,	exp(topRatio-bottomRatio)) 
	Ax<-sample(c(1,0), 1, prob=c(MH,1-MH))
	return(Ax)
					}


					
 #' gibbsHMM_PT
#'
#' parallel tempering with a column prior - option to mix over column or stick to j=1
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

ReplicateSimer2<-function(  N, n, Kfit=10, SimID, ITERATIONS,BURN,  AMAX,  PRIOR_TYPE, PTchain=20){
		#  STORE SIMULATIONS in a list
		simFunctionMorpher<-function(SimNumber){
			if(	SimNumber==1){ 	return(FunkSim1)
			}else if (SimNumber==2){	return(FunkSim3)
			}else if (SimNumber==3){	return(FunkSim4)
			}	}
		MorphingSIMULATE<-simFunctionMorpher(SimID)
		SIMS<-lapply( rep(n,N),  MorphingSIMULATE )
		# Compute density for L1 norm and store in a list
		simDensityMorpher<-function(SimNumber){
			if(	SimNumber==1){ 	return( SimDensity1)
			}else if (SimNumber==2){	return(SimDensity3)
			}else if (SimNumber==3){	return(SimDensity4)
			}	}
		MorphineDENSITY<-simDensityMorpher(SimID)

		SIM_DENSITY_TRUE<-lapply(SIMS, 	MorphineDENSITY)

	NumCores<-min(parallel::detectCores(), N)

		# Clean up Gibbs for lyra...
		library(parallel)
		Result<-mclapply(c(1:N), function(x)  {
		gibbsHMM_PT_wDist_LYRAfinally(YZ=SIMS[[x]],K=Kfit, densTrue=SIM_DENSITY_TRUE[[x]],  M=ITERATIONS,  alphaMAX=AMAX, type= PRIOR_TYPE, alphaMin=0.001, J=PTchain, SuppressAll="TRUE")

		} , mc.cores = NumCores)
		print(NumCores)

# combine results!
#Alive<-sapply(Result, function(x)  median(x$K0[-c(1:BURN)]))
Fin<-data.frame("Run"=rep(0,N), "MedianK0"=rep(0, N), "MeanDist"=rep(0,N))
for(i in 1:N){
	Fin[i,1]<-i
	.result<-Result[[i]]
#	print(head(.result))
	Fin[i,2]<- median(.result$K0[-c(1:BURN)])
	Fin[i,3]<-mean(.result$f2Dist[-c(1:BURN)])
}

 #Alive<-sapply(Result, function(x)  median( x[[]]$K0[-c(1:BURN)]))
# Alive<-sapply(Result, function(x)  median( x$K0[-c(1:BURN)]))
# L1norm<-sapply(Result, function(x)  mean(x[['f2Dist']][-c(1:BURN)]))
# SmallResults<-data.frame("AliveStates"=Alive, "L1norm"=L1norm)


		return(Fin)
		}




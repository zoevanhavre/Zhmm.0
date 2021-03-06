 #' gibbsHMM_PT 
#'
#' parallel tempering with a column prior - option to mix over column or stick to j=1
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

ReplicateSimer<-function(   N, n, SimID, ITERATIONS, BURN, AMAX,  PRIOR_TYPE, PTchain=20){
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
		
		
		# Clean up Gibbs for lyra...
		library(parallel)
		Result<-mclapply( c(1:N), function(x)  {
		gibbsHMM_LYRA( SIMS[[x]], SIM_DENSITY_TRUE[[x]],  M=ITERATIONS, burn=BURN, alphaMAX=AMAX, type= PRIOR_TYPE, alphaMin=0.001, J=PTchain)
		
		} )
		return(Result)
		}



	
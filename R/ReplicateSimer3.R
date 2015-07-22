 #' ReplicateSimer2
#'
#' parallel tempering with a column prior - option to mix over column or stick to j=1
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

ReplicateSimer3<-function(  N, n, Kfit=10, SimID, ITERATIONS, BURN,  AMAX, AMIN, PRIOR_TYPE,GridSize=50,  PTchain=30){
		#  STORE SIMULATIONS in a list
		simFunctionMorpher<-function(SimNumber){
			if(	SimNumber==1){ 	return(FunkSim1)
			}else if (SimNumber==2){	return(FunkSim2)
			}else if (SimNumber==3){	return(FunkSim4)
			}	}
		MorphingSIMULATE<-simFunctionMorpher(SimID)
		SIMS<-lapply( rep(n,N),  MorphingSIMULATE )

		# Clean up Gibbs for lyra...
# TODO Update inputs to match
Result.store<-data.frame("Replicate"=c(1:N), "SimID"=SimID, "n"=n,"AlphaMax"=AMAX,"alphaMin"=AMIN, "Prior"=PRIOR_TYPE, "ModeK0"=0, "MeanfDist"=0,  "WorstMixedMean"=0, "WorstMixedMin"=0)
BestModel<-vector("list", N)
for (.rep in 1:N){

#TODO # TV values for each simulation
if(	SimID==1){ 
	TVsim<-list(
		"z"=SIMS[[.rep]]$S[1:2],
		"q0"=ALTERNATEq0(matrix(  c(.2,0.3,0.5,    0.5,0.25,0.25,    0.25, 0.65, 0.1), nrow=3, byrow=T)),
		"Q"= c(  0.2,0.3,0.5,    0.5,0.25,0.25,    0.25, 0.65, 0.1))
}else if (SimID==2){	
	TVsim<-list(
		"z"=SIMS[[.rep]]$S[1:2],
		"q0"=ALTERNATEq0(matrix( c(  0.8,0.1,0.1,    0.2,0.75,0.05,    0.1, 0.3, 0.6), nrow=3, byrow=T)),
		"Q"= c(  0.8,0.1,0.1,    0.2,0.75,0.05,    0.1, 0.3, 0.6))
}else if (SimID==3){
	q1<-matrix(0.1, ncol=5, nrow=5) ; diag(q1)<-.6 ; q1[5,]<-c(    .225,.225,.225, .225, 0.1)
	TVsim<-list(
		"z"=SIMS[[.rep]]$S[1:2],
		"q0"=ALTERNATEq0(q1),
		"Q"= c(t(q1)))
			}

My.Result<-gibbsHMM_Main(YZ=SIMS[[.rep]],K=Kfit,   M=ITERATIONS,  alphaMAX=AMAX, type= PRIOR_TYPE, alphaMin=AMIN, J=PTchain, SuppressAll="FALSE", TV=TVsim, gridN=GridSize)
My.Result.PP<-Zhmm_PP( My.Result , burn=BURN, prep=100, isSim=TRUE, simlabel=paste( "Sim" ,SimID,"n",n, "Prior", PRIOR_TYPE, "Amax", AMAX,"Rep",.rep, "of", N, sep=""))

# pars bext model 
BestModel[[.rep]]<- list("Parameters"=My.Result.PP[[2]][which.max(My.Result.PP[[1]][,2])],
"States"=My.Result.PP[[3]][which.max(My.Result.PP[[1]][,2])])


Result.store$ModeK0[.rep]<-as.numeric(names(sort(table(factor(My.Result$Track$K0[-c(1:BURN)])),decreasing=TRUE)[1]))
Result.store$MeanfDist[.rep]<-mean(My.Result$Track$f2Dist[-c(1:BURN)])
Result.store$WorstMixedMean[.rep]<-mean(My.Result$Track$WorstMixProp[-c(1:BURN)])
Result.store$WorstMixedMin[.rep]<-min(My.Result$Track$WorstMixProp[-c(1:BURN)])

write.csv(Result.store[1:.rep,], file=paste( "Sim" ,SimID,"n",n, "Rep", N, "Prior", PRIOR_TYPE, "Amax", AMAX,"Amin", AMIN,".csv", sep=""))
save(Result.store, file=paste( "Sim" ,SimID,"n",n, "Rep", N, "Prior", PRIOR_TYPE, "Amax", AMAX,"Amin", AMIN, ".RDATA", sep="_"))

Sys.sleep(0.1)
print(Result.store[1:.rep,])
Sys.sleep(0.1)
}

pdf( file= paste( "Sim" ,SimID,"n",n, "Prior", PRIOR_TYPE, "MaxAlpha", AMAX,"Amin", AMIN,"Iters",ITERATIONS, ".pdf", sep="") ,width=6, height=3, pointsize=8)
	 		print( wq::layOut(
	 		list(ggplot(data=Result.store, aes(x=ModeK0))+geom_histogram(binwidth=1)+theme_bw()+xlab("K_0")+
	 			ggtitle(paste( "Sim" ,SimID ,"(n=" ,n, ") Prior", PRIOR_TYPE,  "Alpha", round( AMAX,3) , ": K_0"))  ,	1 ,	 1:2),
		    list(ggplot(data=Result.store, aes(y=MeanfDist, x= factor(1)))+geom_boxplot()+theme_bw()+ylab("Distance") ,	1,	 3)  )    )

dev.off()

	return(list(Result.store, BestModel))
		}




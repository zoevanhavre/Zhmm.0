 #' UpdateStates
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))


UpdateStatesNewINIT<-function( Y=Y, trans=qok, .mu,m, SumPerK, iniDist=c("PropObs", "UniIni") ){

  K<-sqrt(length(trans))
  n<-length(Y)
  T<-n+1  # since we want 0-T states estimated
NumberAlive<-sum(SumPerK>0)
PropPerGroup<-SumPerK/n
                 sampleZ<-rep(0, T)
                  trans<- matrix(c(trans), ncol=K, byrow=TRUE)
                  const<-rep(0,n)
                  OneStep<-matrix(nrow=T, ncol=K)
                  FilterP<-matrix(nrow=T-1, ncol=K)
                  Smooth<-matrix(nrow=T, ncol=K)

            # 2.1 Run Forward Filter to get P(St=j| yt, .) for j=1,...K and ti=1,...,T obs  (calling time var ti as t is used)
            #     Obtain P(Z(t)=l|Y(t-1),.) 1 step ahead pred of Z(ti)
            #for t=1 
                #compute P(Z(1)=l|t(0),.) using initial distribution of statesp
                #P(Z(1)=1|Y(0), .) 								  #   for (i in 1:K) OneStep[1,i]<- sum(initS[1:K] * trans[i,])
              
# Uniform?
                  if(iniDist=="UniIni"){
                     initS<-rep(1/NumberAlive, K) 
                     initS[PropPerGroup==0]<-0 } else if(iniDist=="PropObs"){
                   # equal to proportion n allocated
                  initS<-PropPerGroup}
###############            
# FOR t=1:

              OneStep[1,]<- initS

            #  OneStep[1,]<- apply( trans, 1, function(x)  sum( initS[1:K]*x)  )
              const[1]<-sum(dnorm(Y[1],mean=.mu)*OneStep[1,])
              
              #const[1]<-sum(dnorm(Y[1],mean=.mu[1,1:K])*OneStep[1,])

              FilterP[1,]<- (dnorm(Y[1],mean=.mu)*OneStep[1,]) /const[1]   #Obtain filtered probs P(Z(1)|Y(1))
                    
                  
              #for ti=2:...,
              for (ti in 2:n){					  #  for (i in 1:K) OneStep[ti,i]<-  sum(trans[i,]*FilterP[ti-1, ])  					
                OneStep[ti,]<- apply( trans, 1, function(x)  sum( x*FilterP[ti-1, ])  )
                const[ti]<-sum( dnorm(Y[ti],mean=.mu)*OneStep[ti,]) 		                    #Filter:  P(Z(t)|Y(t))
                FilterP[ti,]<- (dnorm(Y[ti],mean=.mu)*OneStep[ti,]) /const[ti]   
                }

              # 2.2 Sample final state Z(T) from P(Z(T)=j| y(T),.) (filter from t=100)
              Smooth[T,]<-FilterP[T-1,]
              sampleZ[T]<-sample(c(1:K), size=1, prob=Smooth[T,])

           for (ti in (T-1):1){
            	Smooth[ti,]<- (trans[sampleZ[ti+1],]*FilterP[ti, ]) / sum(trans[sampleZ[ti+1],]*FilterP[ti,])  	
              sampleZ[ti]<-sample(c(1:K), size=1, prob=Smooth[ti,])
                  		} 
	              
	           Z<-sampleZ
	  	MAP<-sum(log(const))
		
		return(list( "Z"=Z, "MAP"=MAP))	}
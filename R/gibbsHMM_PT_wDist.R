 #' gibbsHMM_PT 
#'
#' parallel tempering with a column prior - option to mix over column or stick to j=1
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))



gibbsHMM_PT_wDist<-function(YZ, M=2000, K=10 ,alphaMAX=1, type= 1, alphaMin=0.001, J=20, lab="sim", densTrue){
    #____SET UP_________________________________________
    ifelse(class(YZ)=='data.frame',    Y<-YZ$Observed, Y<-YZ)
    n=length(Y) # sample size
    varknown<-1 # known variace 
    mu0=0; var0=100
     # INITIALIZE
   # startVal<-makeStart(Y, K);  states0<-startVal$states0   #FUNK
    states0<-replicate(J, list())
    for(j in 1:J){states0[[j]]<-makeStartSimpler(Y, K)}

     # J= number of chains
        TrackParallelTemp<-matrix(nrow=M, ncol=J)
         TrackParallelTemp[1,]<-c(1:J)
      # TO BE INCORPORATED INTO J LISTS
    MU<-replicate(J,  matrix(nrow=M, ncol=K),  simplify=F)
    Q  <-replicate(J, matrix(nrow=M, ncol=K*K) , simplify=F)
    Qold<- replicate(J, diag(K), simplify=F)
    q0 <-replicate(J, matrix(nrow=M, ncol=K), simplify=F)
    Z  <-replicate(J, matrix(nrow=M, ncol=n+1)  ,  simplify=F)  #include 0 for initial state to be estimated too?
    SteadyScore<-data.frame("Iteration"=c(1:M), "K0"=K) ##### THIS IS NEW##
    PTsuccess  <- data.frame("Chain"=c(1:J),"Tries"=0,"Success"=0, "Ratio"=NA) #include 0 for initial state to be estimated too?
    ntSTORE<-replicate(J, list())
    K0Final<-matrix(nrow=M, ncol=J)
    MAP<-c(1:M)  # KEEP TARGET ONLY
    f2now<-vector(length=M)
    fDist<-vector(length=M)
    # ALPHA

    #alphaMAX<-(K-1)*(1+K-2+alphaMin)*(1+1/( (1/2) - alphaMin*(K-1))) -(K-1)*alphaMin+0.1
    Alpha_lows<-c(alphaMAX, exp(seq(log(alphaMAX), log(alphaMin), length=J))[-1])
    #Store alphs for PT
    STORE_Alphas<-replicate(J, list())
    pb <- txtProgressBar(min = 0, max = M, style = 3)
         
      #names(TrackParallelTemp)<-   AllAlphas[,2]
    # functions
    for (m in 1:M){ 
          
          if(m %% 100==0){Sys.sleep(0.1)
    #        print(PTsuccess)
          setTxtProgressBar(pb, m)
                  if(M < 20001){    
      #     par(mfrow=c(1,4))
      #    plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l')
      #    ts.plot(q0[[J]], main='q0 from target posterior', col=rainbow(K))
      #    ts.plot(TrackParallelTemp, main='Track Parallel Tempering', col=rainbow(J), gpars=list(yaxt="n") )
        #  axis(2, at=1:J, tick=1:J, labels=round(AllAlphas[,2],4), las=2) 
      #    axis(2, at=1:J, tick=1:J, labels=round(c(Alpha_lows),4), las=2)
      #    image(Z[[J]][,order(Y)], col=rainbow(K), main="Allocations vs ordered Y")
       
     par(mfrow=c(1,4))
    #        plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l')
    #        plot(as.vector(q0[[J]][c(m-99:m), ]), as.vector(MU[[J]][c(m-99:m), ]))
        slices <- prop.table(table((factor(SteadyScore$K0[-(1:m/5)], levels=c(1:K)))))
      barplot(slices, ylim=c(0,1), main="Number of non-empty states", xlab="Number of non-empty states", ylab="Probability (from MCMC)") 
        abline(h=seq(0, 1, .05), lwd=0.5, col='LightGrey')
         #   plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l')
          plot(as.vector(q0[[J]][c(m-99:m), ]), as.vector(MU[[J]][c(m-99:m), ]), main="Latest Posterior Draws (100 iter)", ylab="Stationary Dist", xlab="Means")

              #, color=rgb(0,0,0,alpha=as.vector(q0[[J]])))
           # image(Z[[J]][,order(Y)], col=rainbow(K), main="Allocations vs ordered Y")
            ts.plot(TrackParallelTemp, main='Track Parallel Tempering', col=rainbow(J), gpars=list(yaxt="n") )
            axis(2, at=1:J, tick=1:J, labels=round(c(Alpha_lows),4), las=2)
            ts.plot(fDist, main='L1 norm distance (y1, y2)')
       
          Sys.sleep(0)}}
          
                        
      for (j in 1:J){ # FOR EACH CHAIN           
          # FOR EACH CHAIN...
             # make matrix of alphas
 AllAlphas<-matrix(Alpha_lows[j],ncol=K, nrow=K)   

if (type==1){
    AllAlphas[,1]<-alphaMAX      # make said column Amax
}else if (type=="diag"){
   diag(AllAlphas)<-alphaMAX 
}else if (type=="mix"){
    if(sample( c(1, 0), size=1, prob=c(0.5, 0.5))==1){   # Put on diagonal
      diag(AllAlphas)<-alphaMAX 
    } else {  
      AllAlphas[,1]<-alphaMAX
    }}
   

            #  if(runif(1)>0.5){        AllAlphas[,1]<-alphaMAX      
           #     }else{                    diag(AllAlphas)<-alphaMAX}
                
          
                              # 1 Parameters given states Z(m-1)
                      # 1.1  Transition matrix Q from conditional posterior
                      if (m==1) {nt<-CountTrans(states0[[j]], K)
                        } else { nt<-CountTrans(Z[[j]][m-1,],K)}   # HERE ACCESS STATES
                          ntSTORE[[j]]<-nt
                        # draw transition probs for state 1:K
                    qnew<-matrix(ncol=K, nrow=K)
                    for(k in 1:K){qnew[k,]<-rdirichlet(par=nt[k,]+AllAlphas[k,])}
STORE_Alphas[[j]]<-AllAlphas
    q0new <-  ALTERNATEq0(qnew)  

                        #METROPOLIS Hastings STEP     
              if (m>1){   
                    q0Previous<-ALTERNATEq0(Qold[[j]])
                    A<-q0new[Z[[j]][m-1,1]]/q0Previous[Z[[j]][m-1,1]]   
                    # A<-q0new[Z[[j]][m-1,1]]/q0[[j]][m-1,Z[[j]][m-1,1]]   
                                                            if(A=='NaN'){A<-0}     
                    U<-runif(1,c(0,0.99))
              if(A>runif(1,c(0,0.99))){  #  Accept new values
                            Q[[j]][m,]<-as.vector(t(qnew))
                            q0[[j]][m,]<-q0new
              } else {  #Reject, chose OLD values of Q
                            #Q[[j]][m,]<-as.vector(t(Q[[j]][m-1,]))
                            Q[[j]][m,]<-as.vector(t(Qold[[j]]))

                            #q0[[j]][m,]<-q0[[j]][m-1,]  
                            q0[[j]][m,]<-q0Previous
                            }
              }else{ Q[[j]][m,]<-as.vector(t(qnew)) # 1st iteration always approved.
                          q0[[j]][m,]<-q0new  }

# new SAVE Qold for next iter        # **NEW**
Qold[[j]]<-  matrix( Q[[j]][m,]  , K,K, byrow=TRUE)                     # **NEW**


                        # 1.2 Update mu's    
                        # compute needed values:  N(k) = number of times state k is visited in chain, and sum(y_k) = sum of y's in state k 
                        if (m==1) {sumNcount<-formu(states0[[j]][-(n+1)],Y,K)
                        } else {sumNcount<-formu(Z[[j]][m-1,-(n+1)], Y,K)}   
                           sumtot<-cbind(sumNcount$sumy, sumNcount$ny)                                                                                                    #       sums<-sumNcount$sumy ;    tots<-sumNcount$ny
                      # new means          
                      mudraw<-apply(sumtot, 1,  function (x)     rnorm(1, mean= ((mu0/var0)+(x[1]/varknown)) / ((1/var0)+(x[2]/varknown)), sd= sqrt( 1/( (1/var0) + (x[2]/varknown)))  ))
                      MU[[j]][m,]<-mudraw  
                        
                    # 2 Update States given parameters
                  newZ<- UpdateStates( Y, Q[[j]][m,], MU[[j]][m,], initS= q0[[j]][m,], m)
                  Z[[j]][m,]<-newZ$Z
                  if(j==J) MAP[m]<-newZ$MAP


             } #end of PT chain loop
             
## PT move

          if(m>1) {TrackParallelTemp[m,]<-TrackParallelTemp[m-1,]}     
          if(m>20){
        #  if( sample(c(1,0),1, 0.9)==1){   # FREQ OF TEMPERING! 
          
## NEW, try each chain each time
    #Chain1<-sample( 1:(J-1), 1)   
    #Chain2<-Chain1+1

       
       if (sample(c(1,0),1, prob=c(0.9,.1))==1){
      if( m%%2==0){chainset<- c(1:(J-1))[c(1:(J-1))%%2==0]   #evens
      } else {chainset<- c(1:(J-1))[c(1:(J-1))%%2!=0] }   #odds

 for( eachChain in 1:length(chainset)){
                Chain1<-chainset[eachChain]  
                Chain2<-Chain1+1
                 PTsuccess[Chain1,"Tries"]<-     PTsuccess[Chain1,"Tries"]+1
                 PTsuccess[Chain2,"Tries"]<-     PTsuccess[Chain2,"Tries"]+1
                 
              Alpha1<-STORE_Alphas[[Chain1]]
              Alpha2<-STORE_Alphas[[Chain2]]

              Qchain1<-matrix(Q[[Chain1]][m,], nrow=K, byrow=TRUE)
              Qchain2<- matrix(Q[[Chain2]][m,]   , nrow=K, byrow=TRUE)
           MHratio<- parallelAcceptHMM_QV(Qchain1, Qchain2, Alpha1 ,Alpha2 , Z[[Chain1]][m,1], Z[[Chain2]][m,1])

# ratio if ergodic only:    # MHratio<- parallelAcceptHMM_ergodic(Qchain1, Qchain2, Alpha1 ,Alpha2)
            # MHratio<- parallelAcceptHMM_QV(Qchain1, Qchain2, Alpha1 ,Alpha2, Z[[Chain1]][m,1], Z[[Chain2]][m,1])
               
          if (MHratio==1){                                 # switch 
          PTsuccess[Chain1,"Success"]<-PTsuccess[Chain1, "Success"]+1
          PTsuccess[Chain2,"Success"]<-PTsuccess[Chain2, "Success"]+1
     
                   #new
                   .tpt1<-  TrackParallelTemp[m,Chain1 ]
                   .tpt2<-  TrackParallelTemp[m,Chain2 ]             
                  TrackParallelTemp[m,Chain1 ]<-.tpt2
                  TrackParallelTemp[m,Chain2 ]<-.tpt1 
                                                                                                         
          .p1<- q0[[Chain1]][m,]
          .p2<- q0[[Chain2]][m,]
          q0[[Chain1]][m,]<-.p2
          q0[[Chain2]][m,]<-.p1
                                                                                                             
          .m1<- MU[[Chain1]][m,]
          .m2<- MU[[Chain2]][m,]
          MU[[Chain1]][m,]<-.m2
          MU[[Chain2]][m,]<-.m1
                                                                                                               
          .s1<- Q[[Chain1]][m,]
          .s2<- Q[[Chain2]][m,]
          Q[[Chain1]][m,]<-.s2
          Q[[Chain2]][m,]<-.s1
          
          #ici
.Qold1<- Qold[[Chain1]]
.Qold2<- Qold[[Chain2]]
          Qold[[Chain1]]<-.Qold2
          Qold[[Chain2]]<-.Qold1
                                                                                                                # Zs
          .z1<- Z[[Chain1]][m,]
          .z2<- Z[[Chain2]][m,]
          Z[[Chain1]][m,]<-.z2
          Z[[Chain2]][m,]<-.z1
          }   }  


         PTsuccess[,"Ratio"]<-apply(PTsuccess[,c(2,3)], 1, function(x) x[2]/x[1])
              }}
          
   SteadyScore$K0[m]<-sum(table(Z[[J]][m,])>0)

# for all chains 
       K0Final[ m, ]<-sapply(   Z  ,  function(x)  sum(table(x[m,])>0))
#
#
#
# NEW
## Compute density f2 at this iteration (L1 norm)
f2now[m]<- density_f2(y=  Y,  .q0=q0[[J]][m ,] , .Q=Q[[J]][m ,], .mu=MU[[J]][m ,] )
fDist[m]<-   abs( f2now[m]-densTrue)

            }  # end of iteration loop
close(pb)
             PTsuccess[,"Ratio"]<-apply(PTsuccess[,c(2,3)], 1, function(x) x[2]/x[1])

allResults<-list("Means"=MU[[J]], "Trans"=Q[[J]], "States"=Z[[J]], "q0"=q0[[J]], "YZ"=YZ, "MAP"=MAP, "K0"=SteadyScore$K0, "PTsuccess"=PTsuccess, "f2dens"=f2now, "f2Dist"=fDist)

			     if(m ==M){Sys.sleep(0.01)
			      #  pdf( file=paste("HmmTracker_",lab, '.pdf', sep="") , height=4, width=12)
			          png( file=paste("HmmTracker_",lab, '.png', sep="") , height=800, width=1200)
			          par(mfrow=c(2,3))
			       #1
			       #   plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l', ylab="K0")
                              slices <- prop.table(table((factor(SteadyScore$K0[-(1:M/5)], levels=c(1:K)))))
                              barplot(slices, ylim=c(0,1), main="Number of non-empty states", xlab="Number of non-empty states", ylab="Probability (from MCMC)") 
                              abline(h=seq(0, 1, .05), lwd=0.5, col='LightGrey')


			      #2
			          MuXqoPlot( allResults, M/5, minq=0.01, plotlab=lab)  
			       #3
			          ts.plot(TrackParallelTemp, main='Track Parallel Tempering', col=rainbow(J), gpars=list(yaxt="n") , ylab="alpha")
			         # axis(2, at=1:J, tick=1:J, labels=round(AllAlphas[,2],4), las=2) 
			          axis(2, at=1:J, tick=1:J, labels=round(c( Alpha_lows),4), las=2) 
			          #ts.plot(Bigmu[[nCh]], main='emptying Mu', col=rainbow(k))
			       #4
			          image(Z[[J]][,order(Y)], col=rainbow(K), main="Allocations vs ordered Y")
			          #image(ZSaved[[nCh]][order(Y),], col=rainbow(K), main="Allocations")
			       #5  
			        ts.plot(fDist, main='L1 norm distance (y1, y2)')
			       #6 
			        plotbyK0<-cbind(  y= fDist[-(1:M/5)], x=SteadyScore$K0[-(1:M/5)])
			        boxplot( y~x, data=plotbyK0, main='L1 norm by K0')
			          dev.off()
			          Sys.sleep(0)}

      return(allResults)
      }


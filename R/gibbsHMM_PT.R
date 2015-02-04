 #' gibbsHMM_PT
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))



gibbsHMM_PT<-function(YZ, M=2000, K=10, mu0=0, var0=100,alphaMAX=1, alphaMin=1e-05, J=20, lab="sim", type="First"){
    #____SET UP_________________________________________
    ifelse(class(YZ)=='data.frame',    Y<-YZ$O, Y<-YZ)
    n=length(Y) # sample size
    varknown<-1 # known variace 
    
     # INITIALIZE
    startVal<-makeStart(Y, K);  states0<-startVal$states0   #FUNK
     
     # J=20
        TrackParallelTemp<-matrix(nrow=M, ncol=J)

         TrackParallelTemp[1,]<-c(1:J)
      # TO BE INCORPORATED INTO J LISTS
    MU<-replicate(J,  matrix(nrow=M, ncol=K),  simplify=F)
    Q  <-replicate(J, matrix(nrow=M, ncol=K*K) , simplify=F)
    q0 <-replicate(J, matrix(nrow=M, ncol=K), simplify=F)
    Z  <-replicate(J, matrix(nrow=M, ncol=n+1)  ,  simplify=F)  #include 0 for initial state to be estimated too?
   SteadyScore<-data.frame("Iteration"=c(1:M), "K0"=K) ##### THIS IS NEW##

    K0Final<-matrix(nrow=M, ncol=J)
    MAP<-c(1:M)  # KEEP TARGET ONLY
    
    # ALPHA
    
    #alphaMAX<-(K-1)*(1+K-2+alphaMin)*(1+1/( (1/2) - alphaMin*(K-1))) -(K-1)*alphaMin+0.1
    Alpha_lows<-c(alphaMAX, exp(seq(log(alphaMAX), log(alphaMin), length=J))[-1])
         
      #names(TrackParallelTemp)<-   AllAlphas[,2]
    # functions
    for (m in 1:M){ 
              
          if(m %% 20==0){Sys.sleep(0.01)
          par(mfrow=c(1,4))
          plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l')
          ts.plot(q0[[J]], main='q0 from target posterior', col=rainbow(K))
          ts.plot(TrackParallelTemp, main='Track Parallel Tempering', col=rainbow(J), gpars=list(yaxt="n") )
        #  axis(2, at=1:J, tick=1:J, labels=round(AllAlphas[,2],4), las=2) 
          axis(2, at=1:J, tick=1:J, labels=round(c( Alpha_lows),4), las=2)
          image(Z[[J]][,order(Y)], col=rainbow(K), main="Allocations vs ordered Y")
          Sys.sleep(0)}
          
                        
      for (j in 1:J){ # FOR EACH CHAIN           
          # FOR EACH CHAIN...
             # make matrix of alphas
                  AllAlphas<-matrix(Alpha_lows[j],ncol=K, nrow=K)
                  if(type=="First"){ 
                      AllAlphas[,1]<-alphaMAX      
                  } else if(type=="Diagonal") {
                      diag(AllAlphas)<-alphaMAX
                  }
                                # 1 Parameters given states Z(m-1)
                      # 1.1  Transition matrix Q from conditional posterior
                      if (m==1) {nt<-CountTrans(states0, K)
                        } else { nt<-CountTrans(Z[[j]][m-1,],K)}   # HERE ACCESS STATES
                          
                        # draw transition probs for state 1:K
                        # for (i in 1:K) qnew[i,]<-rdirichlet(par=  nt[i,]+AllAlphas[j, ])
                     #    for (i in 1:K) print(rdirichlet(par=  nt[i,]+AllAlphas[j, ]))

                   # qnew<- t(apply( nt,1, function(x) rdirichlet( par=x+AllAlphas[j, ])))
                    qnew<-matrix(ncol=K, nrow=K)
                    for(k in 1:K){qnew[k,]<-rdirichlet(par=nt[k,]+AllAlphas[k,])}


if(dim(getq0NEW(qnew))[1]>1){
        q0new <-  ALTERNATEq0(qnew)   
} else {                    
        q0new<-getq0NEW(qnew)  
}
if(sum(q0new<0)>0){ 
        q0new <-  ALTERNATEq0(qnew)   
}
  
                
                        #METROPOLIS Hastings STEP     
                  #      if (m==1){ A<-q0new[states0[1]]/startVal$q0[states0[1]]
                   #           ifelse(A>runif(1,c(0,0.99)), qok<-qnew , qok<-qnew)        
                   #     } else  { A<-q0new[Z[[j]][m-1,1]]/q0[[j]][m-1,Z[[j]][m-1,1]]
                    #                 U<-runif(1,c(0,0.99))
                         #         if (A>U){ qok<-qnew}
                if (m>1){   
                  A<-q0new[Z[[j]][m-1,1]]/q0[[j]][m-1,Z[[j]][m-1,1]]   
                    U<-runif(1,c(0,0.99))
                  if(A>runif(1,c(0,0.99))){ 
                     #    Accept new values
                        Q[[j]][m,]<-as.vector(t(qnew))
                        q0[[j]][m,]<-q0new
                                } else {
                       Q[[j]][m,]<-as.vector(t(Q[[j]][m-1,]))
                       q0[[j]][m,]<-q0[[j]][m-1,]  
                                }}else{ Q[[j]][m,]<-as.vector(t(qnew))
                        q0[[j]][m,]<-q0new  }

                        # 1.2 Update mu's    
                        # compute needed values:  N(k) = number of times state k is visited in chain, and sum(y_k) = sum of y's in state k 
                        if (m==1) {sumNcount<-formu(states0[-(n+1)],Y,K)
                        } else {sumNcount<-formu(Z[[j]][m-1,-(n+1)], Y,K)}   
                           sumtot<-cbind(sumNcount$sumy, sumNcount$ny)                                                                                                    #       sums<-sumNcount$sumy ;    tots<-sumNcount$ny
                      # new means          
                      mudraw<-apply(sumtot, 1,  function (x)     rnorm(1, mean= ((mu0/var0)+(x[1]/varknown)) / ((1/var0)+(x[2]/varknown)), sd= sqrt( 1/( (1/var0) + (x[2]/varknown)))  ))

                  #SAVE sampled parameters
                  #Q[[j]][m,]<-as.vector(t(qok))
                 # q0[[j]][m,]<-getq0NEW(qok)  
                  MU[[j]][m,]<-mudraw  
                        
                    # 2 Update States given parameters
                  newZ<- UpdateStates( Y, Q[[j]][m,], MU[[j]], initS= q0[[j]][m,], m)
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

       
       if (sample(c(1,0),1, prob=c(0.4,.6))==1){
      if( m%%2==0){chainset<- c(1:(J-1))[c(1:(J-1))%%2==0]   #evens
      } else {chainset<- c(1:(J-1))[c(1:(J-1))%%2!=0] }   #odds

 for( eachChain in 1:length(chainset)){
          Chain1<-chainset[eachChain]  
          Chain2<-Chain1+1
        # DOUBLE CHECK THEORY HERE
          # HOW TO DO DIAG PROPERLY??
          a1<-c(Alpha_lows[1],rep(Alpha_lows[Chain1], K-1))
          a2<-c(Alpha_lows[1],rep(Alpha_lows[Chain2], K-1))
         # MHratio<- parallelAccept(q0[[Chain1]][m,], q0[[Chain2]][m,], AllAlphas[ Chain1,] , AllAlphas[Chain2,] )

          MHratio<- parallelAccept(q0[[Chain1]][m,], q0[[Chain2]][m,], a1 ,a2 )
          if (MHratio==1){                                 # switch 
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
          
                                                                                                                # Zs
          .z1<- Z[[Chain1]][m,]
          .z2<- Z[[Chain2]][m,]
          Z[[Chain1]][m,]<-.z2
          Z[[Chain2]][m,]<-.z1
          }   }   }}
          
   SteadyScore$K0[m]<-sum(table(Z[[J]][m,])>0)

# for all chains 
       K0Final[ m, ]<-sapply(   Z  ,  function(x)  sum(table(x[m,])>0))

            }  # end of iteration loop
            if(m ==M){Sys.sleep(0.01)
         pdf( file=paste("HmmTracker_",lab, '.pdf', sep="") , height=4, width=12)
          par(mfrow=c(1,4))
          plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l', ylab="K0")
          ts.plot(q0[[J]], main='q0 from target posterior', col=rainbow(K))
          ts.plot(TrackParallelTemp, main='Track Parallel Tempering', col=rainbow(J), gpars=list(yaxt="n") , ylab="alpha")
         # axis(2, at=1:J, tick=1:J, labels=round(AllAlphas[,2],4), las=2) 
          axis(2, at=1:J, tick=1:J, labels=round(c( Alpha_lows),4), las=2)
          
          #ts.plot(Bigmu[[nCh]], main='emptying Mu', col=rainbow(k))
          image(Z[[J]][,order(Y)], col=rainbow(K), main="Allocations vs ordered Y")
          #image(ZSaved[[nCh]][order(Y),], col=rainbow(K), main="Allocations")
          dev.off()
          Sys.sleep(0)}
          
      return(list("Means"=MU[[J]], "Trans"=Q[[J]], "States"=Z[[J]], "q0"=q0[[J]], "YZ"=YZ, "MAP"=MAP, "K0"=SteadyScore$K0))
      }


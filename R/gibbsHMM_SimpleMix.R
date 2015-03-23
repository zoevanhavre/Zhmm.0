 #' gibbsHMM_PT
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))



gibbsHMM_SimpleMix<-function(YZ, M=2000, K=10,alphaMAX=1,alphaMin=0.001 ,PrbDiag=c("half", "fair"), lab="sim"){
    #____SET UP_________________________________________
    mu0=0; var0=100
    ifelse(class(YZ)=='data.frame',    Y<-YZ$Observed, Y<-YZ)
    n=length(Y) # sample size
    varknown<-1 # known variace 
    
     # INITIALIZE
    states0<-list()
    states0<-makeStartSimpler(Y, K)


    #states0<-makeStart(Y, K)
    #;  states0<-startVal$states0   #FUNK
     
    MU<-matrix(nrow=M, ncol=K)
    Q  <- matrix(nrow=M, ncol=K*K) 
          Qold<-  diag(K)
    q0 <-matrix(nrow=M, ncol=K)
    Z  <- matrix(nrow=M, ncol=n+1)    #include 0 for initial state to be estimated too?
   SteadyScore<-data.frame("Iteration"=c(1:M), "K0"=K) ##### THIS IS NEW##

    MAP<-c(1:M)  # KEEP TARGET ONLY
    
    # ALPHA
    
    #alphaMAX<-(K-1)*(1+K-2+alphaMin)*(1+1/( (1/2) - alphaMin*(K-1))) -(K-1)*alphaMin+0.1
    #Alpha_lows<-c(alphaMAX, exp(seq(log(alphaMAX), log(alphaMin), length=J))[-1])
    pb <- txtProgressBar(min = 0, max = M, style = 3)
         
      #names(TrackParallelTemp)<-   AllAlphas[,2]
    # functions
    for (m in 1:M){ 
          
          if(m %% 100==0){Sys.sleep(0.1)
          setTxtProgressBar(pb, m)
                  if(M < 20001){    
                    par(mfrow=c(1,3))
          plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l')
          ts.plot(q0, main='q0 from target posterior', col=rainbow(K))
        #  axis(2, at=1:J, tick=1:J, labels=round(AllAlphas[,2],4), las=2) 
          image(Z[,order(Y)], col=rainbow(K), main="Allocations vs ordered Y")
          Sys.sleep(0)}}
          

             # make matrix of alphas
AllAlphas<-matrix(alphaMin,ncol=K, nrow=K)
if (PrbDiag=="half"){
    if(sample( c(1, 0), size=1, prob=c(0.5, 0.5))==1){   # Put on diagonal
      diag(AllAlphas)<-alphaMAX 
    } else { 
ChooseColumn<-sample( c(1:K), size=1, prob=rep(1/K, K))    # draw non-diag position
AllAlphas[,ChooseColumn]<-alphaMAX      # make said column Amax
}
}else if (PrbDiag=="fair"){
  if(sample( c(1, 0), size=1, prob=c(1/K,(K-1)/K) )==1 ){     # Put on diagonal
      diag(AllAlphas)<-alphaMAX 
   } else { 
ChooseColumn<-sample( c(1:K), size=1, prob=rep(1/K, K))    # draw non-diag position
AllAlphas[,ChooseColumn]<-alphaMAX      # make said column Amax}
}   } 
  


                                # 1 Parameters given states Z(m-1)
                      # 1.1  Transition matrix Q from conditional posterior
                      if (m==1) {nt<-CountTrans(states0, K)
                        } else { nt<-CountTrans(Z[m-1,],K)}   # HERE ACCESS STATES
                          
                        # draw transition probs for state 1:K
                        # for (i in 1:K) qnew[i,]<-rdirichlet(par=  nt[i,]+AllAlphas[j, ])
                     #    for (i in 1:K) print(rdirichlet(par=  nt[i,]+AllAlphas[j, ]))

                   # qnew<- t(apply( nt,1, function(x) rdirichlet( par=x+AllAlphas[j, ])))
                    qnew<-matrix(ncol=K, nrow=K)
                    for(k in 1:K){qnew[k,]<-rdirichlet(par=nt[k,]+AllAlphas[k,])}

    q0new <-  ALTERNATEq0(qnew)  
                        if (m>1){   
                          q0Previous<-ALTERNATEq0(Qold)
                          A<-q0new[Z[m-1,1]]/q0Previous[Z[m-1,1]]   
                          U<-runif(1,c(0,0.99))
                  if(A>runif(1,c(0,0.99))){ 
                     #Accept new values
                        Q[m,]<-as.vector(t(qnew))
                        q0[m,]<-q0new
                                } else {
                       Q[m,]<-as.vector(t(Qold))

                     #  q0[[j]][m,]<-q0[[j]][m-1,]  
                         q0[m,]<-q0Previous
                                }}else{ Q[m,]<-as.vector(t(qnew))
                        q0[m,]<-q0new  }
Qold<-  matrix( Q[m,]  , K,K, byrow=TRUE)                     # **NEW**

                        # 1.2 Update mu's    
                        # compute needed values:  N(k) = number of times state k is visited in chain, and sum(y_k) = sum of y's in state k 
                         if (m==1) {sumNcount<-formu(states0[-(n+1)],Y,K)
                        } else {sumNcount<-formu(Z[m-1,-(n+1)], Y,K)}  
                           sumtot<-cbind(sumNcount$sumy, sumNcount$ny)                                                                                                    #       sums<-sumNcount$sumy ;    tots<-sumNcount$ny
                      # new means          
                      mudraw<-apply(sumtot, 1,  function (x)     rnorm(1, mean= ((mu0/var0)+(x[1]/varknown)) / ((1/var0)+(x[2]/varknown)), sd= sqrt( 1/( (1/var0) + (x[2]/varknown)))  ))
                      MU[m,]<-mudraw  
                        
                    # 2 Update States given parameters
                  newZ<- UpdateStates( Y, Q[m,], MU[m,], initS= q0[m,], m)
                  Z[m,]<-newZ$Z
                MAP[m]<-newZ$MAP


              #end of PT chain loop
             

          
   SteadyScore$K0[m]<-sum(table(Z[m,])>0)



            }  # end of iteration loop
close(pb)
            
            if(m==M){Sys.sleep(0.01)
      #   pdf( file=paste("HmmTracker_",lab, '.pdf', sep="") , height=4, width=12)
          png( file=paste("HmmTracker_",lab, '.png', sep="") , height=400, width=1200)
          par(mfrow=c(1,3))
          plot(SteadyScore$K0~SteadyScore$Iteration, main='#non-empty groups', type='l', ylab="K0")
          ts.plot(q0, main='q0 from target posterior', col=rainbow(K))
         # axis(2, at=1:J, tick=1:J, labels=round(AllAlphas[,2],4), las=2) 
          #ts.plot(Bigmu[[nCh]], main='emptying Mu', col=rainbow(k))
          image(Z[,order(Y)], col=rainbow(K), main="Allocations vs ordered Y")
          #image(ZSaved[[nCh]][order(Y),], col=rainbow(K), main="Allocations")
          dev.off()
          Sys.sleep(0)}
          
      return(list("Means"=MU, "Trans"=Q, "States"=Z, "q0"=q0, "YZ"=YZ, "MAP"=MAP, "K0"=SteadyScore$K0))
      }


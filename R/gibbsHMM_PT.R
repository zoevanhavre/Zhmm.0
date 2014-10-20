 #' gibbsHMM_PT
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))



gibbsHMM_PT<-function(YZ, M=2000, K=5, mu0=0, var0=100, alphaMin=0.5,  p=1){
    #____SET UP_________________________________________
    if (class(YZ)=='matrix')   Y<-YZ$O
    n=length(Y) # sample size


    MU<-matrix(nrow=M, ncol=K); colnames(MU)<-paste("mu", c(1:K), sep="")
    Q<-matrix(nrow=M, ncol=K*K) 
    q0<-matrix(nrow=M, ncol=K)
    Z<-matrix(nrow=M, ncol=n+1)  #include 0 for initial state to be estimated too?
  
    MAP<-c(1:M) 
    varknown<-1 # known variace 
     #ALPHA COMPUTATION
    alphaMAX<-(K-1)*(1+K-2+alphaMin)*(1+1/( (1/2) - alphaMin*(K-1))) -(K-1)*alphaMin+0.1
    alphas<-c(rep(alphaMAX,p), rep(alphaMin, K-p))    
    
   # INITIALIZE
    startVal<-makeStart(Y, K);  states0<-startVal$states0   #FUNK

    # functions
    for (m in 1:M){ 
         # Sys.sleep(0.1);   setTxtProgressBar(pb, m)

        # 1 Parameters given states Z(m-1)
        # 1.1  Transition matrix Q from conditional posterior
          if (m==1) {nt<-CountTrans(states0, K)
          } else { nt<-CountTrans(Z[m-1,],K)}   # REMEMBER you specified [m-1] from state storage matrix Z
            
          # draw transition probs for state 1:K
          qnew<-matrix(nrow=K, ncol=K)
          for (i in 1:K) qnew[i,]<-rdirichlet(par=  nt[i,]+alphas)
          q0new<-getq0(qnew)  

          #METROPOLIS Hastings STEP     
          if (m==1){ 
              A<-q0new[states0[1]]/startVal$q0[states0[1]]
              U<-runif(1,c(0,0.99))
                ifelse(A>U, qok<-qnew , qok<-qnew)        
          } else  { 
              A<-q0new[Z[m-1,1]]/q0[m-1,Z[m-1,1]]
              U<-runif(1,c(0,0.99))
                if (A>U){
                qok<-qnew
                }}
          

          #save transitions at iteration (m) in Q
          Q[m,]<-as.vector(t(qok))
          #compute new initial dist as ergodic distribution: and store
          q0[m,]<-getq0(qok)    ## recompute if met hastings step taken

        # 1.2 Update mu's    
          # compute needed values:  N(k) = number of times state k is visited in chain, and sum(y_k) = sum of y's in state k 
          if (m==1) {sumNcount<-formu(states0[-(n+1)],Y,K)
          } else {sumNcount<-formu(Z[m-1,-(n+1)], Y,K)}   
          sums<-sumNcount$sumy ;    tots<-sumNcount$ny

              
          # new means          
          mudraw<-matrix(rep(0,K))
          for(i in 1:K)  mudraw[i]<-rnorm(1, mean= ((mu0/var0)+(sums[i]/varknown)) / ((1/var0)+(tots[i]/varknown)), sd= sqrt( 1/( (1/var0) + (tots[i]/varknown)))  )
          MU[m,]<-mudraw  
            
        # 2 Update States given parameters
          ## set up:
            T<-n+1  # since we want 0-T states estimated
            trans<-qok
            initS<-q0[m,]
            OneStep<-matrix(nrow=T, ncol=K)
            const<-rep(0,n)
            FilterP<-matrix(nrow=T-1, ncol=K)
            Smooth<-matrix(nrow=T, ncol=K)

          # 2.1 Run Forward Filter to get P(St=j| yt, .) for j=1,...K and ti=1,...,T obs  (calling time var ti as t is used)
          #     Obtain P(Z(t)=l|Y(t-1),.) 1 step ahead pred of Z(ti)
            #for t=1 
              #compute P(Z(1)=l|t(0),.) using initial distribution of states
              #P(Z(1)=1|Y(0), .)
              for (i in 1:K) OneStep[1,i]<- sum(initS[1:K] * trans[i,])
              const[1]<-sum(dnorm(Y[1],mean=MU[1,1:K])*OneStep[1,])
              FilterP[1,]<- (dnorm(Y[1],mean=MU[1,])*OneStep[1,]) /const[1]   #Obtain filtered probs P(Z(1)|Y(1))
                  
                
            #for ti=2:100,
            for (ti in 2:n){
              # for states 
              for (i in 1:K) OneStep[ti,i]<-  sum(trans[i,]*FilterP[ti-1, ])  
              #Filter:  P(Z(t)|Y(t))
              const[ti]<-sum( dnorm(Y[ti],mean=MU[m,])*OneStep[ti,]) 
              FilterP[ti,]<- (dnorm(Y[ti],mean=MU[m,])*OneStep[ti,]) /const[ti]   
              }

            MAP[m]<-sum(log(const))

            # 2.2 Sample final state Z(T) from P(Z(T)=j| y(T),.) (filter from t=100)
              Smooth[T,]<-FilterP[T-1,]
              Z[m,T]<-sample(c(1:K), size=1, prob=Smooth[T,])

            # 2.3 Sample Z(t) for t=T-1, T-2,..., 
            for (ti in (T-1):1){    # smooth
              Smooth[ti,]<- (trans[Z[m,ti+1],]*FilterP[ti, ]) / sum(trans[Z[m,ti+1],]*FilterP[ti,])
              Z[m,ti]<-sample(c(1:K),size=1, prob=Smooth[ti,])  }       # Draw state & STORE Z(t)'s in a T*1 matrix
              }
            #output
      return(list("Means"=MU, "Trans"=Q, "States"=Z, "q0"=q0, "YZ"=YZ, "Smooth"=Smooth, "MAP"=MAP))
      }


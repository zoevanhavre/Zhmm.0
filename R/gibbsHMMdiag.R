 #' gibbsHMMdiag
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))



gibbsHMMdiag<-function(YZ, M=2000, K=5, mu0=0, var0=100, alphaMin=0.5,   p=1){
  Y<-YZ$Obs
    #____SET UP_________________________________________
    n=length(Y) # n = 100 = T,  sample size
    MU<-matrix(nrow=M, ncol=K); colnames(MU)<-paste("mu", c(1:K), sep="")
    Q<-matrix(nrow=M, ncol=K*K) 
    q0<-matrix(nrow=M, ncol=K)
    S<-matrix(nrow=M, ncol=n+1)  #include 0 for initial state to be estimated too?
    MAP<-c(1:M)
     makeStart<-function(Y,k=K){
        n<-length(Y)
        states0<-rep(0,n)

        # split data by quantiles &  compute allocations
        states0<-as.numeric(cut2(Y, m=5, g=k))
        #, combine with data
        YX0<-cbind(Y, states0)
        # compute means 
        means0<-gapply(as.data.frame(YX0), FUN=function(x) mean(x$Y), groups=YX0[,2])
        # compute transition matrix
        
        trans0 <- matrix(nrow = k, ncol = k, 0)
        for (t in 1:(length(states0) - 1)) trans0[states0[t], states0[t + 1]] <- trans0[states0[t], states0[t + 1]] + 1
        for (i in 1:k) trans0[i, ] <- trans0[i, ] / sum(trans0[i, ])

        #stationary dist q0 for this:
        q0<-getq0(trans0)

        #initial unobserved state from stationary dist
        initstates0<-sample(c(1:k), size=1, prob=q0)
        states0<-c(initstates0, states0)
        return(list(means0=means0, states0=states0, trans0=trans0, q0=q0))  }

    startVal<-makeStart(Y, K);  states0<-startVal$states0 
    varknown<-1


    #    alphaMIN<-1/(4*(K-1))
    alphaMAX<-(K-1)*(1+K-2+alphaMin)*(1+1/( (1/2) - alphaMin*(K-1))) -(K-1)*alphaMin+0.1
    alphas<-diag(K)*alphaMAX
    alphas<-apply(alphas,c(1,2), function(x) ifelse(x==0, x<-alphaMin, x<-x) )

      # functions
    CountTrans<-function(stateChain,K=K){
      nt <- matrix(nrow = K, ncol = K, 0)
      for (t in 1:(length(stateChain) - 1)) nt[stateChain[t], stateChain[t + 1]] <- nt[stateChain[t], stateChain[t + 1]] + 1
      return(nt)}
    formu<-function(stateChain, Y, K=K){
      stateChain <-factor(stateChain, levels=c(1:K))
      # number in each state and sum of Ys
      xy<-as.data.frame(cbind(Y, stateChain))
      ny<-gapply(xy, FUN=function(x) length(x$Y), groups=stateChain)
      sumy<-gapply(xy, FUN=function(x) sum(x$Y), groups=stateChain)
      return(list(sumy=sumy, ny=ny))}  
    nk<-function(States, K=K){
      nkOrd<-matrix(nrow=dim(States)[1], ncol=K)
      for (i in 1:dim(States)[1]){
      rowS<-factor(States[i,], levels=1:K)
      nkOrd[i,]<-table(rowS)
      }
      return(nkOrd)
      }
    getq0<-function(Q){
      K<-dim(Q)[1]
      U<-matrix(rep(1/K, K*K), K,K)
      u<-rep(1/K, K)
      I<-diag(K)
      ipu<-I-Q+U
      u%*%solve(ipu)
      }      
    rdirichlet<-function(m=1,par){
          k=length(par);  mat=matrix(0,m,k)       
          for (i in 1:m)  {sim=rgamma(k,shape=par,scale=1); mat[i,]=sim/sum(sim)}
          mat             
          }
  
    ##################
    # START MAIN LOOP HERE for m in 1:M
          pb <- txtProgressBar(min = 0, max = M, style = 3)
        for (m in 1:M){ 
        Sys.sleep(0.1);   setTxtProgressBar(pb, m)

        # 1 Parameters given states S(m-1)
        # 1.1  Transition matrix Q from conditional posterior
          if (m==1) {nt<-CountTrans(states0, K)
          } else { nt<-CountTrans(S[m-1,],K)}   # REMEMBER you specified [m-1] from state storage matrix S
            
          # draw transition probs for state 1:K
          qnew<-matrix(nrow=K, ncol=K)
          #for (i in 1:K) qnew[i,]<-rdirichlet(par=  nt[i,]+alphas)
          for (i in 1:K) qnew[i,]<-rdirichlet(par=  nt[i,]+alphas[i,])


      #q0new<-getq0(qnew)  
               # HERE IS TRICK TO GET Q0...
      altq0<-function(qnew){ estq0<-round(qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew%*%qnew, 4)
      okid<-c(1:K)[apply((estq0 == 1), 1, sum)==0]  ;      return( estq0[okid[1], ])}
      .q0new <- tryCatch({ getq0(qnew)  }, error = function(e) { altq0(qnew)     })
      if(sum(.q0new<0)>0){q0new<-altq0(qnew)
   } else {q0new<-.q0new}


          #METROPOLIS Hastings STEP     
          if (m==1){ 
              A<-q0new[states0[1]]/startVal$q0[states0[1]]
              U<-runif(1,c(0,0.99))
                ifelse(A>U, qok<-qnew , qok<-qnew)        
          } else  { 
              A<-q0new[S[m-1,1]]/q0[m-1,S[m-1,1]]
              U<-runif(1,c(0,0.99))
                if (A>U){
                qok<-qnew
                }}
          

          #save transitions at iteration (m) in Q
          Q[m,]<-as.vector(t(qok))
          #compute new initial dist as ergodic distribution: and store
          #q0[m,]<-compdelta(qok)    ## recompute if met hastings step taken
          .q0 <- tryCatch({ getq0(qok)  }, error = function(e) { altq0(qok)     })
         if(sum(.q0<0)>0){ q0[m,]<- altq0(qok)
          }else{ q0[m,]<- .q0}

        # 1.2 Update mu's    
          # compute needed values:  N(k) = number of times state k is visited in chain, and sum(y_k) = sum of y's in state k 
          if (m==1) {sumNcount<-formu(states0[-(n+1)],Y,K)
          } else {sumNcount<-formu(S[m-1,-(n+1)], Y,K)}   
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
          #     Obtain P(S(t)=l|Y(t-1),.) 1 step ahead pred of S(ti)
            #for t=1 
              #compute P(S(1)=l|t(0),.) using initial distribution of states
              #P(S(1)=1|Y(0), .)
              for (i in 1:K) OneStep[1,i]<- sum(initS[1:K] * trans[i,])
              const[1]<-sum(dnorm(Y[1],mean=MU[1,1:K])*OneStep[1,])
              FilterP[1,]<- (dnorm(Y[1],mean=MU[1,])*OneStep[1,]) /const[1]   #Obtain filtered probs P(S(1)|Y(1))
                  
                
            #for ti=2:100,
            for (ti in 2:n){
              # for states 
              for (i in 1:K) OneStep[ti,i]<-  sum(trans[i,]*FilterP[ti-1, ])  
              #Filter:  P(S(t)|Y(t))
              const[ti]<-sum( dnorm(Y[ti],mean=MU[m,])*OneStep[ti,]) 
              FilterP[ti,]<- (dnorm(Y[ti],mean=MU[m,])*OneStep[ti,]) /const[ti]   
              }

            MAP[m]<-sum(log(const))

            # 2.2 Sample final state S(T) from P(S(T)=j| y(T),.) (filter from t=100)
              Smooth[T,]<-FilterP[T-1,]
              S[m,T]<-sample(c(1:K), size=1, prob=Smooth[T,])

            # 2.3 Sample S(t) for t=T-1, T-2,..., 
            for (ti in (T-1):1){    # smooth
              Smooth[ti,]<- (trans[S[m,ti+1],]*FilterP[ti, ]) / sum(trans[S[m,ti+1],]*FilterP[ti,])
              S[m,ti]<-sample(c(1:K),size=1, prob=Smooth[ti,])  }       # Draw state & STORE S(t)'s in a T*1 matrix
              }
            #output
      return(list("Means"=MU, "Trans"=Q, "States"=S, "q0"=q0, "YZ"=YZ, "Smooth"=Smooth, "MAP"=MAP))
      }
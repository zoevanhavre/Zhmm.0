 #' gibbsHMMfreeVar
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

gibbsHMMfreeVar<-function(YZ, M=2000, K=5, mu0=0, var0=100, alphaMin=0.5,  p=1){
    #____SET UP_________________________________________
    Y<-YZ$O
    
           a=2.5; b=2/var(Y); tau=1  #NEW
    SIGMA<-matrix(nrow=M, ncol=K); colnames(SIGMA)<-paste("Sigma", c(1:K), sep="") # NEW

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
            return(list(means0=means0, states0=states0, trans0=trans0, q0=q0))
          }
    startVal<-makeStart(Y, K);  states0<-startVal$states0 
    varknown<-1
    
    #FIX THIS
    #alphaMIN<-1/(4*(K-1))
    #alphaMAX<-(4*(K-1)*(K-1)) + 0.75 + A
    alphaMAX<-(K-1)*(1+K-2+alphaMin)*(1+1/( (1/2) - alphaMin*(K-1))) -(K-1)*alphaMin+0.1
    alphas<-c(rep(alphaMAX,p), rep(alphaMin, K-p))    
    

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

    forsigma<-function(stateChain, Y, K=K, muNow){
      stateChain <-factor(stateChain, levels=c(1:K))
      # number in each state and sum of Ys
      xy<-data.frame("Y"=Y, "stateChain"= stateChain, "muMatch"=0, "Ymu"=0)
      for (i in 1:length(stateChain)){  xy$muMatch[i]<-muNow[stateChain[i]]}
      xy$Ymu<-(xy$Y-xy$muMatch)^2
      sumy<-gapply(xy, FUN=function(x) sum(x$Ymu), groups=stateChain)
      return(sumy)}  


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
   
 
    for (m in 1:M){ 
         # Sys.sleep(0.1);   setTxtProgressBar(pb, m)

        # 1 Parameters given states S(m-1)
        # 1.1  Transition matrix Q from conditional posterior
          if (m==1) {nt<-CountTrans(states0, K)
          } else { nt<-CountTrans(S[m-1,],K)}   # REMEMBER you specified [m-1] from state storage matrix S
            
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
              A<-q0new[S[m-1,1]]/q0[m-1,S[m-1,1]]
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
                # if (m==1) {sumNcount<-formu(states0[-(n+1)],Y,K)
                # } else {sumNcount<-formu(S[m-1,-(n+1)], Y,K)}   

          if (m==1) {sumNcount<-formu(states0[-1],Y,K)
          } else {sumNcount<-formu(S[m-1,-1], Y,K)}   

          sums<-sumNcount$sumy ;    tots<-sumNcount$ny

              
          # new means          
          if(m==1){  MU[m,]<- rnorm(K, mean=(mu0*var0+sums)/(var0+tots),  sd=sqrt(var0/(var0+tots))) # must be sqrt as r takes in sd not var
          } else { MU[m,]<-    rnorm(K, mean=(mu0*var0+sums)/(var0+tots),  sd=sqrt(SIGMA[m-1,]/(var0+tots)))} # must be sqrt as r takes in sd not var
          

          if (m==1) {sv<-forsigma(states0[-1], Y, K , MU[m,]) 
          } else {sv<-forsigma(S[m-1,-1], Y, K , MU[m,]) }  
           # changes, added /ns
        
         # 7 Generate Sigma's (and save)
         SIGMA[m,]<-  rinvgamma(K, a+(sumNcount$ny+1)/2,  b+0.5*tau*(MU[m,]-mu0)^2+0.5*sv)


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
      return(list("Means"=MU, "Trans"=Q, "States"=S,"SIGMA"=SIGMA, "q0"=q0, "YZ"=YZ, "Smooth"=Smooth, "MAP"=MAP))
      }



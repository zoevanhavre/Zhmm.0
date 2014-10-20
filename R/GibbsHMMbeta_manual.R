#' GibbsHMMbeta_manual
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))
GibbsHMMbeta_manual<-function(Y, M,   alpha1=aMax, alpha2=alphaMin, beta1=alphaMin, beta2=aMax ){
      K<-2 
      

      #SETUP 
      n=length(Y) # n = 100 = T,  sample size
      K=2   # number of states  (Known)
        # functions:                
          makeStart<-function(Y){
            n<-length(Y) ; states0<-rep(0,n+1)
            # split data by quantiles &  compute allocations
              for (i in 1:n){ ifelse(Y[i]>quantile(Y,0.5),  states0[i+1]<-2, states0[i+1]<-1) } #lower values group 1, larger values group 2
              YX0<-cbind(Y, states0[-n])#, combine with data
              means0<-c( mean(subset(YX0, YX0[,2]==1)[,1]),mean(subset(YX0, YX0[,2]==2)[,1]))   # compute means     
              trans0<-CountTrans(states0) ; trans0[1,]<-trans0[1,]/sum(trans0[1,]); trans0[2,]<-trans0[2,]/sum(trans0[2,])  # compute transition matrix
              # stationary dist
              q1<-trans0[1,1]; q2<-trans0[2,2]; q0<-c((1-q1)/(1+(1-q2)-q1), (1-q1)/(1+(1-q2)-q1))
              states0[1]<-sample(c(1,2), size=1, prob=q0) #initial unobserved state from stationary dist
              return(list(means0=means0, states0=states0, trans0=trans0, q0=q0))
            }
          CountTrans<-function(stateChain){
                 nt11<-0;nt12<-0;nt21<-0;nt22<-0
                 for (i in 2:length(stateChain)){
                  if    (stateChain[i-1]==1 & stateChain[i]==1) { nt11<-nt11+1
                  } else if (stateChain[i-1]==1 & stateChain[i]==2) { nt12<-nt12+1
                  } else if (stateChain[i-1]==2 & stateChain[i]==1) { nt21<-nt21+1
                  } else if (stateChain[i-1]==2 & stateChain[i]==2) { nt22<-nt22+1 }}
                  return(rbind(c(nt11,nt12), c(nt21, nt22)))}
          formu<-function(stateChain, Y){
                #subset
                s1<-subset(cbind(stateChain, Y), stateChain==1)
                s2<-subset(cbind(stateChain, Y), stateChain==2)
                # time in each state
                sumy<-c(sum(s1[,2]), sum(s2[,2]))
                # sum of ys in each state
                ny<-c(dim(s1)[1], dim(s2)[1])
                return(c(sumy, ny))}
            
      startVal<-makeStart(Y)
      states0<-startVal$states0
      #empty storage matrices for pars to be estimated
        MU<-matrix(nrow=M, ncol=2); colnames(MU)<-c("mu1", "mu2")
        Q<-matrix(nrow=M, ncol=K) ;colnames(Q)<-c("q1", "q2")
        q0<-matrix(nrow=M, ncol=K)
        S<-matrix(nrow=M, ncol=n+1)  #include 0 for initial state to be estimated too

    
      ##################
      # START MAIN LOOP HERE for m in 1:M

      for (m in 1:M){
          # 1 Parameters given states S(m-1)
      #    Starting with some state process S(0)

       # 1.1  Transition matrix Q from conditional posterior
      if (m==1) {nt<-CountTrans(states0)
      } else { nt<-CountTrans(S[m-1,])}   # REMEMBER you specified [m-1] from state storage matrix S
        
      # draw transition probs for state k=1 and k=2
      q1new<-rbeta(1,alpha1+ nt[1,1] , alpha2+nt[1,2])
      q2new<-rbeta(1,beta1+ nt[2,2] , beta2+nt[2,1]) 

      if(q1new+q2new==2){ q0new<-c(1,0)
      }else{ 
        q0new<-c((1-q2new)/(1+(1-q2new)-q1new), (1-q1new)/(1+(1-q2new)-q1new))
      }

      # METROPOLIS Hastings STEP
      
      if (m==1){ 
         # A<-q0new[states0[1]]/startVal$q0[states0[1]]
         # U<-runif(1,c(0,1))
         # if (A>U){
              q1<-q1new; q2<-q2new 
           #  }   else {   q1<-  ; q2<- }       
      } else  { 
          A<-q0new[S[m-1,1]]/q0[m-1,S[m-1,1]]
          U<-runif(1,c(0,1))
            if (A>U){
            q1<-q1new; q2<-q2new 
            }  else {   q1<- Q[m-1,1]; q2<- Q[m-1,2] }}
      
      #save transitions at iteration (m) in Q
      Q[m,]<-c(q1,q2)  
      #compute new initial dist as ergodic distribution: and store
      if (q1+q2==2){q0[m,]<-c(1,0)
      } else {
      q0[m,]<-c((1-q2)/(1+(1-q2)-q1), (1-q1)/(1+(1-q2)-q1))  }
      
      # 1.2 Update mu's
            
          # compute needed values:  N(k) = number of times state k is visited in chain, and sum(y_k) = sum of y's in state k 
          if (m==1) {sumNcount<-formu(states0[-(n+1)],Y)
          } else {sumNcount<-formu(S[m-1,-(n+1)], Y)}   
          sums<-sumNcount[1:2] ;    tots<-sumNcount[3:4]

          # sample mu's from full conditionals:
          mu1<-rnorm(1, mean= sums[1]/(tots[1]+1), sd= sqrt( 1/(tots[1]+1)))
          mu2<-rnorm(1, mean= sums[2]/(tots[2]+1), sd= sqrt( 1/(tots[2]+1)))

          #save at iteration (m)
          MU[m,]<-c(mu1, mu2)
          
      # 2 States given parameters
          ## set up:
            T<-n+1  # since we want 0-T states estimated
            #need transition matrix at last iteration:
            q1<-Q[m,1]; q2<-Q[m,2]
            trans<-matrix(c(q1, 1-q1, 1-q2, q2), 2,2, byrow = TRUE)
            #need initial dist of states 
            initS<-q0[m,]
            #new storage for each iteration for specific probs:
            OneStep<-matrix(nrow=T, ncol=K) ;const<-rep(0,100)  ;  FilterP<-matrix(nrow=T-1, ncol=K);Smooth<-matrix(nrow=T, ncol=K)

            # 2.1 Run Forward Filter to get P(St=j| yt, .) for j=1,...K and ti=1,...,T obs  (calling time var ti as t is used)
            #     Obtain P(S(t)=l|Y(t-1),.) 1 step ahead pred of S(ti)
              #for t=1 
                #compute P(S(1)=l|t(0),.) using initial distribution of states
                #P(S(1)=1|Y(0), .)
                OneStep[1,1]<-  initS[1] * trans[1,1]+initS[2]*trans[2,1]
                OneStep[1,2]<-  initS[1] * trans[1,2]+initS[2]*trans[2,2]
                
              #  Obtain filtered probs P(S(1)|Y(1))
                const[1]<-dnorm(Y[1],mean=MU[1,1])*OneStep[1,1] + dnorm(Y[1],mean=MU[1,2])*OneStep[1,2]
                FilterP[1,1]<- (dnorm(Y[1],mean=MU[1,1])*OneStep[1,1]) /const[1]
                FilterP[1,2]<- (dnorm(Y[1],mean=MU[1,2])*OneStep[1,2])  /const[1]
                  
              #for ti=2:100,
              for (ti in 2:n){
              # for states l=1,2
                OneStep[ti, 1]<- trans[1,1]*FilterP[ti-1, 1] + trans[2,1]*FilterP[ti-1, 2]
                OneStep[ti, 2]<-trans[1,2]*FilterP[ti-1, 1] + trans[2,2]*FilterP[ti-1, 2]
                
              #Filter:  P(S(t)|Y(t))
                const[ti]<-dnorm(Y[ti],mean=MU[m,1])*OneStep[ti,1] + dnorm(Y[ti],mean=MU[m,2])*OneStep[ti,2]
                FilterP[ti,1]<- (dnorm(Y[ti],mean=MU[m,1])*OneStep[ti,1]) /const[ti]
                FilterP[ti,2]<- (dnorm(Y[ti],mean=MU[m,2])*OneStep[ti,2])  /const[ti]   }

            # 2.2 Sample final state S(T) from P(S(T)=j| y(T),.) (filter from t=100)
            Smooth[T,]<-FilterP[T-1,]
            S[m,T]<-sample(c(1,2), size=1, prob=Smooth[T,])

            # 2.3 Sample S(t) for t=T-1, T-2,..., 
            for (ti in (T-1):1){
            # smooth
            Smooth[ti,1]<- (trans[1,S[m,ti+1]]*FilterP[ti, 1]) / (trans[1,S[m,ti+1]]*FilterP[ti, 1]+trans[2,S[m,ti+1]]*FilterP[ti, 2])
            Smooth[ti,2]<- (trans[2,S[m,ti+1]]*FilterP[ti, 2]) / (trans[1,S[m,ti+1]]*FilterP[ti, 1]+trans[2,S[m,ti+1]]*FilterP[ti, 2])
            # Draw state & STORE S(t)'s in a T*1 matrix
            S[m,ti]<-sample(c(1,2),size=1, prob=Smooth[ti,])
            }
            rm(q1,q2)

      }

      #output
      return(list("Means"=MU, "Trans"=Q, "States"=S, "q0"=q0, "Y"=Y))
      }


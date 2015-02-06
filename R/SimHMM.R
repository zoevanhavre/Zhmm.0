 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

SimHMM<-function(Q,Mu, q0=NA, n=100){
        k<-dim(Q)[1]
        #stationary dist q0 for this:
      if(is.na(q0)){ q0<-getq0NEW(Q)}
        #state chain & ys
        X<-c(rep(0,n))
        Y<-c(rep(0,n))
        X[1]<-sample(c(1:k), size=1,prob =q0) 

        for (i in 2:n){ 
        #states
        X[i]<-sample(c(1:k), size=1,prob =Q[X[i-1],]) }

                                #y's
                        #  for (i in 1:n){          Y[i]<-rnorm(1, Mu[X[i]], 1)}
         Y<-sapply( X, function(x)   rnorm(1, Mu[x], 1) ) 

        return(data.frame("States"=X, "Observed"=Y))
        }

 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))



#FunkSim1<-function(n){  data.frame( "States"=1, "Observed"=rnorm(n, mean=3, sd=1)) }

 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))


FunkSim1<-function(n){   SimHMM( Q=matrix( c(  0.2,0.3,0.5,    0.5,0.25,0.25,    0.25, 0.65, 0.1), nrow=3, byrow=T),  Mu= c(0,3,10),q0=NA, n)  }

 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))


FunkSim2<-function(n){   SimHMM( Q=matrix( c(  0.2,0.3,0.5,    0.5,0.25,0.25,    0.25, 0.65, 0.1), nrow=3, byrow=T),  Mu= c(0,2,4),q0=NA, n)  }

#' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))


FunkSim3<-function(n){  
q1<-matrix(0.1, ncol=5, nrow=5)
diag(q1)<-.6
SimHMM( Q=q1,  Mu= c(-10, -5, 0, 5, 10),q0=NA, n)  }


#' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

FunkSim4<-function(n){
        nrep<-n/100
        X_100<-c( rep(1, 20), rep(2, 15), rep(1, 5), rep(3, 3), rep( 1, 25), rep(2, 10), rep(1, 20), rep(3, 2))
         X<-rep(X_100, nrep)
         Y<-sapply( X, function(x)   rnorm(1, c(1,4,-10)[x], 1) ) 
     
       return(data.frame("States"=X, "Observed"=Y))
        }



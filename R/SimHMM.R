 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

SimHMM<-function(Q,Mu, n=100){
        k<-dim(Q)[1]
        #stationary dist q0 for this:
        q0<-getq0(Q)

        #state chain & ys
        X<-c(rep(0,n))
        Y<-c(rep(0,n))
        X[1]<-sample(c(1:k), size=1,prob =q0) 

        for (i in 2:n){ 
        #states
        X[i]<-sample(c(1:k), size=1,prob =Q[X[i-1],]) }

        #y's
        for (i in 1:n){
          Y[i]<-rnorm(1, Mu[X[i]], 1)}

        return(data.frame("States"=X, "Observed"=Y))
        }

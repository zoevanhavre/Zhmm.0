 #' Function to compute stationary distribution of HMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))


SimHMM2K<-function(q1,q2,mu1,mu2, n=100){
    #trans matrix
    Q<- matrix(c(q1, 1-q1, 1-q2, q2), 2,2, byrow = TRUE)
    #stationary dist q0 for this:
    a<-q1;b<-1-q2
    q0<-c(b/(1+b-a), (1-a)/(1+b-a))

    #state chain & ys
    X<-c(rep(0,n))
    Y<-c(rep(0,n))
    X[1]<-sample(c(1,2), size=1,prob =q0) 

    for (i in 2:n){ 
    #state
    X[i]<-sample(c(1,2), size=1,prob =Q[X[i-1],]) 
    #y's
      if (X[i]==1) { Y[i]<-rnorm(1, mu1, 1)
      } else {Y[i]<-rnorm(1,mu2, 1) }}

    return(data.frame("States"=X, "Observed"=Y))
    }      
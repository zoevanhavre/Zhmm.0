#' gibbsHMM
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

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
#' unswStates
#'
#' This function draws samples from a Wishart dist
#' @param mydata, alpha, Krange
#' @keywords bayes factor, bayesm
#' @export
#' @examples
#' # 
#'  #

unswStates<-function( Out, K=K){
        # matrices to store 
            permutations <- function(n){
              if(n==1){
                  return(matrix(1))
              } else {
                  sp <- permutations(n-1)
                  p <- nrow(sp)
                  A <- matrix(nrow=n*p,ncol=n)
                  for(i in 1:n){
                      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
                  }
                  return(A)
              }}

          perm<-permutations(n =K, r = K, v = 1:K)      #compute all permutations
          nperm<-dim(perm)[1]
          niter<-length(Out$MAP)              
          SoldSnow<-matrix(nrow=2, ncol=nperm)      # states to compare 2xnS
          perdist<-rep(0,nperm)             # stored distances  1xnperm

          #to store final unswitched results
          finMu<-Out$Means
          finQ<-Out$Trans
          finq0<-Out$q0
          finS<-Out$States

          Zmap<-factor(Z[,which.max(Out$MAP)], levels=1:K, ordered=FALSE)


          for(i in 1:dim(Out$States)[1]){   #for each iteration
                Snow<-factor(Out$States[i,], levels=1:K, ordered=FALSE); Sold<-Zmap; Snew<-Snow
                PermScore<-function(x) {  Zw<-Snow;   levels(Zw)<-x; sum(as.numeric(as.character(Zw))!=Zmap)   }
                perdist<- apply(perm, 1,PermScore)       
     
            minPerm <-perm[which.min(perdist),]   
                levels(Snew)<-minPerm     ; finS[i,]<-as.numeric(as.character(Snew))  #format states
                permPARS<-c(1:K)[order(minPerm)]
                finMu[i,]<-Out$Means[i, permPARS]   
                finq0[i,]<-Out$q0[i, permPARS]

                newmat<-matrix(Out$Trans[i,], nrow=K, byrow=TRUE);  copymat<-matrix(Out$Trans[i,], nrow=K, byrow=TRUE)
                for(Row in 1:K) { for(Col in 1:K){ newmat[Row,Col]<-copymat[ permPARS[Row], permPARS[Col]]  }}
                finQ[i,]<-as.vector(t(newmat))  
                  }   
          
          # posterior Z (MAP)                  
          maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
          Zhat<- factor( apply(t(finS), 2,maxZ), levels=1:K)
          Mhat<-apply(finMu, 2, mean)
          Qhat<-apply(finQ, 2, mean)
          q0hat<-apply(finq0, 2, mean)


          return(list("Means"=finMu, "Trans"=finQ,"q0"=finq0, "States"=finS, "Y"=Out$Y, "hat_mu"=Mhat, "hat_Q"=Qhat, "hat_q0"=q0hat, "hat_Z"=Zhat))
            }


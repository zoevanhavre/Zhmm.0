 #' processGibbsFREEVAR
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

processGibbsFREEVAR<-function(Result, K=5,plotname="", simMeans=c(5,2,8)){    
        YZ<-Result$YZ
        filelabel<-plotname
        n<-length(YZ$Observed)

        TrimThin<-function(out, burn=200, newN=100){
                ids<- c(1:length(out$MAP))   #make dummy id var
                ids<-ids[-c(0:burn)]        # take out burn
                ids<-seq(from=min(ids), to=max(ids), by=round(length(ids)/newN))    # random sample or take every X id?
                #save selected iterations
                MU<-out$Means[ids,]
                SIGMA<-out$SIGMA[ids,]
                Q<-out$Trans[ids,]
                q0<-out$q0[ids,]
                S<-out$States[ids,]
                MAP<-out$MAP[ids]
                return(list("Means"=MU, "SIGMA"=SIGMA, "Trans"=Q,"q0"=q0, "States"=S, "Y"=out$Y, "MAP"=MAP))
                }
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

            perm<-permutations(K)      #compute all permutations
            nperm<-dim(perm)[1]
            niter<-length(Out$MAP)              
            SoldSnow<-matrix(nrow=2, ncol=nperm)      # states to compare 2xnS
            perdist<-rep(0,nperm)             # stored distances  1xnperm

            #to store final unswitched results
            finMu<-Out$Means;finQ<-Out$Trans;finq0<-Out$q0;finS<-Out$States; finSig<-Out$SIGMA
            Zmap<-factor(Out$States[which.max(Out$MAP),], levels=1:K, ordered=FALSE)

            for(i in 1:dim(Out$States)[1]){   #for each iteration
                  Snow<-factor(Out$States[i,], levels=1:K, ordered=FALSE); Sold<-Zmap; Snew<-Snow
                  PermScore<-function(x) {  Zw<-Snow;   levels(Zw)<-x; sum(as.numeric(as.character(Zw))!=Zmap)   }
                  perdist<- apply(perm, 1,PermScore)       
       
              minPerm <-perm[which.min(perdist),]   ;  levels(Snew)<-minPerm     ; finS[i,]<-as.numeric(as.character(Snew))  #format states
              permPARS<-c(1:K)[order(minPerm)] 
              finMu[i,]<-Out$Means[i, permPARS]
              finq0[i,]<-Out$q0[i, permPARS]
              finSig[i,]<-Out$SIGMA[i, permPARS]
              newmat<-matrix(Out$Trans[i,], nrow=K, byrow=TRUE);  copymat<-matrix(Out$Trans[i,], nrow=K, byrow=TRUE)
                  for(Row in 1:K) { for(Col in 1:K){ newmat[Row,Col]<-copymat[ permPARS[Row], permPARS[Col]]  }}
                  finQ[i,]<-as.vector(t(newmat))  
                    }   
            return(list("Means"=finMu, "Trans"=finQ,"q0"=finq0, "Var"=finSig, "States"=finS, "Y"=Out$Y))}

        res.fit<-TrimThin(Result,burn=50, newN=200)
        res.us<-unswStates(Out=res.fit, K=K)
      
        # reorder by final group size? or mean?
        #MedianZ<-factor(apply(res.us$States, 2, median), levels=1:K)
        MedianZ<-apply(res.us$q0, 2, mean)

        finalReord<-order(MedianZ, decreasing=TRUE)
        finalRpars<-c(1:K)[order(finalReord)]

        omeans<- res.us$Means[,finalReord]
        ovariance<-res.us$Var[,finalReord] ## NEW
        oq0<-res.us$q0[,finalReord]
        ostates<-res.us$States 
        for(i in 1:dim(ostates)[1]){           
           Snow2<-factor(ostates[i,], levels=1:K)
           levels(Snow2)<-finalRpars                  
           ostates[i,]<-as.numeric(as.character(Snow2))    }

          oTrans<-res.us$Trans
        for (i in 1:dim(res.us$Trans)[1]){ #UNSWITCH Q 
          newmat<-matrix(res.us$Trans[i,], nrow=K, byrow=TRUE);  copymat<-matrix(res.us$Trans[i,], nrow=K, byrow=TRUE)
          for(Row in 1:K) { for(Col in 1:K){ newmat[Row,Col]<-copymat[ finalRpars[Row], finalRpars[Col]]  }}
          oTrans[i,]<-as.vector(t(newmat)) }


          # rename true Zs too to compare  # reorder by final group means
        truSize<-table(YZ$States)/length(YZ$S)
        ReordY<-order(truSize, decreasing=TRUE); Reorpars<-c(1:K)[order(ReordY)]          
        yStates<-factor(YZ$Sta, levels=1:K);        levels(yStates)<-c(Reorpars, K-1,K)                  
        yStates<-as.numeric(as.character(yStates))     
       
       
       #compute means of q0, sd and ratio  
        q0Oversd<-apply(oq0,2,mean)/apply(oq0, 2, sd)
        q0Ratio<-sum( q0Oversd >1) # changing it to 1?
        
        # group size
        groupSize<-nk(ostates[,-1], K=K)/n
        trueGsize<-table(YZ$States)[]/n
        
        # predicted:
          ordX<-function (x) names(x[which.max(x)])
          predZ<-  unlist(lapply(apply(ostates, 2,table), ordX))
          RScore<-sum(yStates==predZ[-length(predZ)])/length(yStates) 
          YZpred<-data.frame("Y"= YZ$Obs, "Zpred"=as.numeric(predZ[-length(predZ)]),"trueZ"=yStates )  # is this right or shud be -1?
    
          Zhat<-predZ
          Mhat<-apply(omeans, 2, mean)
          Vhat<-apply(ovariance, 2, mean)
          Qhat<-apply(oTrans, 2, mean)
          q0hat<-apply(oq0, 2, mean)
          all<-list("Zhat"=Zhat, "Mhat"=Mhat,"q0hat"=q0hat, "Qhat"=Qhat , "Vhat"=Vhat)

        return(list(  "Zhat"=Zhat, "Mhat"=Mhat,"q0hat"=q0hat, "Qhat"=Qhat , "Vhat"=Vhat,"RandScore"=RScore, "YZpred"=YZpred$Zpred))
          }

#' PostProcessing function univ
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope



Zhmm_PostProc<-function( Grun, mydata, burn=1000, Thin=1, prep=1000, isSim=TRUE, simlabel="sim"){	
	#	maxZ<-function (x)  {as.numeric(names(which.max(table( x ))))}
	
	 Grun<-TrimThin(Grun, burn, Thin)		
	 #extract Y's from sim or data:
	 ifelse(isSim==TRUE, Y<-mydata$Observed,  Y<-mydata)
	 
	targetK0<-Grun$K0
	K0<-as.numeric(names(table(targetK0 )))
	n<-length(Y)  
	K<-dim(Grun$q0)[2]	
	p_vals<-data.frame("K0"=K0, "PropIters"=as.numeric(table(Grun$K0)/dim(Grun$q0)[1]), "RAND"=NA, "MAE"=NA, "MSE"=NA,"Pmin"=NA, "Pmax"=NA, "Concordance"=NA, "MAPE"=NA, "MSPE"=NA)
	## 1. split by K0
	K0estimates<-vector("list", length(K0))
	Zestimates<-vector("list", length(K0))
	USfullrun<-vector("list", length(K0))
	#for each K0:
	for ( .K0 in 1:length(K0)){
		grunK0<-Grun
		# split data by K0
		.iterK0<- c(na.omit(c(1:dim(Grun$q0) [1])[ targetK0  ==K0[.K0]]))
		if(length(.iterK0)>100){
		grunK0$Mu<-	Grun$Means[.iterK0,]
		grunK0$Q<-	Grun$Trans[.iterK0,]
		grunK0$q0<-	Grun$q0[.iterK0,]
		grunK0$MAP<-	Grun$MAP[.iterK0]
		grunK0$Z<-	Grun$States[.iterK0,]
		grunK0$K0<-	Grun$K0[.iterK0]

		## 2. unswitch
		grunK0us<-Zswitch_hmm(grunK0,0.05 )			
		Zetc<-Zagg_HMM(grunK0us, Y)
		USfullrun[[.K0]]<-grunK0us

		# PLOTS
		p1<-ggplot(data=grunK0us$Pars, aes(x=q0, fill=factor(k))) + geom_density( alpha=0.4)+ggtitle("Stationary distribution ")+ylab("")+xlab("")  +  theme(legend.position = "none")
		p2<-ggplot(data=grunK0us$Pars, aes(x=mu, fill=factor(k))) + geom_density( alpha=0.4)+ggtitle("Means")+ylab("")+xlab("") +  theme(legend.position = "none")
		p3<-ggplot(data=grunK0us$Pars, aes(x=X5, fill=factor(k))) +geom_density(alpha=0.4)+ggtitle("One column of Q")+ylab("")+xlab("") +  theme(legend.position = "none")
		grobframe <- arrangeGrob(p1, p2, p3, ncol=3, nrow=1,main = textGrob(paste(simlabel,": posterior parameter estimates for", K0[.K0]," groups"), gp = gpar(fontsize=8, fontface="bold.italic", fontsize=14)))
		ggsave(plot=grobframe, filename= paste("PosteriorParDensities_",simlabel,"_K0", K0[.K0],".pdf", sep="") , width=20, height=7, units='cm')

		## 3. RAND, MSE	
		if(isSim==TRUE){	Zhat<- Zetc$Zpred
					p_vals$RAND[.K0]<-(sum(mydata$States==Zhat)/n) 
		} else {			p_vals$RAND[.K0]<-'NA'}

		K0estimates[[.K0]]<-cbind(Zetc$theta, "K0"=K0[.K0])
		Zestimates[[.K0]]<-Zetc$Zpred
		## 4. Predict replicates
		postPredTests<-PostPredFunk( grunK0us,Zetc, Y, prep, simlabel)
		# store output in p_vasl
		p_vals$MAE[.K0]<- Zetc$MAE
		p_vals$MSE[.K0]<- Zetc$MSE
		p_vals$Pmin[.K0]<-postPredTests$MinP
		p_vals$Pmax[.K0]<-postPredTests$MaxP
		p_vals$MAPE[.K0]<-postPredTests$MAPE
		p_vals$MSPE[.K0]<-postPredTests$MSPE
		p_vals$Concordance[.K0]<-1-postPredTests$Concordance	}	}
	return(list(p_vals, K0estimates,Zestimates, USfullrun))		      }
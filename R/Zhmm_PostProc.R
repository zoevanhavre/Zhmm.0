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
		if(length(.iterK0)>200){
		grunK0$Mu<-	Grun$Means[.iterK0,]
		grunK0$Q<-	Grun$Trans[.iterK0,]
		grunK0$q0<-	Grun$q0[.iterK0,]
		grunK0$MAP<-Grun$MAP[.iterK0]
		grunK0$Z<-	Grun$States[.iterK0,]
		grunK0$K0<-	Grun$K0[.iterK0]

		## 2. unswitch
		grunK0us<-Zswitch_hmm(grunK0,0.05 )			
		Zetc<-Zagg_HMM(grunK0us, Y)
		USfullrun[[.K0]]<-grunK0us

		# PLOTS
		p1<-ggplot(data=grunK0us$Pars, aes(x=q0, fill=factor(k))) + geom_density( alpha=0.4)+ylab("")+xlab("")  +  theme(legend.position = "none")+
		ggtitle(bquote(atop(italic(paste( "p(K=", .(K0[.K0]), ")=", .(round(p_vals$PropIters[.K0],2)), sep="")), atop("Stationary distribution"))))

		p2<-ggplot(data=grunK0us$Pars, aes(x=mu, fill=factor(k))) + geom_density( alpha=0.4)+ylab("")+xlab("") +  theme(legend.position = "none")+
		ggtitle( bquote( atop(italic( paste(.(simlabel) , ", K=", .(K0[.K0]) ) ), atop("Means"))))
	
	#	p3<-ggplot(data=grunK0us$Pars, aes(x=X5, fill=factor(k))) +geom_density(alpha=0.4)+ggtitle("One column of Q")+ylab("")+xlab("") +  theme(legend.position = "none")
	#	grobframe <- arrangeGrob(p1, p2, p3, ncol=3, nrow=1,main = textGrob(paste(simlabel,": posterior parameter estimates for", K0[.K0]," groups"), gp = gpar(fontsize=8, fontface="bold.italic", fontsize=14)))
	#	ggsave(plot=grobframe, filename= paste("PosteriorParDensities_",simlabel,"_K0", K0[.K0],".pdf", sep="") , width=20, height=7, units='cm')


		## 3. RAND, MSE	
		if(isSim==TRUE){	Zhat<- Zetc$Zpred
					p_vals$RAND[.K0]<-(sum(mydata$States==Zhat)/n) 
		} else {			p_vals$RAND[.K0]<-'NA'}

		K0estimates[[.K0]]<-cbind(Zetc$theta, "K0"=K0[.K0])
		Zestimates[[.K0]]<-Zetc$Zpred
		
#NEW
		# clust probabilities:
p3<-HmmAllocationPlot(outZ=grunK0us$Z[,-(n+1)], myY=Y)
		## 4. Predict replicates
		postPredTests<-PostPredFunk( grunK0us,Zetc, Y, prep, simlabel)
		p4<-postPredTests$ggp
		# store output in p_vasl
		p_vals$MAE[.K0]<- Zetc$MAE
		p_vals$MSE[.K0]<- Zetc$MSE
		p_vals$Pmin[.K0]<-postPredTests$MinP
		p_vals$Pmax[.K0]<-postPredTests$MaxP
		p_vals$MAPE[.K0]<-postPredTests$MAPE
		p_vals$MSPE[.K0]<-postPredTests$MSPE
		p_vals$Concordance[.K0]<-1-postPredTests$Concordance		


# clusters:
plotyz<-data.frame("X"=1:n, "Y"=Y, "Post_Z"=Zetc$Zpred[-(n+1)])
p5<-ggplot(plotyz, aes(y=Y,x=X ))+geom_point(aes(colour=Post_Z),)+ geom_line(alpha = 1/4)+theme_bw()+  theme(legend.position = "none")+ggtitle("Posterior Allocations")+xlab("Time") 

		pdf( file= paste("pp_", simlabel, "K_ ",K0[.K0] , ".pdf",sep='') ,width=10, height=5)
	 		print( wq::layOut(	list(p1, 	1, 1:2),  
		        	list(p2, 	1, 3:4),   
		         	list(p3,	1,5:6),
		         	list(p5, 	2,1:3),  
		          	list(p4, 	2,4:6)))
		dev.off()
		} }


	return(list(p_vals, K0estimates,Zestimates, USfullrun))		      }
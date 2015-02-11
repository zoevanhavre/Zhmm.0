#' PostProcessing Predictive function
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope

PostPredFunk<-function(.GrunK0us, .Zetc, .Y, .prep , .simlabel){
			#Y<-.GrunK0us$Y
				n<-length(.Y)
				K<- max(.GrunK0us$Pars$k)
			   .GrunK0us$Pars$k<-factor(.GrunK0us$Pars$k, levels=c(1:max(.GrunK0us$Pars$k)))
				#swq0<- reshape(.GrunK0us$Pars, v.names="q0", idvar="Iteration", timevar="k", direction='wide')[,-1]
				#swMeans<- reshape(.GrunK0us$Pars, v.names="mu", idvar="Iteration", timevar="k", direction='wide')[,-1]
				

			#	swq0<-swq0[,(dim(swq0)[2]-K+1) :dim(swq0)[2]]
			#	swMeans<-swMeans[,(dim(swMeans)[2]-K+1) :dim(swMeans)[2]]

				#swQ<- reshape(.GrunK0us$Pars, v.names="Sig", idvar="Iteration", timevar="k", direction='wide', drop=c("Mu", "P"))[,-1]
			
				#DrawIters<-function(x) sample(c(1:max(.GrunK0us$Pars$Iteration)), size=x, replace = T, prob = NULL)
				.iters<-sample(c(1:max(.GrunK0us$Pars$Iteration)), size=.prep, replace = T, prob = NULL)
								
				# apply to .iters :   draw Z and do rnorm
				DrawRepY<-function(x){ 
 					drawPars<-subset(.GrunK0us$Pars, Iteration==x)
 					Q<-as.matrix(drawPars[, c((4+1):(K+4))])
 					rownames(Q)<-NULL
					SimHMM(Q, drawPars$mu, drawPars$q0,n)}
				
				.yrep<-matrix( nrow=.prep, ncol=n)
				.zrep<-matrix( nrow=.prep, ncol=n)
				for(i in 1:.prep){
					xy<-DrawRepY(.iters[i])
					.yrep[i,]<-xy[,2]
					.zrep[i,]<-xy[,1]
				}
				#.yzrep<-sapply(.iters, DrawRepY)
				#.yrep<-matrix(.yzrep[1,],nrow=.prep, byrow=T)
				#.zrep<-matrix(.yzrep[2,],nrow=.prep, byrow=T)

				## calculate various values
				# min/max
				MinP<-sum(apply(.yrep, 1, min) < min(.Y))/.prep
				MaxP<-sum(apply(.yrep, 1, max) > max(.Y))/.prep

				# Prediction Concordance 
				ComputePredConcordance<-function(x){sum( (x< quantile(.Y, .025)) | (x > quantile(.Y, 1-.025))  ) /n}
				.pc<-apply(.yrep, 1, ComputePredConcordance)
				#p_vals$Concordance[.K0]<-paste(mean(.pc), " (",quantile(.pc, .025), ",", quantile(.pc, 1-.025), ")", sep="")
				Concordance<-mean(.pc)
		 		
		 		# 4.2 MSPE
				# take Z matrix and replace with estimated mean
				Zemu<-.zrep
		     	.PosteriorMeans<-.Zetc$theta$value[.Zetc$theta$variable=="mu"]
			    Zemu<-apply( Zemu, c(1,2), function(x) {return(.PosteriorMeans[x])} )
			        
					MSPE_dist<-apply((.yrep-Zemu)^2, 1, sum)
					MAPE_dist<-apply(abs(.yrep-Zemu), 1, sum)

				MSPE<-mean(MSPE_dist)
				MAPE<-mean(MAPE_dist)

				### 4.3 Plot data VS replicates	
				predplot<-ggplot(data.frame("Y"=.Y, "n"=c(1:n)), aes(x=Y))  + 
				#geom_histogram(aes(y=..density..),  colour="red", fill="white")+
				geom_line(data=melt(.yrep),stat="density", aes(x=value,group=Var1), size=0.5, color="blue", alpha=0.1)+
				geom_density(color="green", size=1, linetype="dashed")+ geom_hline(yintercept=0, colour="white", size=1)+
				theme_bw()+ggtitle("Predicted densities")
				#ggsave(plot=predplot, filename= paste("PredictiveDensities_",.simlabel,"_K0",K,".pdf", sep="") ,width=10, height=10, units='cm' )

				#ggsave(plot=predplot, filename= paste("PredictiveDensities_",.simlabel,"_K0",K,".bmp", sep="") )
				
				return(list( "MinP"=MinP, "MaxP"=MaxP, "MAPE"=MAPE,  "MSPE"=MSPE, "Concordance"=Concordance, "ggp"=predplot))
			}

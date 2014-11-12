 #' Function to plot the number of non-empty groups
#'
#' ...
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))


plotK0_HMM<-function(runResult,burn=1000 , Ptitle=" ", PlotName="nonemptyK0.pdf"){

          alphaMAX<-(K-1)*(1+K-2+alphaMin)*(1+1/( (1/2) - alphaMin*(K-1))) -(K-1)*alphaMin+0.1
          AllAlphas<-matrix(nrow=J, ncol=K)
          AllAlphas[,1]<-alphaMAX
         # AllAlphas[,2:K]<-seq(alphaMAX, alphaMin, length=J)
          AllAlphas[,2:K]<-      c( alphaMin+(alphaMAX-alphaMin)/c(1:(J-1))^3, alphaMin)



	runK0<-sapply(  runResult[["Zs"]],  function(x)  { apply( x, 2,function(x) sum(table(x )>0) )}   ) 
	.rk0<-melt(runK0[-1:-burn,])
	names(.rk0)<-c("Iteration", "Chain", "K0"); 
	.rk0$Chain<-as.factor(.rk0$Chain)

	alphas<-as.character(signif( c(60, 50, 40,30, 20, 10, 5, 3, 2, 1, 0.5, 1/2^(c(2:10, 12, 14, 18, 20, 25, 30, 35, 40, 45, 50 ))), 2))
	alphas<-data.frame( "Alpha"=alphas, "Chain"=c(1:length(alphas)))

	gp<-ggplot(.rk0, aes(y=K0,x=Chain ))+geom_boxplot()+ ggtitle(Ptitle)+xlab("Prior Tempering Chain")+ylab("# Non-Empty Groups")+theme_bw()+scale_y_continuous(breaks=c(1:10))+  theme(panel.grid.minor=element_blank())+geom_text(data = alphas, aes(x = Chain, y = 0, label = Alpha),  size = 3, angle = 90, colour='red') + annotate("text", x = 0.4, y = 0,  size = 3,angle=90, label = "Alpha", colour="red")
	ggsave(gp, file=PlotName)
	}

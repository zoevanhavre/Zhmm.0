 #' MuXqoPlot
#'
#' Plot means by stationary distribution
#' @param ...
#' @keywords ...
#' @export
#' @examples ...

MuXqoPlot<-function(gibbsout, burn, trueValues="NA", minq=0.01, plotlab=""){
		both<-cbind(melt(gibbsout$Means[-c(1:burn),])[,3],
		   			melt(gibbsout$q0[-c(1:burn),])[,3])
		
		color <- as.factor(melt(gibbsout$Means[-c(1:burn),])[,2])
		raincol<-rainbow(length(levels(color)))
		levels(color)<-raincol
		trancol<-sapply(c(1:length(color)), function(x) adjustcolor(color[x], alpha.f=both[x,2]))
		minmaxMEANS<-c(min(both[both[,2]>minq ,1]), max(both[both[,2]>minq ,1]))

		plot(  melt(gibbsout$Means[-c(1:burn),])[,3],
		 	   melt(gibbsout$q0[-c(1:burn),])[,3],
		           # col=rgb(0,0,0,alpha=melt(gibbsout$q0[-c(1:burn),])[,3]),
		           col=trancol, xlim=minmaxMEANS,
		            xlab="Mean", ylab="Stationary Distribution", bg='grey',
		            main=paste(plotlab, ",zoomed, with transparency ~ q0 ", sep=" "))

		if(trueValues!="NA"){

			points(trueValues, pch=7, cex=2)
		}}


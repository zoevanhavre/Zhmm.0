 #' MuXqoPlot
#'
#' Plot means by stationary distribution
#' @param ...
#' @keywords ...
#' @export
#' @examples ...

MuXqoPlotALL<-function(gibbsout, burn, trueValues=NA, plotlab=""){
		plot(  melt(gibbsout$Means[-c(1:burn),])[,3],
		 	   melt(gibbsout$q0[-c(1:burn),])[,3],
		            #col=rgb(0,0,0,alpha=melt(gibbsout$q0[-c(1:burn),])[,3]),
		            xlab="Mean", ylab="Stationary Distribution", bg='grey',
		             main=paste(plotlab, ", Means X q0 (All)", sep=" "))

		if(trueValues!="NA"){
			points(trueValues,col='green', pch=7, cex=1)
		}}

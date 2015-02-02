#' Trim n Thin
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope

 TrimThin<-function(out, burn){
	                ids<- burn:length(out$MAP)  # TRIM
	              if(Thin>1){  ids<-seq(from=min(ids), to=max(ids), by=Thin)   } #THIN
	                #save selected iterations
	                MU<-out$Means[ids,]
	                Q<-out$Trans[ids,]
	                q0<-out$q0[ids,]
	                S<-out$States[ids,]
	                MAP<-out$MAP[ids]
	                K0<-out$K0[ids]
	                return(list("Means"=MU, "Trans"=Q,"K0"=K0,"q0"=q0, "States"=S, "Y"=out$Y, "MAP"=MAP))
	                }

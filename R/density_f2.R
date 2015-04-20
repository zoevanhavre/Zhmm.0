 #' Function to compute, at 1 iteration, 
#'
#' density of dirichlet
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))



density_f2<-function( y,  .mu ,.Q , .q0=NA){

y1<-y[1]
y2<-y[2]
k<-length(.mu)
.Q<-matrix(.Q,  ncol=k,  byrow=TRUE)

if(is.na(.q0)){ .q0<-ALTERNATEq0(.Q)}  
gridpoints<-expand.grid(c(1:k), c(1:k))

return(sum(apply(gridpoints, 1, function(x)   density_cell_f2(y1,  y2, q01=  .q0[x[1]],  .Q[x[1], x[2]],  .mu[x[1]], .mu[x[2]] )  ) ))
}


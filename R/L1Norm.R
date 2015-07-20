 #' gibbsHMM_PT
#'
#' parallel tempering with a column prior - option to mix over column or stick to j=1
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))


L1Norm<-function(y,TV, NOW, gridN=10, Plots=FALSE ){
	#TV list $z $q0 $Q $mu
	#NOW list  $z_now $q0_now $Q_now $mu_now
	y1<-y[1]
	y2<-y[2]
	kTRUE<-length(TV$mu)
	k<-length(NOW$mu)

	# # TRUE values
	z1<-TV$z[1]
	z2<-TV$z[2]
	Q<-matrix(TV$Q,  ncol=kTRUE,  byrow=TRUE)
	Q12<-Q[z1, z2]
	q01<-TV$q0[z1]
	mu1<-TV$mu[z1]
	mu2<-TV$mu[z2]

	# # NOW values (T = iteration t)
	z1_t<-NOW$z[1]
	z2_t<-NOW$z[2]
	Q_t<-matrix(NOW$Q,  ncol=k,  byrow=TRUE)
	Q12_t<-Q_t[z1_t, z2_t]
	q01_t<-NOW$q0[z1_t]
	mu1_t<-NOW$mu[z1_t]
	mu2_t<-NOW$mu[z2_t]

	# grid HERE
	griddy<-expand.grid(seq(min(y), max(y), length =gridN), seq(min(y), max(y), length =gridN))

	# Values given truth
	.truth<-apply(griddy, 1, function(x)   density_cell_f2(x[1],  x[2], q01,  Q12, mu1, mu2)  )

	# Values given Iteration
	.post<-apply(griddy, 1, function(x)   density_cell_f2(x[1],  x[2], q01_t,  Q12_t, mu1_t, mu2_t)  )
	if(Plots==FALSE){
	return(sum(abs(.post-.truth)))
} else{
try<-cbind(griddy, "Posterior"=.post, "Truth"=.truth, "Diff"=abs(.post-.truth))
t1<-ggplot(try, aes(x=Var1, y=Var2, fill=Truth))+geom_tile()+ggtitle("Truth")
t2<-ggplot(try, aes(x=Var1, y=Var2, fill=Posterior))+geom_tile()+ggtitle("Post.")
t3<-ggplot(try, aes(x=Var1, y=Var2, fill=Diff))+geom_tile()+ggtitle("| Post.-Truth |")

	return(list(sum(abs(.post-.truth)), t1, t2, t3))
}
	}

#' Processing output of label switching function
#'
#' explain here
#' @param output of HMM label switching function, computes parameter estimates and some error values. 
#' @keywords postprocessing
#' @export
#' @examples
#' #nope

Zagg_HMM<-function(USout, .Y = Y) {
.par <- melt(USout$Pars, id.vars = c("Iteration", "k"))
theta <- aggregate(value ~ variable + factor(k), mean, data = .par)
K <- max(.par$k)
Zhat<- factor( apply(USout$Z, 2,maxZ))[-((length(.Y))+1)]
Zemu <- as.numeric(Zhat)
.Mus <- theta$value[theta$variable == "mu"]
for (i in 1:length(Zemu)) { Zemu[i] <- .Mus[as.numeric(Zhat[i])]  }
MSE <- sum((.Y - Zemu^2) )
MAE <- sum(abs(.Y - Zemu))
list(theta = theta, Zpred = Zhat, MSE = MSE, MAE = MAE)
				}
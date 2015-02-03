 #' Compute median state from unswitched MCMC samples from Zmix
#'
#' Compute median state
#' @param x
#' @keywords states
#' @export
#' @examples #something

maxZ<-function (x)  as.numeric(names(which.max(table( x ))))

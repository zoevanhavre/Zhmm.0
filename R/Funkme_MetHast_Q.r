 #' Funkme_MetHast_Q 
#'
#' text...
#' @param x, alpha, log=False
#' @keywords dirichlet
#' @export
#' @examples dDirichlet(c(.1, .9), c(0.1,0.1))

# Funkme_MetHast_Q<-function(.qnew, .qold,  .z1){          # NOT DONE! - pick right Qold and Z
#                     K<-dim(.qnew)[1]
#                     q0Previous<-ALTERNATEq0(matrix(.qold, nrow=K, byrow=TRUE))              
#                     .q0new<-ALTERNATEq0(.qnew)           
#                     A<-.q0new[.z1]/q0Previous[.z1]  
#                     if(A=='NaN'){A<-0}  #0 or 1?
#                     if(A>runif(1,c(0,0.99))){
#                     return(.qnew)
#                     } else {
#                   return(matrix(.qold, nrow=K, byrow=TRUE))
#                 }} 
Funkme_MetHast_Q<-function(.q0new, .q0old,  .z1){          # NOT DONE! - pick right Qold and Z        
                    A<-.q0new[.z1]/.q0old[.z1]  
                    if(A=='NaN'){A<-0}  #0 or 1?
                    if(A>runif(1,c(0,0.99))){
                    return("TRUE")
                    } else {
                  return("FALSE")
                }} 
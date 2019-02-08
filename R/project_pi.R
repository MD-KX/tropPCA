#' Porjected points
#'
#' @param D_s The distance matrix used to build tropical space
#' @param D Distance vector (matrix)
#'
#' @return
#' @export
#'
#' @examples
project_pi<-function(D_s,D){
  if(is.null(dim(D_s))){
    lambda <- min(D - D_s)
    pi_D <- c(t(lambda + t(D_s)))
  }else{
    lambda <- apply(D - D_s, 2, min)#D_s by row
    pi_D <- apply(t(lambda + t(D_s)),1,max)
  }

  # pi_D <- ifelse(rep(is.null(dim(D_s)),28), c(t(min(D - D_s)  + t(D_s))), apply(t(apply(D - D_s, 2, min) + t(D_s)),1,max))

  # pi_D <- ifelse(rep(is.null(dim(D_s)),28), c(t(min(D[[1]] - D_s)  + t(D_s))), apply(t(apply(D[[1]] - D_s, 2, min) + t(D_s)),1,max))

  return(pi_D)
}

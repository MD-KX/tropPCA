#' Tropical distance between two trees.
#'
#'
#' @description
#' function to get tropical distance of two points
#' @export
#'
#' @examples

tropical_dist<-function(D_1,D_2){
  e <- length(D_1)
  t_dist <- 0
  for(i in 1:(e-1)){
    for(j in (i+1):e){
      if(abs(D_1[i]-D_2[i]-D_1[j]+D_2[j])>t_dist){
        t_dist<-abs(D_1[i]-D_2[i]-D_1[j]+D_2[j])
      }
    }
  }
  t_dist
}

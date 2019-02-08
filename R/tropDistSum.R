#' The sum of tropical distance of multiple trees
#'
#' @param pc_base
#' Three disntance vetors used to build tropical space
#' @param distVec_all
#' All of the distance vectors.
#' @return
#' @export
#'
#' @examples
tropDistSum <- function(pc_base, distVec_all){

  proj_points <- parLapply(cl, distVec_all, project_pi , D_s = pc_base)
  tropical_dist_vec <- mapply(tropical_dist, distVec_all, proj_points)
  sum_dist <- sum(tropical_dist_vec)


  return(sum_dist)
}

#' Tropical MCMC
#'
#' @param distVect_all
#' All of the distance vectors.
#' @param N
#'  Number of points in tropical space.
#' @param pcs
#' Number of pricinpal components
#' @param nr
#' Number of repeat times.
#' @param env
#' Parameter for changing environment.
#' @return
#' @export
#'
#' @examples
tropMCMC <- function(distVect_all, N, pcs, nr = 2, env = .GlobalEnv){
  env$sumsValue <- rep(NA, N)
  env$comb_list <- list()
  env$points_list <- list()
  env$trop_list <- list()
  D_all <- matrix(unlist(distVec_all), ncol=N)
  for(j in 1:nr){
    sample_init <- sample(N,pcs)
    best <- 100000
    out <- c(1:N)[-sample_init]
    pc_base_init <- D_all[,sample_init]
    init_value <- tropDistSum(pc_base_init, distVec_all)
    for(i in 1:N){
      change_ind <- sample(pcs,1)
      out_change <- sample(out, 1)
      comb_set <- c(sample_init[-change_ind], out_change)
      # Here, if we have some other constrains, then just add them
      u <- runif(1)

      new_base <- D_all[,comb_set]

      update_value <- tropDistSum(new_base, distVec_all)
      r <- init_value/update_value

      if(u < min(r, 1)){
        sample_init <- comb_set
        if(update_value < best){
          proj_points<-lapply(distVec_all,project_pi,D_s=new_base)
          tropicalDist_vec<-mapply(tropical_dist,distVec_all,proj_points)
          label_output<-list(sample_init,proj_points,tropicalDist_vec,update_value)
          best <- update_value
        }
      }
      out <- out[-which(out==out_change)]
      init_value <- update_value
      if(length(out)==0) break

      print(c(i,u < min(r, 1), best ,length(out)))
    }
    env$comb_list[[j]] <- label_output[[1]]
    env$points_list[[j]] <- label_output[[2]]
    env$trop_list[[j]] <- label_output[[3]]
    env$sumsValue[j] <- label_output[[4]]
  }
}

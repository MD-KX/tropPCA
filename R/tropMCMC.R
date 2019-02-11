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
  env$sumsValue <- rep(NA, nr)
  env$comb_list <- list()

  D_all <- matrix(unlist(distVec_all), ncol=N)
  for(j in 1:nr){
    sample_init <- sample(N,pcs)
    best <- 100000
    out <- c(1:N)[-sample_init]
    pc_base_init <- D_all[,sample_init]
    init_value <- tropDistSum(pc_base_init, distVec_all)
    while(length(out)!=0){
      change_ind <- sample(pcs,1)
      out_change <- sample(out, 1)
      comb_set <- c(sample_init[-change_ind], out_change)

      new_base <- D_all[,comb_set]

      update_value <- tropDistSum(new_base, distVec_all)
      r <- init_value/update_value

      if(runif(1) < min(r, 1)){
        sample_init <- comb_set

        best <- ifelse(update_value < best, update_value, best)
      }
      out <- out[-which(out==out_change)]
      init_value <- update_value

    }
    env$comb_list[[j]] <- sample_init
    env$sumsValue[j] <- best
  }
  return(list(comb_list, sumsValue))
}

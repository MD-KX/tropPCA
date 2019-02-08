#' Distance matrix
#'
#' @param trees
#' Phylogeny trees
#' @param tipOrder
#' The names of tip label. You cold give an order to them.
#' @description
#' This function is to get the distance vector of a phylogeny tree in tropical space
#' @return
#' @export
#'
#' @examples
distMat <- function(trees, tipOrder){ # Here trees should be a list
  if(class(trees)=="multiPhylo"){
    trees_root <- root(trees, outgroup = tipOrder[1],resolve.root=TRUE)

    chronotrees <- parLapply(cl, trees_root, chronos)
    dist_chrono <- parLapply(cl, chronotrees,cophenetic)

    dist_ordered <- parLapply(cl, dist_chrono, function(x) x[tipOrder, tipOrder])
    distVec_all <- parLapply(cl, dist_ordered,vec_fun)

    # chronotrees <- lapply( trees_root, chronos)
    # dist_chrono <- lapply(chronotrees,cophenetic)
    #
    # dist_ordered <- lapply( dist_chrono, function(x) x[tipOrder, tipOrder])
    # distVec_all <- lapply( dist_ordered,vec_fun)

  }else {
    treeOne <- root(trees, outgroup = tipOrder[1],resolve.root=TRUE)
    chronoTree <- chronos(treeOne)
    dist_chrono_one <- cophenetic(chronoTree)

    dist_ordered_one <- dist_chrono_one[tipOrder, tipOrder]
    distVec_all <- vec_fun(dist_ordered_one)
  }

  return(distVec_all)
}

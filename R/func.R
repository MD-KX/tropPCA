#function to generate random trees
dist_vector<-function(k,m){
  vecTreesArray<-array(NA,dim = c(k,choose(m,2))) 
  tree_store<-list()
  for(p in 1:k){
    t <- rcoal(m) ### generate an equidistant tree 
    tree_store[[p]]<-t
    t$tip.label<-seq(1,m,by=1)  ###change the labels
    
    D1 <- cophenetic(t) ### Compute the distant matrix
    
    ###vectorize the distance matrix
    for(row.num in 1:(m-1)){
      for(col.num in (row.num+1):m){
        vecTreesArray[p,col.num-row.num+(m-1+(m-1-row.num+2))*(row.num-1)/2]<-D1[row.num,col.num]  
      }
    }   
  }
  return(list(vecTreesArray,tree_store))
}

#function to get tropical distance of two points
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

#function to project tree D onto tconv(D_s) with s PC's and m leaves
# project_pi<-function(s,m,D_s,D){
#   lambda<-rep(NA,s)
#   for(i in 1:s){
#     lambda[i]<-min(D-D_s[i,])#D_s by row
#   }
#   
#   pi_D<-rep(NA,choose(m,2))
#   for(i in 1:choose(m,2)){
#     pi_D[i]<-max(lambda+D_s[,i])
#   }
#   return(pi_D)
# }

project_pi<-function(D_s,D){
  if(is.null(dim(D_s))){
    lambda <- min(D - D_s) 
    pi_D <- c(t(lambda + t(D_s))) 
  }else{
    lambda <- apply(D - D_s, 2, min)#D_s by row
    pi_D <- apply(t(lambda + t(D_s)),1,max)
  }
  return(pi_D)
}

vec_fun<-function(x){
  m<-dim(x)[1]
  vecTreesVec<-rep(NA,choose(m,2)) 
  for(row.num in 1:(m-1)){
    for(col.num in (row.num+1):m){
      vecTreesVec[col.num-row.num+(m-1+(m-1-row.num+2))*(row.num-1)/2]<-x[row.num,col.num]  
    }
  }
  vecTreesVec
}
# Function ----------------------------------------------------------------

tropDist <- function(pc_base, distVec_all){

  proj_points <- lapply(distVec_all, project_pi , D_s = pc_base)
  tropical_dist_vec <- mapply(tropical_dist, distVec_all, proj_points)
  sum_dist <- sum(tropical_dist_vec)


  return(sum_dist)
}

# tropDist <- function(pc_base, distVec_all){
#   
#   proj_points <- lapply(distVec_all, project_pi , D_s = pc_base)
#   tropical_dist_vec <- mapply(tropical_dist, distVec_all, proj_points)
#   sum_dist <- sum(tropical_dist_vec)
#   
#   label_output<-list(pc_base,proj_points,tropical.dist.vec,sum_dist)
# 
#   return(label_output)
# }  

# Build a new function
distMat <- function(phyTrees, tipOrder, outgp){ # Here trees should be a list
  if(class(phyTrees)=="multiPhylo"){
    # for(i in 1:length(phyTrees)){
    #   pickUp <- try(trees_root_ori[[i]] <- root(phyTrees[[i]], outgroup = outgp,resolve.root=TRUE), TRUE)
    #   if(class(pickUp)=="try-error") { next } else { trees_root_ori[[i]] <- pickUp }
    # }
    # null_ind <- sapply(trees_root_ori, is.null)
    # trees_root <- trees_root_ori[!null_ind]
    trees_root_ori <- lapply(phyTrees, function(x) tryCatch(root(x, outgroup = outgp, resolve.root=TRUE), error = function(x) NULL))
    null_ind <- sapply(trees_root_ori, is.null)
    trees_root <- trees_root_ori[!null_ind]
    
    chronotrees <- lapply(trees_root, chronos)
    dist_chrono <- lapply(chronotrees,cophenetic)
    
    dist_ordered <- lapply(dist_chrono, function(x) x[tipOrder, tipOrder])
    distVec_all <- lapply(dist_ordered,vec_fun)
    
  }else {
    treeOne <- root(phyTrees, outgroup = outgp,resolve.root=TRUE)
    chronoTree <- chronos(treeOne)
    dist_chrono_one <- cophenetic(chronoTree)

    dist_ordered_one <- dist_chrono_one[tipOrder, tipOrder]
    distVec_all <- vec_fun(dist_ordered_one)
  }

  return(distVec_all)
}


tropMCMC <- function(distVect_all, N, pcs, nr, env = .GlobalEnv){
  env$sumsValue <- rep(NA, nr)
  env$comb_list <- list()
  env$points_list <- list()
  env$trop_list <- list()
  D_all <- matrix(unlist(distVec_all), ncol=N) 
  for(j in 1:nr){
    sample_init <- sample(N,pcs)
    best <- 100000
    out <- c(1:N)[-sample_init]
    pc_base_init <- D_all[,sample_init]
    init_value <- tropDist(pc_base_init, distVec_all)
    for(i in 1:N){
      change_ind <- sample(pcs,1)
      out_change <- sample(out, 1)
      comb_set <- c(sample_init[-change_ind], out_change)
      # Here, if we have some other constrains, then just add them
      u <- runif(1)
      
      new_base <- D_all[,comb_set]
      
      update_value <- tropDist(new_base, distVec_all)
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

# Tropical tree triangle
plot.trop.triangle <- function(D){
  k <- ncol(D)
  pdf("Tropical_Triangle.pdf")
  plot(D[1,],D[2,])
  for(i in 1:(k - 1)){
    for(j in (i + 1):k){
      tseg1 <- tropical.geodesic.dim.2(D[,i],D[,j])
      tseg2 <- tropical.geodesic.dim.2(D[,i],D[,j],flag=1)
      if(tseg1[[2]] < tseg2[[2]]) tseg <- tseg1
      else tseg <- tseg2
      segments(tseg[[1]][1,1],tseg[[1]][2,1],tseg[[1]][1,2],tseg[[1]][2,2],col= 'black')
      segments(tseg[[1]][1,2],tseg[[1]][2,2],tseg[[1]][1,3],tseg[[1]][2,3],col= 'black')
    }
  }
  points(x=proj_2D_plot_m[,2],y=proj_2D_plot_m[,3],pch=16,cex=0.75,col=freq)
  dev.off()
}
normalize.ultrametrices <- function(D){
  k <- ncol(D)
  new.D <- matrix(rep(0, 2*k), nrow=2, ncol=k)
  for(i in 2:3)
    new.D[i-1, ] <- D[i, ] - D[1, ]
  return(new.D)
}

# input: matrix D of size (s x e) whose rows are vertices of tropical polytope, point P in rowspan(D)
#output: point Q corresponding to P in column span of D
polytope_iso<-function(D, P){
  e = length(P)
  s = dim(D)[[1]]
  Q = mat.or.vec(1, s)
  for (i in seq(s)){
    maxvalue = D[i,1] - P[[1]]
    for (j in seq(e)){
      maxvalue = max(maxvalue, D[i,j] - P[[j]])
    }
    Q[[i]]=maxvalue
  }
  return(Q)
}
normalize.proj<-function(D){
  r<-length(D)
  D.new<-rep(NA,r)
  for(i in 1:r){
    D.new[i]<-D[i] - D[1]
  }
  return(D.new)
}
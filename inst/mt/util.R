#' @param n0 nombre de jeux de données non-informatif
#' @param n1 nombre de jeux de données avec information
#' @param n_grp  nombre de groupe
#' @param grp_size  taille des groupes
#' @param separability coefficient de séparabilité des groupes
scenario1 <- function(n0, n1, n_grp, grp_size, separability){

  get_one_data <- function(i) {
    noise <- rnorm(grp_size * n_grp, sd = 1)
    mu <- rep(i * (1:n_grp), each = grp_size)
    scale(mu + noise, TRUE, FALSE)
  }
  
  grp_levels <- rep(c(0, separability), c(n0, n1))
  
  lapply(grp_levels, get_one_data)
}

#' @param n0 nombre de jeux de données non-informatif
#' @param n1 nombre de jeux de données avec information
#' @param n_grp  nombre de groupe
#' @param n_grp_per_data  nombre de groupes propres à chaque dataset
#' @param grp_size  taille des groupes
#' @param separability coefficient de séparabilité des groupes
scenario2 <- function(n0, n1, n_grp, n_grp_per_data, grp_size, separability){

  get_one_data <- function(separability_) {
    
    noise <- rnorm(grp_size * n_grp, sd = 1)
  
    ## draw n.dif group per feature
    grp_means <- sample(1:n_grp, n_grp_per_data)
    
    ## how many times per feature (at least one)
    grp_means_sizes <- as.vector(rmultinom(1, size = n_grp - n_grp_per_data, prob = rep(1/n_grp_per_data, n_grp_per_data)) + 1)

    ## group means in all data sets    
    means <- sample(rep(grp_means, grp_means_sizes))
    
    mu <- rep(separability_ * means, each = grp_size)
    scale(mu + noise, TRUE, FALSE)
  }

  grp_levels <- rep(c(0, separability), c(n0, n1))
  
  lapply(grp_levels, get_one_data)
}

directClustering <- function(dataSets) {
  hclust(dist(do.call("cbind", dataSets), method = "euclidean"), method = "ward.D2")
}

averagedClustering <- function(dataSets) {
  AD <- Reduce("+", lapply(dataSets, dist, method = "euclidean")) / length(dataSets)
  hclust(AD, method = "ward.D2")
}

mergeTreesWard <- function(dataSets) {
  hc_list <- lapply(dataSets, FUN = function(x) {
    univarclust::ward_1d(x)
  })
  mergeTrees::mergeTrees(hc_list)
} 

oneDimClustering <- function(dataSet) {
  hclust(dist(dataSet, method = "euclidean"), method = "ward.D2")
}

## run approaches
oneRun <- function(sim_label, dataSets, reference, score = NID, k = 10){
  
  ## Merge Tree with Ward 1D
  time_MTW <- system.time(MTW_res <- mergeTreesWard(dataSets))[3]
    
  ## One dimension
  time_OD <- system.time(OD_res <- oneDimClustering(dataSets[[length(dataSets)]]))[3]
  
  # Direct Clustering
  time_DC <- system.time(DC_res <- directClustering(dataSets))[3]

  # Average Distance
  time_AD <- system.time(AD_res <- averagedClustering(dataSets))[3]

  nb_ind <- length(dataSets[[1]])

  ## spectral versions with random SVD
  time_rSVD <- system.time(
    {
      rSVD <- rsvd(do.call("cbind", dataSets), k = k) 
      dataSets_rspectral <- as.list(as.data.frame(rSVD$u %*% diag(rSVD$d)))  
    }
  )[3]


  # Spectral Merge trees
  time_rSMTW <- system.time(rSMTW_res <- mergeTreesWard(dataSets_rspectral))[3] + time_rSVD
  
  NID_1    <- apply(cutree(OD_res  , 1:nb_ind), 2, score, c2 = reference)
  NID_MTW  <- apply(cutree(MTW_res  , 1:nb_ind), 2, score, c2 = reference)
  NID_DC   <- apply(cutree(DC_res  , 1:nb_ind), 2, score, c2 = reference)
  NID_AD   <- apply(cutree(AD_res  , 1:nb_ind), 2, score, c2 = reference)
  NID_rSMTW <- apply(cutree(rSMTW_res, 1:nb_ind), 2, score, c2 = reference)
  
  do.call(rbind, 
          list(OD  = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="One"  , NID=NID_1    , time = time_OD),
               MTW  = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="MTW"   , NID=NID_MTW   , time = time_MTW),
               DC  = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="DC"   , NID=NID_DC   , time = time_DC),
               AD  = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="AD"   , NID=NID_AD   , time = time_AD),
               SMTW = data.frame(nb_grp=1:nb_ind, Sim=sim_label, method="SMTW"  , NID=NID_rSMTW  , time = time_rSMTW)
          ))

}

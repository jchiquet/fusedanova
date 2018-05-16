#####################################################-
# Transformation en matrice merge type hclust$merge #
#####################################################-
CreationMatriceMerge <- function(aggreg_mat){

  n <- nrow(aggreg_mat)
  aggreg_mat$down[aggreg_mat$down == 0] <- NA
  aggreg_mat <- aggreg_mat[order(aggreg_mat$actif, decreasing = TRUE), ]

  merge_mat = matrix(NA, ncol = 2, nrow = n - 1)
  iter <- 1
  while (iter < (n-1)) {

    which_group <- which(aggreg_mat$actif[1] == aggreg_mat$element) # label of current guy
    which_merge <- which(aggreg_mat$up[1]    == aggreg_mat$actif)   # index of guy with which it is merging
    
    merge_element <- aggreg_mat[which_merge, ]
    merge_mat[iter, ] <- c(aggreg_mat$element[1], merge_element$element)
    
    aggreg_mat$down[which_group] <- aggreg_mat$actif[which_merge]
    aggreg_mat[which_merge, ] <- c(iter + n, merge_element[c("actif", "up", "down")])
    aggreg_mat <- aggreg_mat[-1, ]

    iter <- iter + 1
  }

  # fusion des deux dernieres colonnes (il est cense n'en rester que deux)
  merge_mat[iter,] <- aggreg_mat$element[1:2]
  
  # codage singleton/fusion pour hclust (-/+)
  merge_mat[merge_mat <= n] <- -merge_mat[merge_mat <= n] 
  merge_mat[merge_mat > n]  <-  merge_mat[merge_mat > n] - n
  
  return(merge_mat)
}

#########################################
# Recuperation de l'ordre des individus #
#########################################

OrdreIndividus <- function(aggreg_mat){

  aggreg_mat$down[aggreg_mat$down == 0] <- NA
  aggreg_mat <- aggreg_mat[order(aggreg_mat$actif, decreasing = TRUE), ]
  
  group_parents <- split(aggreg_mat$actif, aggreg_mat$up)
  parents <- sort(unique(aggreg_mat$up))
  
  # Ordonner les groupes parents et leurs enfants
  first_order = c(0,group_parents[[1]])# la base, engendre tous les groupes/individus
  # deux fois 0, normal. c'est voulu.
  
  for (i in 2:length(group_parents)) {

    children <- group_parents[[i]]
    parent <- parents[i]
    
    indice_replace <- which(first_order == parent)
    
    if (length(indice_replace) > 0){
      first_order <- c(first_order[1:(indice_replace)], children,
                       first_order[(indice_replace + 1):length(first_order)]) 
    } else {
      first_order = c(first_order, parent, children)
    }
  }
  
  Order <- unique(first_order[!is.na(first_order)])
  Order <- rev(aggreg_mat$element)[Order + 1]
  Order
}

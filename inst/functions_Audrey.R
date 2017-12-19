library("fusedanova")

##################################################-
# Recuperation des regles a partir d'un objet FA #
##################################################-

# Parametres : 
# - n : nombre d'individus
# - FA :objet de type Fused-ANOVA

RecuperationRegles <- function(n, FA){
  iupCol <- 5 ; idownCol <- 4 ; splopCol <- 3
  lambdaCol <- 2 ; betaCol <- 1 ; lengthCol <- 6
  
  FATable <- FA@result[[1]]$table
  FATable <- as.matrix(FATable, ncol = 5)
  Length <- FATable[,iupCol] - FATable[,idownCol] + 1
  FATable <- cbind(FATable, Length)
  FATable <- FATable[order(FATable[,lambdaCol],FATable[,lengthCol], 
                           decreasing = FALSE),]
  FATable <- FATable[FATable[,lengthCol]!=1,]
  MatriceRegles <- matrix(0,nrow(FATable),6)
  MatriceVoisinDroite <- matrix(1:n,n,2,byrow = FALSE)
  for(i in 1:nrow(FATable)){
    x <- FATable[i,]
    idown <- x[idownCol]
    iup <- x[iupCol]
    Split <- MatriceVoisinDroite[idown,2]
    MatriceVoisinDroite[idown,2] <- iup
    MatriceRegles[i, ] <- c(1,x[lambdaCol],idown,Split, Split+1, iup)
  }
  
  MatriceRegles <- MatriceRegles[order(MatriceRegles[,2],FATable[,lengthCol], decreasing = TRUE),]
  return(MatriceRegles)
  
}

########################################################################-
# Recuperation du plus petit paquet entre les deux groupes d'une regle #
########################################################################-

# Parametre :
# - Vecteur : une ligne d'une matrice de regles

PlusPetitPaquet <- function(Vecteur){
  dimCol <- 1 ; lambdaCol <- 2 ; MinusCol <- 3
  SplitCol <- 4 ; Splitp1Col <- 5 ; MaxCol <- 6
  if((Vecteur[SplitCol]-Vecteur[MinusCol]+1)<(Vecteur[MaxCol]-Vecteur[Splitp1Col])){
    return(Vecteur[MinusCol]:Vecteur[SplitCol])
  }
  else{
    return(Vecteur[Splitp1Col]:Vecteur[MaxCol])
  }
  
}

#################################################################-
# Fonction permettant de savoir si une regle declenche un split #
#################################################################-

# Parametres :
# - RegleCourante : la regle que l'on considere
# - NumGroupeCourant : le numero du groupe cree s'il y en a un
# - NumeroGroupes : une matrice 2xn contenant pour chaque element
#                   son numero de groupe
# - Order : une liste comprenant les ordres des individus par dimension

TraitementRegle <- function(RegleCourante, NumGroupeCourant, NumeroGroupes,
                            Order, GroupesActuels){
  # Indices dans tableRegle
  dimCol <- 1 ; lambdaCol <- 2 ; MinusCol <- 3
  SplitCol <- 4 ; Splitp1Col <- 5 ; MaxCol <- 6
  
  # Indices dans Groupes  
  NumCol <- 1 ; LongueurCol <- 2 ; LambdaCol <- 3
  NbSortantCol <- 4 ; SplitONCol <- 5 ; ParentCol <- 6
  
  # Recuperation des elements concernes par la regle a partir de l'ordre de la 
  # dimension
  Elements <- Order[[RegleCourante[1]]][PlusPetitPaquet(RegleCourante)]
  # Recuperation des numeros de groupes concernes et du nombre d'element par groupe
  # voulant sortir
  b <- NumeroGroupes[Elements]
  a <- unique(b)
  Table <- numeric(length(a))
  for(i in 1:length(a)){
    Table[i] <- sum(b==a[i])
  }
  NumerosGroupesConcernes <- a
  Split <- FALSE
  
  for(i in 1:length(NumerosGroupesConcernes)){
    GroupesActuels[NumerosGroupesConcernes[i],][NbSortantCol] <-  Table[i]
    if(Table[i]< GroupesActuels[NumerosGroupesConcernes[i],][LongueurCol]){
      GroupesActuels[NumerosGroupesConcernes[i],][SplitONCol] <- 1
      NumGroupeCourant <- NumGroupeCourant+1
      NouveauGroupe <- c(NumGroupeCourant,Table[i],
                         RegleCourante[lambdaCol],
                         0,0,
                         GroupesActuels[NumerosGroupesConcernes[i],][NumCol])
      GroupesActuels <- rbind(GroupesActuels, NouveauGroupe)
      NumeroGroupes[Elements[NumeroGroupes[Elements] ==
                               NumerosGroupesConcernes[i]]] <- NumGroupeCourant
      Split <- TRUE
      # Actualisation de la longueur des groupes :
      GroupesActuels[NumerosGroupesConcernes[i],][LongueurCol] <- 
        GroupesActuels[NumerosGroupesConcernes[i],][LongueurCol]- Table[i]
    }
  }
  return(list(Split = Split, GroupesActuels = GroupesActuels, 
              NumGroupeCourant = NumGroupeCourant, NumeroGroupes = NumeroGroupes))
}

#####################################################-
# Transformation en matrice merge type hclust$merge #
#####################################################-

# Parametres :
# - n : le nombre d'individus
# - Groupes : une liste des etats successifs des groupes a chaque activation
#             d'une regle impliquant un split
# - NumeroGroupes : une matrice 2xn contenant pour chaque element
#                   son numero de groupe

CreationMatriceMerge <- function(n,Groupes, NumeroGroupes){
  MatriceIndGroupe <- matrix(1:n, nrow=n, ncol=2,byrow = FALSE)
  
  Final <- Groupes[[length(Groupes)]][n:1,]
  
  # Indices dans Groupes/Final  
  NumCol <- 1 ; LongueurCol <- 2 ; LambdaCol <- 3
  NbSortantCol <- 4 ; SplitONCol <- 5 ; ParentCol <- 6
  
  NumeroGroupeMerge <- n+1
  MatriceMerge <- matrix(0,(length(Final[,1])-1),3)
  
  # Premiere ligne : l'individu
  # deuxieme ligne : le groupe dans lequel est l'individu
  NumerosGroupes <- as.matrix(rbind(1:n,NumeroGroupes), ncol = n)
  
  NumerosGroupes <- NumerosGroupes[,order(NumerosGroupes[2,])]
  
  for(i in 1:(length(Final[,1])-1)){
    IndividuDansLeGroupeFin <- NumerosGroupes[1,Final[i,1]]
    DansLeGroupeFin <- MatriceIndGroupe[IndividuDansLeGroupeFin,2]
    Parent <- Final[i,6]
    DansleGroupeParent <- MatriceIndGroupe[NumerosGroupes[1,Final[i,6]],2]
    MatriceIndGroupe[NumerosGroupes[1,Final[i,6]],2] <- NumeroGroupeMerge
    MatriceIndGroupe[NumerosGroupes[1,Final[i,1]],2] <- NumeroGroupeMerge
    MatriceIndGroupe[MatriceIndGroupe[,2]==DansleGroupeParent,2] <- NumeroGroupeMerge
    MatriceMerge[i,] <- c(NumeroGroupeMerge,
                          DansLeGroupeFin,
                          DansleGroupeParent)
    NumeroGroupeMerge <- NumeroGroupeMerge+1
  }
  
  MatriceMerge <- MatriceMerge[,2:3]
  
  ind.Inf.np1 <- MatriceMerge < n+1
  MatriceMerge[ind.Inf.np1]  <- -MatriceMerge[ind.Inf.np1]
  MatriceMerge[!ind.Inf.np1] <- MatriceMerge[!ind.Inf.np1] -n
  
  Height <- Final[,LambdaCol][Final[,LambdaCol]!=0]
  ##Height <- Final[-1,LambdaCol] # Modification : 30/11/17
  return(list(MatriceMerge= MatriceMerge, Height = Height))
}

# Sortie :
# - MatriceMerge : la matrice de type hclust$merge associee a l'arbre
# - Height : les lambdas, qui definissent les hauteurs des fusions

######################################################################-
# Recuperation de l'ordre des individus a partir de la liste Groupes #
######################################################################-

# A peut-?tre inclure dans la fonction qui fabrique la matrice Merge

# Parametres :
# - n : le nombres d'individus
# - Groupes : une liste des etats successifs des groupes a chaque activation
#             d'une regle impliquant un split
# - NumeroGroupes : une matrice 2xn contenant pour chaque element
#                   son numero de groupe

OrdreIndividus <- function(n,Groupes, NumeroGroupes){
  GroupesFinaux <- Groupes[[length(Groupes)]]
  # print(GroupesFinaux)
  
  # Indices dans Groupes/Final  
  NumCol <- 1 ; LongueurCol <- 2 ; LambdaCol <- 3
  NbSortantCol <- 4 ; SplitONCol <- 5 ; ParentCol <- 6
  
  #   NumerosGroupes <- as.matrix(rbind(1:n,NumeroGroupes), ncol = n)
  #   NumerosGroupes <- NumerosGroupes[,order(NumerosGroupes[2,])]
  
  Order <- numeric(2)
  Order[1] <- 1
  for(i in 2:n){
    IndiceParent <- match(GroupesFinaux[i,][ParentCol],Order)
    Order1 <- Order[(1:IndiceParent)]
    Order2 <- Order[((IndiceParent+1):length(Order))]
    Order <- c(Order1, GroupesFinaux[i,][NumCol], Order2)
  }
  
  Order <- Order[Order!=0]
  for(i in 1:n){
    Order[i] <- match(Order[i],NumeroGroupes)
    # Order[i] <- NumerosGroupes[1,Order[i]]
  }
  
  return(Order)
}

# Sortie :
# - Order : l'ordre des individus, tel qu'il n'y aura pas de coupure dans hclust


#######################################################-
###### Fonctions de similarite (calcul laplacien) #####
#######################################################-

# Found online, cbind for a list sparse vectors
# https://stackoverflow.com/questions/8843700/creating-sparse-matrix-from-a-list-of-sparse-vectors
# sv.cbind <- function (...) {
#   input <- lapply( list(...), as, "dsparseVector" )
#   thelength <- unique(sapply(input,length))
#   stopifnot( length(thelength)==1 )
#   return( sparseMatrix( 
#     x=unlist(lapply(input,slot,"x")), 
#     i=unlist(lapply(input,slot,"i")), 
#     p=c(0,cumsum(sapply(input,function(x){length(x@x)}))),
#     dims=c(thelength,length(input))
#   ) )
# }

s <- function(x1, x2, alpha=1) {
  # exp(- alpha * norm(as.matrix(x1-x2), type="F"))
  exp(-alpha*sum((x1-x2)^2))
  # -alpha*sum((x1-x2)^2)
}

# make.sparse.similarity <- function(my.data, similarity){
#   # my.data = t(data_matrix)
#   # similarity = s
#   library(Matrix)
#   N <- nrow(my.data)
#   res = lapply(1:N, FUN = function(x){
#     # x = 2
#     row_x = rep(0, N)
#     for(j in x:N) {
#       row_x[j] =  similarity(my.data[x,], my.data[j,])
#     }
#     return(as(row_x, "sparseVector"))
#   })
#   S = do.call(sv.cbind, res)
#   return(S)
# }

# make.sparse.similarity <- function(my.data, similarity, cutoff = 0.2){
#   N <- nrow(my.data)
#   S <- matrix(rep(NA,N^2), ncol=N)
#   res = lapply(1:N, FUN = function(x){
#     row_x = rep(NA, N)
#     for(j in x:N) {
#       row_x[j] =  similarity(my.data[x,], my.data[j,])
#     }
#     return(row_x)
#   })
#   S = do.call(rbind, res)
#   S[lower.tri(S, diag = FALSE)] <- t(S)[lower.tri(S, diag = FALSE)]
#   S[abs(S)<cutoff] <- 0
#   S = as(S, "sparseMatrix")
#   return(S)
# }
# 
# make.sparse.similarity <- function(my.data, similarity, alpha, cutoff = 0.2){
#   N <- nrow(my.data)
#   # S <- matrix(rep(0,N^2), ncol=N)
#   res = lapply(1:N, FUN = function(x){
#     row_x = rep(0, N)
#     for(j in x:N) {
#       row_x[j] =  similarity(my.data[x,], my.data[j,], alpha)
#     }
#     return(row_x)
#   })
#   S = do.call(rbind, res)
#   # S[lower.tri(S, diag = FALSE)] <- t(S)[lower.tri(S, diag = FALSE)]
#   S[abs(S)<cutoff] <- 0
#   S = as(S, "sparseMatrix")
#   return(S)
# }

make.similarity <- function(my.data, similarity, alpha, lower = FALSE, sparse = TRUE, cutoff = NULL) {
  N <- nrow(my.data)
  S <- matrix(rep(0,N^2), ncol=N)
  alpha_new = alpha/sd(unlist(my.data))^2
  res = lapply(1:N, FUN = function(x){
    row_x = rep(0, N)
    for(j in x:N) {
      row_x[j] =  similarity(my.data[x,], my.data[j,], alpha_new)
    }
    return(row_x)
  })
  S = do.call(rbind, res)
  if(!is.null(cutoff)) S[abs(S)<cutoff] <- 0
  if(lower) S[lower.tri(S, diag = FALSE)] <- t(S)[lower.tri(S, diag = FALSE)]
  if(sparse) S = as(S, "sparseMatrix")
  return(S)
}

##############################-
# Fonction qui regroupe tout # 
##############################-

# Parametres :
# - Data : matrice des donnees, avec en ligne ce que l'on veut clusteriser
# - verb : verbose
# - esp : si verbose = TRUE, affiche toutes les esp regles, le nombre de groupes 
# et la regle regardee.
# - spectral : utilisation d'une SVD et des vecteurs propres pour faire le clustering
# - k spectral : le nombre de vecteurs propres a retenir pour le clustering si spectral
# - data_type : a enlever et a mettre en auto
# - gamma : vecteur applique a chacune des dimensions (colonnes) pour le clustering
# - parallel. :
# - nu : pour la SVD, combien de vecteurs propres sont a calculer. PAr defaut le nombre de
# colonnes de DATA (cense etre plus petit que n...)
#


# FA <- DataSets[[1]]
# FAC <- AgregationArbres(list(FA),1000,1)

AgregationArbres.wrapper <- function(Data, verb = TRUE, esp = 1000, 
                                     spectral = TRUE,
                                     k_spectral = NULL,
                                     data_type = "counts",  # data_type = counts ou continuous
                                     gamma = NULL,# puisque spectral TRUE par defaut... 
                                     # cutoff = NULL,
                                     # random_svd = TRUE,
                                     use_eigs = TRUE,
                                     # use_sym = FALSE,
                                     parallel. = TRUE,
                                     alpha = 1/2,
                                     nu = NULL
                                     # nv = 0
                                     
){
  library("plyr") ; library("fusedanova"); library("parallel") ; library(rARPACK); library(PRIMME);
  # library(rsvd)
  Data = t(data_matrix)
  # cutoff = NULL
  p = ncol(Data)
  n = nrow(Data)
  # alpha = 1/2
  
  if(is.null(nu)) nu = n
  if(is.null(k_spectral)){k_spectral = ceiling(p/10)}
  
  if(spectral){
    # Data = t(data_matrix)
    S <- make.similarity(Data, s, alpha, lower = TRUE, sparse = FALSE, cutoff = NULL)  
    eS = eigen(S, only.values = TRUE)$values
    to_delete = sum(abs(eS)<=1e-16)
    if(use_eigs){
      # if(abs(det(diag(rowSums(S))-S))<=1e-5 || det(diag(rowSums(S))-S)==Inf||det(diag(rowSums(S))-S)==-Inf){
      #   message("det(diag(rowSums(S))-S) is either too small or too high, trying something else")
      #   # S <- chol(S)
      #   # Z = eigs_sym((crossprod(diag(rowSums(S))-S)), k = k_spectral, which = "SM")$vectors
      #   # Z = eigs_sym(crossprod(diag(rowSums(S))-S), k = k_spectral, which = "SM")$vectors
      #   # S = crossprod(t(Data))
      #   # k_spectral = k_spectral+10
      #   Z = eigs_sym((diag(rowSums(S))-S),  k = p, which = "SM")$vectors
      # }else{
        # S <- make.sparse.similarity(Data, s, alpha, cutoff = 0)
        # Z = eigs_sym(diag(rowSums(S+t(S))-S), k = k_spectral, which = "SM")$vectors
        Z = rARPACK::eigs_sym((diag(rowSums(S))-S), k = p, which = "SM")$vectors
      # compute p vectors instead of k_spectral because some of the smallest might not converge. 
      # Z = PRIMME::eigs_sym((diag(rowSums(S))-S), NEig = p, which = "SA")$vectors 
        # PRIMME = plus long que rARPACK et donne de moins bons resultats
      # }
        # if(to_delete>0) Z = Z[,-((ncol(Z)-to_delete+1):ncol(Z))]
      Z = Z[, (ncol(Z)-k_spectral+1):ncol(Z)]
    }
    else{
      L_SVD = svd(diag(rowSums(S))-S, nu = nu, nv = 0)$u
      # Z = L_SVD$u[,(ncol(L_SVD$u)-k_spectral+1):ncol(L_SVD$u)]
      # L_SVD = PRIMME::svds(diag(rowSums(S))-S, p, which = "S")$u # random... 
      # if(to_delete>0) L_SVD = L_SVD[,-((ncol(L_SVD)-to_delete+1):ncol(L_SVD))]
      Z = L_SVD[,(ncol(L_SVD)-k_spectral+1):ncol(L_SVD)]
    }
    
    rm(S)
    
    gamma = rep(0.5, k_spectral) 
    p_old = p
    Data_old = Data
    p = ncol(Z)
    Data = Z
  }
  
  if(parallel.){
    DataSets <- mclapply(1:p, function(i) {
      ## fusedanova(Data[,i])
      fusedanova(Data[,i], weights = "laplace", gamma = gamma[i], standardize = TRUE)
    }, mc.cores=4, mc.preschedule = TRUE)
  }else{
    DataSets = lapply(1:p, FUN = function(i){
      fusedanova(Data[,i], weights = "laplace", gamma = gamma[i], standardize = TRUE)
    } )
  }
  
  
  ListeTablesRegles <- llply(DataSets, RecuperationRegles, n = n)
  
  dimCol <- 1 ; lambdaCol <- 2 ;
  
  for(i in 1:p){
    ListeTablesRegles[[i]][,dimCol] <- i
  }
  
  Order <- lapply(1:p, FUN=function(i) DataSets[[i]]@result[[1]]$order)
  
  ReunionRegles <- do.call("rbind",ListeTablesRegles)
  ReunionRegles <- ReunionRegles[order(ReunionRegles[,lambdaCol], 
                                       decreasing = TRUE),]
  
  NumeroGroupes <- rep(1,n)
  NumGroupeCourant <- 1
  Groupes <- list()
  NumCol <- 1 ; LongueurCol <- 2 ; LambdaCol <- 3 ;
  NbSortant <- 4 ; SplitONCol <- 5 ; ParentCol <- 6 ;
  Groupes[[1]] <- matrix(c(1,n,0,0,0,0), 1,6)
  IndiceDeRegle <- 1
  ReglesActives <- matrix(0, ncol = 6, nrow = n-1)
  compteur <- 1
  
  while(NumGroupeCourant < n && IndiceDeRegle <= nrow(ReunionRegles)){
    if(IndiceDeRegle %% esp == 0 && verb == TRUE){
      print(paste0(IndiceDeRegle, " / ", nrow(ReunionRegles)))
      print(NumGroupeCourant)
    }
    GroupesActuels <- Groupes[[length(Groupes)]]
    RegleCourante <- ReunionRegles[IndiceDeRegle,]
    Resultat <- TraitementRegle(RegleCourante, NumGroupeCourant, NumeroGroupes,
                                Order, GroupesActuels)
    if(Resultat$Split){
      Groupes[[length(Groupes)+1]] <- Resultat$GroupesActuels
      NumGroupeCourant <- Resultat$NumGroupeCourant
      NumeroGroupes <- Resultat$NumeroGroupes
      ReglesActives[compteur,] <- RegleCourante
      compteur <- compteur+1
    }
    IndiceDeRegle <- IndiceDeRegle+1
  }
  
  if(verb){
    print(paste0(IndiceDeRegle, " / ", nrow(ReunionRegles)))
    print(NumGroupeCourant)
  }
  
  ReglesActives <- ReglesActives[ReglesActives[,dimCol]!=0,]
  Ordre <- OrdreIndividus(n, Groupes, NumeroGroupes)
  out <- CreationMatriceMerge(n, Groupes,NumeroGroupes)
  
  MatriceMerge  <- out$MatriceMerge
  Height <- out$Height
  
  Cluster <- list(merge = MatriceMerge, height = Height, order = Ordre)
  Cluster <- unclass(Cluster)
  class(Cluster) <- "hclust"
  # plot(Cluster)
  return(list(Cluster = Cluster,Regles = ReglesActives,
              Orders = Order, ReunionRegles = ReunionRegles))
}

####################################################-
#### Recuperation hierarchie sur une dimension #####
####################################################-

# DataVect = vecteur des donnees de la dimension concernee
# resFA : un resultat de fusedANOVA pour la dimension concernee

Recup_hier = function(DataVect, resFA, verb = FALSE, esp = 1000){
  p = 1
  n = length(DataVect)
  
  DataSets <- list(resFA)
  
  # Partie du code identique a celle du wrapper.
  ListeTablesRegles <- llply(DataSets, RecuperationRegles, n = n)
  
  dimCol <- 1 ; lambdaCol <- 2 ; MinusCol <- 3
  SplitCol <- 4 ; Splitp1Col <- 5 ; MaxCol <- 6
  
  for(i in 1:p){ListeTablesRegles[[i]][,dimCol] <- i}
  
  Order <- lapply(1:p, FUN=function(i) DataSets[[i]]@result[[1]]$order)
  
  ReunionRegles <- do.call("rbind",ListeTablesRegles)
  ReunionRegles <- ReunionRegles[order(ReunionRegles[,lambdaCol], 
                                       decreasing = TRUE),]
  NumeroGroupes <- rep(1,n)
  NumGroupeCourant <- 1
  Groupes <- list()
  NumCol <- 1 ; LongueurCol <- 2 ; LambdaCol <- 3 ;
  NbSortant <- 4 ; SplitONCol <- 5 ; ParentCol <- 6 ;
  Groupes[[1]] <- matrix(c(1,n,0,0,0,0), 1,6)
  IndiceDeRegle <- 1
  ReglesActives <- matrix(0, ncol = 6, nrow = n-1)
  compteur <- 1
  
  while(NumGroupeCourant < n && IndiceDeRegle <= nrow(ReunionRegles)){
    if(IndiceDeRegle %% esp == 0 && verb == TRUE){
      print(paste0(IndiceDeRegle, " / ", nrow(ReunionRegles)))
      print(NumGroupeCourant)
    }
    GroupesActuels <- Groupes[[length(Groupes)]]
    RegleCourante <- ReunionRegles[IndiceDeRegle,]
    Resultat <- TraitementRegle(RegleCourante, NumGroupeCourant, NumeroGroupes,
                                Order, GroupesActuels)
    if(Resultat$Split){
      Groupes[[length(Groupes)+1]] <- Resultat$GroupesActuels
      NumGroupeCourant <- Resultat$NumGroupeCourant
      NumeroGroupes <- Resultat$NumeroGroupes
      ReglesActives[compteur,] <- RegleCourante
      compteur <- compteur+1
    }
    IndiceDeRegle <- IndiceDeRegle+1
  }
  
  if(verb){
    print(paste0(IndiceDeRegle, " / ", nrow(ReunionRegles)))
    print(NumGroupeCourant)
  }
  
  ReglesActives <- ReglesActives[ReglesActives[,dimCol]!=0,]
  Ordre <- OrdreIndividus(n, Groupes, NumeroGroupes)
  out <- CreationMatriceMerge(n, Groupes,NumeroGroupes)
  
  MatriceMerge  <- out$MatriceMerge
  Height <- out$Height
  
  Cluster <- list(merge = MatriceMerge, height = Height, order = Ordre)
  Cluster <- unclass(Cluster)
  class(Cluster) <- "hclust"
  
  cut_tree <- cutree(Cluster,1:n)
  
  BIC_dim_all =  lapply(1:n, FUN = function(x){
    BIC(LogLikelihoodDim(DataVect, group =cut_tree[,x]), nbGroupes = x, n = n, p=1)
  })
  
  # Logiquement le nombre de groupes correspond au cut dans l'arbre. 
  
  best_cut = which.min(BIC_dim_all)  
  best_group = cut_tree[, best_cut]
  
  return(list(Cluster = Cluster, BIC_dim_all = BIC_dim_all, best_cut = best_cut, 
              best_group = best_group, Regles = ReglesActives,
              Orders = Order, ReunionRegles = ReunionRegles))
} 

#####################################################-
#### Fonction de selection du gamma base sur le BIC #
#####################################################-

# gamma_to_try : vecteur de valeurs de gamma
# Data : une table n*p ? ou juste une dimension ?
Gamma_BIC = function(Data, gamma_to_try, parallel = TRUE, verb = FALSE){
  # Data = t(data_matrix)
  p = ncol(Data)
  n = nrow(Data)
  
  best_gamma = rep(NA, p)
  for(i in 1:p){
    if(verb) print(i)
    if(parallel){
      test_gamma = unlist(mclapply(1:length(gamma_to_try), FUN = function(x){
        res <- Recup_hier(Data[,i], fusedanova(Data[,i], weights = "laplace", gamma = gamma_to_try[x], standardize = TRUE))
        return(res$BIC_dim_all[which.min(res$BIC_dim_all)])
      }, mc.cores=4, mc.preschedule = TRUE))
    }else{
      test_gamma = unlist(lapply(1:length(gamma_to_try), FUN = function(x){
        res <- Recup_hier(Data[,i], fusedanova(Data[,i], weights = "laplace", gamma = gamma_to_try[x], standardize = TRUE))
        return(res$BIC_dim_all[which.min(res$BIC_dim_all)])
      }))
    }
    best_gamma[i] = gamma_to_try[which.min(test_gamma)]
  }
  
  return(best_gamma)
}



############################################################-
#### Fonctions likelihood et BIC pour validation gamma #####
############################################################-

# recuperees du stage de M2

LogLikelihoodDim <- function(data, group){
  Y <- data
  n  <- length(Y)
  p  <- 1
  K  <- length(unique(group))
  nk <- tabulate(group)
  
  betak <- rowsum(Y, group)/nk
  RSS <- sum((Y - betak[group])^2)
  
  sigma2 <- RSS/(n-K)
  
  return(-.5 * (n*p*log(2*pi) + n*sum(log(sigma2)) + (n-K)))
  
}

BIC <- function(loglikelihood, nbGroupes, n, p){
  return(-2*loglikelihood+log(n)*p*nbGroupes)
}


hclust.one.dim <- function(n, fa1) {
  
  Order <- list(fa1@result[[1]]$order)
  fusion_rules <- RecuperationRegles(n, fa1)

  NumeroGroupes <- rep(1,n)
  NumGroupeCourant <- 1
  Groupes <- list()
  Groupes[[1]] <- matrix(c(1,n,0,0,0,0), 1,6)
  IndiceDeRegle <- 1
  ReglesActives <- matrix(0, ncol = 6, nrow = n-1)
  compteur <- 1
  
  while(NumGroupeCourant < n && IndiceDeRegle <= nrow(fusion_rules)){
    GroupesActuels <- Groupes[[length(Groupes)]]
    RegleCourante <- fusion_rules[IndiceDeRegle,]
    Resultat <- TraitementRegle(RegleCourante, NumGroupeCourant, NumeroGroupes,
                                Order, GroupesActuels)
    if(Resultat$Split){
      Groupes[[length(Groupes)+1]] <- Resultat$GroupesActuels
      NumGroupeCourant <- Resultat$NumGroupeCourant
      NumeroGroupes <- Resultat$NumeroGroupes
      ReglesActives[compteur,] <- RegleCourante
      compteur <- compteur+1
    }
    IndiceDeRegle <- IndiceDeRegle+1
  }
  
  order <- OrdreIndividus(n, Groupes, NumeroGroupes)
  out <- CreationMatriceMerge(n, Groupes,NumeroGroupes)
 
  height <- sort(unique(fa1@result[[1]]$table$lambda))[-1]
  # order <- fa1@result[[1]]$order
  Cluster <- list(merge = out$MatriceMerge, height = height, order = order)
  class(Cluster) <- "hclust"
  Cluster
}
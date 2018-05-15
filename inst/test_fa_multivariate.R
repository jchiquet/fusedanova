rm(list=ls())
library(fusedanova)
library(tidyverse)
data(iris)

fa <- fusedanova(iris %>% select(-Species))

# print(fa$path)
# par(mfrow=c(1,3))
# plot(fa, main="fused-ANOVA new")
# plot(AIC(fa), main="fused-ANOVA new - AIC", type="l")
# plot(BIC(fa), main="fused-ANOVA new - BIC", type="l")
# par()

# source("inst/agregation_to_hc.R")
# n <- nrow(iris)
# mat_element = rbind(1:n, fa$currentGroup)
# 
# # Matrice des parents et des enfants                                                                                                                                                                               
# mat_groupes = rbind(1:n, fa$groupsParent, fa$groupsChildCurrent, fa$groupsIRule)
# 
# # Remise dans l'ordre : faire correspondre la matrice des groupes a la matrice des elements                                                                                                                        
# mat_element2 = mat_element[,order(mat_element[2,])]
# 
# mat_aide = rbind(mat_element2, mat_groupes)
# mat_aide2 = mat_aide[-3,]
# 
# l_element = 1 ; l_groupeAct = 2; l_parent = 3 ; l_enfant = 4 ;
# 
# # Ceux qui ont pour enfants 0 n'ont en fait pas d'enfant                                                                                                                                                           
# mat_aide2[l_enfant, which(mat_aide2[l_enfant,]==0)] <- NA
# mat_aide3 = mat_aide2[,order(mat_aide2[l_groupeAct,], decreasing = TRUE)]
# matrice_aide2 = (mat_aide3[-5,])
# MatriceMerge = CreationMatriceMerge(n, matrice_aide2)
# 
# # Ordre :                                                                                                                                                                                                          
# Order = OrdreIndividus(matrice_aide2)

# # Height :                                                                                                                                                                                                         
# # dans le cas ou les individus sont degroupes artificiellement, completer les hauteurs de coupures.                                                                                                                
# Height <- rev(lambdaRules)
# if(length(Height)!=(n-1)){Height = c(rep(0, n-1-length(Height)), Height)}
# 
# Cluster <- list(merge = MatriceMerge, height = Height, order = Order)
# Cluster <- unclass(Cluster)
# class(Cluster) <- "hclust"


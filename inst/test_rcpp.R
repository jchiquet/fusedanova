library(fusedanova)
source("inst/functions_Audrey.R")
data(aves)

fa1 <- fusedanova(aves$weight, aves$family)
rule_fa1 <- RecuperationRegles(40, fa1)
fa1_path <- fa1@result[[1]]$table[fa1@result[[1]]$table$lambda != 0, ]
fa1_path <- fa1_path[order(fa1_path$lambda, decreasing = TRUE), ]

fa2 <- fusedanova2(aves$weight, aves$family)
fa2_path <- fa2$path
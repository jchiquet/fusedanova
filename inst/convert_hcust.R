source("inst/functions_Audrey.R")
library(fusedanova)
data(aves)
n <- length(aves$family)
fa <- fusedanova(aves$weight, aves$family)

fa2 <- fusedanova2(aves$weight, aves$family)

regles <- RecuperationRegles(n,fa)

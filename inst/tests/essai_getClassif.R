library(fusedanova)
library(dplyr)
data(aves)
fa0 <- fusedanova(aves$weight, aves$family) 
fa <- fusedanova2(aves$weight, aves$family) 

path <- fa$path %>% filter(iup != idown) %>% arrange(desc(lambda))

if (is.null(lambdas))
  lambdas <- unique(path$lambda)

K <- nlevels(fa$cl0)

## fils
fils <- CalculeFils(as.integer(path$idown), K)	
classif <- vCalculeClassif(lambdas,
                           path$lambda,
                           as.integer(path$idown),
                           as.integer(path$iup),
                           as.integer(fils), as.integer(fa$cl0), K)
res <- matrix(classif, nrow = K)

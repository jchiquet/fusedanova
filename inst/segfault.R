library(fusedanova)
data(aves)
out <- fusedanova(aves$weight, aves$family, gamma = 1)

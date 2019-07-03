rm(list = ls())
library(Rmergetrees)
library(fusedanova)
library(aricode) 
library(rsvd)
library(parallel)
library(dplyr)
library(ggplot2)
library(viridis)
theme_set(theme_bw())
source("util.R")

## separability
## 10:  everyboby is working

## SIMULATION PARAMETERS
n_grp <- 20
grp_size <- 50
separability <- 1
n0 <- 30
n1 <- 50
n_grp_per_data <- 2
n_sim <- 10

# SCENARIO 1: CHAQUE JEUX DE DONNÉES N1 PORTE LA MÊME INFORMATION
# 
# DÉFAVORABLE/AUCUN INTÉRÊT POUR MERGETREES
# 
res_scenario1 <- mclapply(1:n_sim, function(sim_label) {
  dataSets  <- scenario1(n0, n1, n_grp, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size) 
  oneRun(sim_label, dataSets, reference, k = 5)
}, mc.cores = 10)

res_plot1 <- Reduce("rbind", res_scenario1) %>% 
  group_by(method, nb_grp) %>% 
  summarize(NID = mean(NID), NID_sd = sd(NID)) %>% 
  filter(method %in% c("One", "MT", "MTW", "DC", "AD", "rSMT", "SMTW"))

p1 <- ggplot(res_plot1,  aes(nb_grp, NID, color = method)) + 
  geom_line() + geom_point(aes(shape = method)) + xlim(0, 100) + scale_color_viridis(discrete = TRUE)

# SCENARIO 2: INFORMATION PORTÉE PAR PLUSIEURS JEUX DE DONNÉES
# 
# FAVORABLE À MERGETREES
# 
res_scenario2 <- mclapply(1:n_sim, function(sim_label) {
  dataSets  <- scenario2(n0, n1, n_grp, n_grp_per_data, grp_size, separability)
  reference <- rep(1:n_grp, each = grp_size) 
  oneRun(sim_label, dataSets, reference)
}, mc.cores = 10)

res_plot2 <- Reduce("rbind", res_scenario2) %>% 
  group_by(method, nb_grp) %>% 
  summarize(NID = mean(NID), NID_sd = sd(NID)) %>% 
  filter(method %in% c("One", "MT", "MTW", "DC", "AD", "rSMT", "SMTW"))

p2 <- ggplot(res_plot2, aes(nb_grp, NID, color = method)) + 
  geom_line() + geom_point(aes(shape = method)) + xlim(0, 100) + scale_color_viridis(discrete = TRUE)


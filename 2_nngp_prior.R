rm(list = ls())
source("./src/NNMatrix.R")
source("./src/functions.R")
load("./data/lgcp_data_full.RData")
getwd()
library(foreach)
library(doParallel)
library(parallel)
library(INLA)

## Build nearest neighbor matrix from spNNGP ##
M = 10

NN.matrix_full <- NNMatrix(coords = loc_n_N, n.neighbors = M,
                      n.omp.threads = 2, search.type = "cb", ord = 1:nrow(loc_n_N))

## set estimates for phi and deltasq from cross-validation "spConjNNGP" ##
rho = sqrt((win$xrange[2]-win$xrange[1])^2 + (win$yrange[2]-win$yrange[1])^2)
phi = 2/rho
N_sum = n+N
ind_x <-c(c(rep(2:M, times = 1:(M - 1)), rep(((M + 1) : N_sum), each = M)), 1:N_sum)
ind_y <- c(c(t(NN.matrix_full$NN_ind))[which(c(t(NN.matrix_full$NN_ind)) > 0)], 1:N_sum)

## obtain A and D using C and N(i)
AD <- getADstan(neardist = NN.matrix_full$NN_dist,  neardistM = NN.matrix_full$NN_distM, 
                N = N_sum, M = M, phi = phi) 
D <- AD[M + 1, ]

X_star <- 
  sparseMatrix(i = ind_x, j = ind_y, 
               x = c(-na.omit(as.vector(AD[-(M + 1), ])), rep(1, N_sum))) / sqrt(D) 

Gamma_precision <- crossprod(X_star, X_star)

save(list = ls(all.names = TRUE), 
     file = "./data/nngp_prior_full.RData", envir = .GlobalEnv)








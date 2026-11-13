getwd()
rm(list = ls())
source('./src/data_pattern_generator.R')

library(FNN)
library(dplyr)

#----------------------- Get lgcp data ---------------------------
# Data pattern generation
domain_max = 1
win <- owin(c(0,domain_max), c(0,domain_max))
spatstat.options(npixel = 256)
# spatstat.options(npixel = 21*21)
y0 = x0 = seq(win$xrange[1], win$xrange[2], length = spatstat.options()$npixel)
d = 2

gridcov_func = function(x, y) {
  
  z <- (
    0.9 * sin(2 * pi * (1.3 * x + 1.7 * y)) +
      0.7 * cos(2 * pi * (2.1 * x - 1.4 * y)) +
      0.5 * sin(2 * pi * (3.2 * x * y)) +
      0.4 * cos(2 * pi * (4.0 * (x - 0.5)^2 + 3.5 * (y - 0.5)^2)) +
      0.3 * sin(2 * pi * (2.8 * x^2 - 2.3 * y^2))
  )
  
  # Normalize smoothly to [-1,1] without saturation
  z / 2.8
}


gridcov <- outer(y0, x0, function(y,x) gridcov_func(x,y))


beta_0 = 7.3
beta_1 = 2.5
nu = 0.5
kappa = 1/4
sigma2 = 0.5
Data_pattern = Data_pattern_generator(beta_0, beta_1, nu, kappa,sigma2, random_seed = 1, max_edge = 0.5)

xy.c_observed = Data_pattern$xy_c_observed
xy.c_test = Data_pattern$xy_c_test
N = nrow(xy.c_observed)

lg.s.c = Data_pattern$lg_s_c
lambda_values <- attr(lg.s.c, 'Lambda')
lambda_values$v = log(lambda_values$v)
log_intensity_value = lambda_values$v

# par(mfrow = c(1,2), oma = c(0, 0, 0, 4))
# image(lambda_values, col = viridis(256), main = "")
# par(mfrow = c(1,1), oma = c(0, 0, 0, 0))



n_test = as.integer(0.2*N)
margin = 0.05

min_dist = 0.04

sampling_x_min = win$xrange[1] + margin
sampling_x_max = win$xrange[2] - margin
sampling_y_min = win$yrange[1] + margin
sampling_y_max = win$yrange[2] - margin

if (sampling_x_min >= sampling_x_max || sampling_y_min >= sampling_y_max) {
  stop("Margin is too large, resulting in an invalid sampling area.")
}

sampling_win <- owin(c(sampling_x_min, sampling_x_max), c(sampling_y_min, sampling_y_max))
cat("Attempting to generate", n_test, "points with minimum distance", min_dist, "...\n")
ppp_test <- rSSI(r = min_dist, n = n_test, win = sampling_win, giveup = 1000, verbose = FALSE)

n_generated <- npoints(ppp_test)
if (n_generated < n_test) {
  warning(paste("Could only generate", n_generated, "points out of the requested",
                n_test, "due to the minimum distance constraint (min_dist =", min_dist, ").",
                "Consider reducing n_test or min_dist if more points are needed."))
} else {
  cat("Successfully generated", n_generated, "points.\n")
}


loc_test <- coords(ppp_test)
loc_test = as.matrix(loc_test)
print(head(loc_test))
cat("Actual number of points generated:", nrow(loc_test), "\n")

# plot(win, main = paste("Generated Test Locations (rSSI, min_dist =", min_dist, ")"))
# plot(sampling_win, border = "blue", add = TRUE)
# points(loc_test, pch = 16, col = "red")
# legend("topright", legend = c("Original Domain", "Sampling Area", "Generated Points"),
#        lty = c(1, 1, NA), pch = c(NA, NA, 16), col = c("black", "blue", "red"))

if (n_generated > 1) {
  min_observed_dist <- min(nndist(ppp_test))
  cat("Minimum distance observed between generated points:", round(min_observed_dist, 4), "\n")
  if(min_observed_dist < min_dist - 1e-9) {
    warning("Observed minimum distance is slightly less than specified min_dist!")
  }
}
n_test = length(loc_test[,1])
#------------------------------ mesh --------------------------------------
pl.dom = cbind(c(0,domain_max,domain_max,0), c(0,0,domain_max,domain_max))
bnd1 = Polygon(pl.dom, hole = FALSE)
mesh_full = inla.mesh.2d(loc_test, pl.dom, boundary = bnd1, max.edge = 0.4)
# plot(mesh_full)
# points(xy.c_observed)

N = nrow(xy.c_observed)
n = mesh_full$n


mesh_vertices <- mesh_full$loc[ ,1:2]
nn_result <- get.knnx(data = mesh_vertices,
                      query = loc_test,
                      k = 1)
mesh_ind_loc_test <- nn_result$nn.index[, 1]

nn_result_test <- get.knnx(data = mesh_vertices,
                           query = xy.c_test,
                           k = 1)
mesh_ind_test <- nn_result_test$nn.index[,1]
mesh_test_count <- tabulate(mesh_ind_test, nbins = nrow(mesh_vertices))

# plot
# par(mar = c(0, 0, 2, 0))
# plot(mesh_full, asp = 1, main = '')
# points(mesh_full$loc[,1:2], pch = 20, col = 'black')
# points(loc_test, pch = 20, col = 'red')
# points(mesh_full$loc[mesh_ind_loc_test,1:2], col = 'blue')

#------------------------------- Dual mesh -------------------------------
alpha_weight = DualMesh_parameter(mesh_full)[[1]]
Dualmesh_loc = DualMesh_parameter(mesh_full)[[2]]

loc_n_N = rbind(mesh_vertices, xy.c_observed)

#------------------------------- covariate ----------------------------------
X = cbind(rep(1, N), gridcov[Reduce('cbind', nearest.pixel(xy.c_observed[,1], xy.c_observed[,2],im(gridcov, x0, y0)))])
X_tilde = cbind(rep(1, n), gridcov[Reduce('cbind', nearest.pixel(mesh_vertices[,1], mesh_vertices[,2],im(gridcov, x0, y0)))])
X_test = cbind(rep(1, n_test), gridcov[Reduce('cbind', nearest.pixel(loc_test[,1], loc_test[,2],im(gridcov, x0, y0)))])


# A_tilde and A
zero_block1 = matrix(0, nrow = n, ncol = N)
zero_block2 = matrix(0, nrow = N, ncol = n)
A_tilde = as(cbind(diag(n), zero_block1), "dgCMatrix")
A = as(cbind(zero_block2, diag(N)), "dgCMatrix")

rm(zero_block1, zero_block2)
rm(x0,y0,sampling_win,sampling_x_max,sampling_x_min,sampling_y_max,sampling_y_min,min_observed_dist,min_dist,margin,ppp_test)
rm(n_generated)

getwd()
save(list = ls(all.names = TRUE), file = "./data/lgcp_data_full.RData", envir = .GlobalEnv)


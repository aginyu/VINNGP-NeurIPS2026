rm(list = ls())
source("./optimization_hyper_diag.R")
getwd()
load("./data/nngp_prior_full.RData")
load("./data/lgcp_data_full.RData")

nN = n+N
beta_ini = matrix(c(5,5), 2, 1)
mu_ini = matrix(0, nN, 1)
desired_rank = 15
decomposition_result <- decompose_gamma(Gamma_precision, rank_r = desired_rank)
Gamma_base_diag_vec <- decomposition_result$Gamma_base_diag_vec 
L_base_factor       <- decomposition_result$L_base_factor      
r_actual            <- decomposition_result$r_actual



if (!exists("N_sum")) N_sum <- n + N 
if (is.null(L_base_factor)) {
  stop("L_base_factor is NULL after decomposition. Check decompose_gamma.")
} else if (ncol(L_base_factor) != r_actual) {
  warning("ncol(L_base_factor) does not match r_actual. Using ncol.")
  r_prior_dim <- ncol(L_base_factor) # Use actual dimension from matrix
} else {
  r_prior_dim <- r_actual # Use the value from decomposition result
}
Sigma_ini= list(diag = 1 / Gamma_base_diag_vec,             
                   W = matrix(0.0, nrow=N_sum, ncol=r_prior_dim), 
                   L = L_base_factor)     
sigma2_ini = 2.5
a = 0.01
b = 0.01


Rprof("my_profile.out", interval = 0.01) 

Gamma_base_precision = Gamma_precision
VINNGP_cg_lowrank2 = VIGGNP(A, A_tilde, 
                                       X, X_tilde, alpha_weight,
                                       mu_ini, Sigma_ini, beta_ini, sigma2_ini, 
                                       Gamma_precision, # Base precision (unscaled)
                                       Gamma_base_diag_vec, L_base_factor, # Decomposed base prior parts
                                       a, b, # IG prior parameters for sigma2
                                       maxiter = 20, tol = c(5,0.05), sigma2_tol = 1e-4, # Outer loop tolerances
                                       inner_maxiter = 5, # Max iterations for inner loops
                                       inner_tol_mu = 1e-6,  # Inner tolerance for mu update
                                       inner_tol_sigma = 1e-6, # Inner tolerance for Sigma update
                                       inner_tol_beta = 1e-6, # Inner tolerance for beta update
                                       mu_update_method = "cg1",
                                       cg_tol = 1e-6, cg_maxiter = 5,cg_step = c(0.01, 0.01),fix_sigma = TRUE)


Rprof(NULL)

summaryRprof("my_profile.out")  
readLines("my_profile.out")

# result
beta_vinngp = VINNGP_cg_lowrank2$beta_hist[[length(VINNGP_cg_lowrank2$beta_hist)]]
mu_vinngp_n = VINNGP_cg_lowrank2$mu[1:n]
mu_vinngp_pred = VINNGP_cg_lowrank2$mu[mesh_ind_loc_test]
mu_vinngp_N = VINNGP_cg_lowrank2$mu[(n+1):(n+N)]
sigma2_vinngp = VINNGP_cg_lowrank2$sigma2

Precision.approx = 1.0/sigma2_vinngp * Gamma_precision
diag(Precision.approx)[1:n] = diag(Precision.approx)[1:n] + exp_WeightedVec_lowrank(alpha_weight, VINNGP_cg_lowrank2$mu, VINNGP_cg_lowrank2$Sigma_lowrank$diag, beta_vinngp, X_tilde)
Sigma_vinngp_nn_diag = diag(inla.qinv(Precision.approx))[1:n]

llk_vinngp_full = as.numeric(sum(as.vector(X%*%beta_vinngp + mu_vinngp_N)) - sum(alpha_weight*exp(as.vector(X_tilde%*%beta_vinngp + mu_vinngp_n))))
llk_vinngp_full_n = as.numeric(sum(mesh_test_count * as.vector(X_tilde%*%beta_vinngp + mu_vinngp_n)) - sum(alpha_weight*exp(as.vector(X_tilde%*%beta_vinngp + mu_vinngp_n))))
print(c(llk_vinngp_full,llk_vinngp_full_n,VINNGP_cg_lowrank2[["running_time"]]))

# #-------------------------------------------- INLA -------------------------------------------------
spatstat.options(npixel = nrow(gridcov))
 y0 = x0 = seq(win$xrange[1], win$xrange[2], length = spatstat.options()$npixel)
spde_sim = inla.spde2.pcmatern(
  mesh = mesh_full,
  alpha = 1.5,
  prior.range = c(0.05, 0.01), ### P(practic.range<0.05)=0.01
  prior.sigma = c(1, 0.01) ### P(sigma>1)=0.01
)
covariate_sim = gridcov[Reduce('cbind', nearest.pixel(
  c(mesh_vertices[,1], xy.c_observed[,1]), c(mesh_vertices[,2], xy.c_observed[,2]),
  im(gridcov, x0, y0)))]
nv = n
n.c = N
y.pp.c = rep(0:1, c(nv, n.c))
e.pp.c = c(alpha_weight, rep(0, n.c))
lmat.c = inla.spde.make.A(mesh_full, xy.c_observed)
imat = Diagonal(nv, rep(1, nv))
A.pp.c = rbind(imat, lmat.c)

# data stack
stk.pp = inla.stack(data = list(y=y.pp.c, e=e.pp.c),
                    A = list(1, A.pp.c), tag = 'pp.c',
                    effects = list(list(b0=1, covariate = covariate_sim),
                                   list(i=1:nv)))
# INLA model training
pp.c.res = inla(y ~ 0 + b0 + covariate + f(i, model=spde_sim),
                family = 'poisson', data = inla.stack.data(stk.pp),
                control.predictor = list(A=inla.stack.A(stk.pp), compute = TRUE),
                E = inla.stack.data(stk.pp)$e)

dnames <- unlist(lapply(strsplit(dimnames(pp.c.res$summary.linear.predictor)[[1]], ".", fixed = T), function(x){x[1]}))
pp.c.res$fitted.values <- pp.c.res$summary.linear.predictor[dnames == "APredictor", ][ , "mean"]
pp.c.res$fitted.values.sd <- pp.c.res$summary.linear.predictor[dnames == "APredictor", ][ , "sd"]
llk_INLA_full = as.numeric(sum(pp.c.res$fitted.values[(n+1):(N+n)])) - as.numeric(sum(alpha_weight*exp(pp.c.res$fitted.values[1:n])))
llk_INLA_full_n = as.numeric(sum(mesh_test_count * pp.c.res$fitted.values[1:n])) - as.numeric(sum(alpha_weight*exp(pp.c.res$fitted.values[1:n])))

#--------------------------------------------- VIFRK -----------------------------------------------
library(scampr)
# model fit
dat_VIFRK = data.frame("x" = c(xy.c_observed[ ,1], mesh_vertices[ ,1]),
                       "y" = c(xy.c_observed[ ,2], mesh_vertices[ ,2]),
                       "pres" = rep(1:0, c(N, n)),
                       "covariates" = gridcov[Reduce('cbind', nearest.pixel(
                         c(xy.c_observed[,1], mesh_vertices[,1]), c(xy.c_observed[,2], mesh_vertices[,2]),
                         im(gridcov, x0, y0)))],
                       "quad.size" = c(rep(0, N), alpha_weight))
# VIFRK model training
VIFRK = scampr(pres~covariates,
               data = dat_VIFRK,
               coord.names = c("x", "y"),
               quad.weights.name = "quad.size",
               model.type = "PO", 
               basis.functions = simple_basis(10, dat_VIFRK))

llk_VIFRK_full = as.numeric(sum(VIFRK$fitted.values[1:N])) - as.numeric(sum(alpha_weight*exp(VIFRK$fitted.values[(N+1): length(VIFRK$fitted.values)])))
llk_VIFRK_full_n = as.numeric(sum(mesh_test_count * VIFRK$fitted.values[(N+1): length(VIFRK$fitted.values)])) - as.numeric(sum(alpha_weight*exp(VIFRK$fitted.values[(N+1): length(VIFRK$fitted.values)])))
numbasis = (length(VIFRK$coefficients)-length(X[1,]))/2
Z = VIFRK$tmb.data$Z_PO_quad
VIFRK_var = Z %*% (as.vector(exp(VIFRK$coefficients[(numbasis+length(X[1,])+1): (2*numbasis+length(X[1,]))])) * t(Z))

#--------------------------------------------- kppm -----------------------------------------------
library(spatstat.geom)
library(spatstat.model)

obs_ppp <- ppp(
  x = xy.c_observed[,1],
  y = xy.c_observed[,2],
  window = win
)
cov_im <- im(gridcov, x0, y0)

mesh_ppp <- ppp(
  x = mesh_vertices[,1],
  y = mesh_vertices[,2],
  window = win
)

time_mincon <- system.time({
  fit_mincon <- kppm(obs_ppp, trend = ~covariate, clusters = "LGCP",
                   covariates = list(covariate = gridcov_func),
                   statistic = "pcf", model = "matern", nu = 0.5, method = 'mincon')
})[["elapsed"]]
print(time_mincon)
pred_mincon_N <- predict(fit_mincon, locations = obs_ppp,type = "lambda")

lambda_mincon_N <- as.numeric(pred_mincon_N)
eta_mincon_N <- log(lambda_mincon_N)
pred_mincon_n <- predict(fit_mincon, locations = mesh_ppp,type = "intensity")
lambda_mincon_n <- as.numeric(pred_mincon_n)
eta_mincon_n <- log(lambda_mincon_n)
llk_mincon_full = as.numeric(sum(eta_mincon_N)) - as.numeric(sum(alpha_weight*lambda_mincon_n))
llk_mincon_full_n = as.numeric(sum(mesh_test_count * eta_mincon_n)) - as.numeric(sum(alpha_weight*lambda_mincon_n))

time_quasi <- system.time({
  fit_quasi <- kppm(obs_ppp, trend = ~covariate, clusters = "LGCP",
                   covariates = list(covariate = gridcov_func),
                   statistic = "pcf", model = "matern", nu = 0.5,
                   method = "mincon", improve.type = 'quasi') 
})[["elapsed"]]
print(time_quasi)

pred_quasi_N <- predict(fit_quasi, locations = obs_ppp, type = "intensity")
lambda_quasi_N <- as.numeric(pred_quasi_N)
eta_quasi_N <- log(lambda_quasi_N)

pred_quasi_n <- predict(fit_quasi, locations = mesh_ppp, type = "intensity")
lambda_quasi_n <- as.numeric(pred_quasi_n)
eta_quasi_n <- log(lambda_quasi_n)

llk_quasi_full <- as.numeric(sum(eta_quasi_N)) - as.numeric(sum(alpha_weight * lambda_quasi_n))
llk_quasi_full_n <- as.numeric(sum(mesh_test_count * eta_quasi_n)) - as.numeric(sum(alpha_weight * lambda_quasi_n))

# ---------------------------------------------- result ---------------------------------------------
logintensity_mean_full_df = data.frame(
  "VINNGP" = as.vector(X_tilde%*%beta_vinngp + mu_vinngp_n),
  "INLA" = as.vector(pp.c.res$fitted.values[1:n]),
  "VIFRK" = as.vector(VIFRK$fitted.values[(N+1):(N+n)])
)
# 1. Extract original Lambda (DO NOT overwrite)
lambda_values <- attr(lg.s.c, "Lambda")

# 2. Grid intensity (true ¦Ë(s))
true_intensity_grid <- lambda_values$v

# 3. Build spatstat image (reuse x0, y0)
lambda_im <- im(true_intensity_grid, x0, y0)

# 4. Evaluate at mesh locations (n points)
true_intensity_mesh <- interp.im(
  lambda_im,
  x = mesh_vertices[,1],
  y = mesh_vertices[,2]
)

# 5. Also get log-intensity
true_log_intensity_mesh <- log(true_intensity_mesh)

l2_INLA = 0
l2_VINNGP = sqrt(sum(alpha_weight * (true_intensity_mesh - exp(logintensity_mean_full_df$VINNGP))^2))/sum(alpha_weight)
l2_INLA = sqrt(sum(alpha_weight * (true_intensity_mesh - exp(logintensity_mean_full_df$INLA))^2))/sum(alpha_weight)
l2_VIFRK = sqrt(sum(alpha_weight * (true_intensity_mesh - exp(logintensity_mean_full_df$VIFRK))^2))/sum(alpha_weight)
l2_mincon = sqrt(sum(alpha_weight * (true_intensity_mesh - lambda_mincon_n)^2))/sum(alpha_weight)
l2_quasi = sqrt(sum(alpha_weight * (true_intensity_mesh - lambda_quasi_n)^2))/sum(alpha_weight)
l2_BLANK = sqrt(sum(alpha_weight * (true_intensity_mesh - 0)^2))/sum(alpha_weight)

l2_log_VINNGP = sqrt(sum(alpha_weight * (true_log_intensity_mesh - logintensity_mean_full_df$VINNGP)^2))/sum(alpha_weight)
l2_log_INLA = sqrt(sum(alpha_weight * (true_log_intensity_mesh - logintensity_mean_full_df$INLA)^2))/sum(alpha_weight)
l2_log_VIFRK = sqrt(sum(alpha_weight * (true_log_intensity_mesh - logintensity_mean_full_df$VIFRK)^2))/sum(alpha_weight)
l2_log_mincon = sqrt(sum(alpha_weight * (true_log_intensity_mesh - eta_mincon_n)^2))/sum(alpha_weight)
l2_log_quasi = sqrt(sum(alpha_weight * (true_log_intensity_mesh - eta_quasi_n)^2))/sum(alpha_weight)
l2_log_BLANK = sqrt(sum(alpha_weight * (true_log_intensity_mesh - 0)^2))/sum(alpha_weight)

logintensity_sd_full_df = data.frame(
  "VINNGP" = sqrt(as.vector(Sigma_vinngp_nn_diag)),
  "INLA" = as.vector(pp.c.res$fitted.values.sd[1:n]),
  "VIFRK" = sqrt(as.vector(diag(VIFRK_var)))
)

table_full = data.frame(
  Method = c("VoGCAM", "INLA", "VIFRK",'MINCON','QUASI'),
  llk_N = c(llk_vinngp_full, llk_INLA_full, llk_VIFRK_full,llk_mincon_full,llk_quasi_full),
  llk_n = c(llk_vinngp_full_n, llk_INLA_full_n, llk_VIFRK_full_n, llk_mincon_full_n,llk_quasi_full_n),
  l2 = c(l2_VINNGP,l2_INLA,l2_VIFRK,l2_mincon,l2_quasi),
  l2_log = c(l2_log_VINNGP,l2_log_INLA,l2_log_VIFRK,l2_log_mincon,l2_log_quasi),
  Time = c(VINNGP_cg_lowrank2[["running_time"]], pp.c.res$cpu.used[4], VIFRK$cpu[1],time_mincon,time_quasi)
)

table_full

save(logintensity_mean_full_df, logintensity_sd_full_df, VIFRK, pp.c.res, VINNGP_cg_lowrank2, mesh_full, table_full,eta_mincon_n,lambda_mincon_n, eta_quasi_n, lambda_quasi_n, true_log_intensity_mesh,true_intensity_mesh, file = "./data/comparsion_fullcase.RData", envir = .GlobalEnv)

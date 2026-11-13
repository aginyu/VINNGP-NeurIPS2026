
#rm(list = ls())

source("./src/dual_mesh_loc.R")

# Load necessary libraries
library(RSpectra) # For decompose_gamma if used later
library(stats)   # For optim, lgamma

# N: the number of data pattern
# n: the number of basis functions (the number of vertices in the mesh's triangles)
# p: the number of data points used for approximate integration (the number of polygons in the dual mesh)
# N_sum: N+n

# @ A (N*N_sum) -> Assumed implicit via indices n, N
# @ A_tidle (n*N_sum) -> Assumed implicit via indices n, N
# @ X (N*m): a matrix composed of the values of m covariates for N data patterns.
# @ X_tilde (n*m): a matrix composed of the values of m covariates for p integration data points.
# @ alpha_weight (n*1): Weights used for approximate integration, obtained when constructing the dual mesh.
# @ beta_c (m*1): the coefficients corresponding to m covariates.
# @ mu (N_sum*1): posterior mean
# @ Sigma_k_list: List representing posterior covariance (low-rank format)
# @ sigma2: Variance hyperparameter
# @ Gamma_base_precision: The base precision matrix (unscaled by sigma2)
# @ Gamma_base_diag_vec: Diagonal of Gamma_base_precision
# @ L_base_factor: Low-rank factor L from Gamma_base_precision = Diag + L L^T
# @ a, b: Hyperparameters for Inverse-Gamma prior on sigma2: IG(a, b)

##################################################################################

# trace of sparse matrix
sparsemat_tr = function(mat){
  return(sum(diag(mat)))
}

# log(det)
log_det = function(mat){
  return(Matrix::determinant(mat, logarithm = TRUE)$modulus[1])
}

Sigma_process = function(Sigma_k_list, sigma2_k,n,update.sigma.only = FALSE){
  L_base_factor = Sigma_k_list$L
  r_prior = ncol(L_base_factor)
  I_r = Matrix::Diagonal(r_prior)
  D_inv_vec = Sigma_k_list$diag
  
  # Ensure correct dimensions
  if (length(D_inv_vec) != nrow(L_base_factor)) stop("Dimension mismatch: Sigma_k_list$diag")
  if (sigma2_k <= 0) stop("sigma2_k must be positive in exp_WeightedVec_lowrank")
  
  # Calculate diagonal of Sigma_k = D_inv - D_inv L' (I+M)^{-1} L'^T D_inv
  # M = L'^T D_inv L' = (1/sigma2_k) * L_base^T D_inv L_base
  if(update.sigma.only == FALSE){
    Lt <- t(L_base_factor)
    M0 <- Lt %*% (D_inv_vec * L_base_factor)
  } else {M0 = Sigma_k_list$M0}
  M = (1/sigma2_k) * M0
  
  # Middle_inv = solve(I + M)
  Middle_inv = tryCatch({
    Matrix::solve(I_r + M)
  }, error = function(e){
    warning("Matrix solve failed for (I+M) in exp_WeightedVec_lowrank. Check M's condition number. Using pseudo-inverse.", call. = FALSE)
    # Fallback or error handling, e.g., using MASS::ginv or stopping
    MASS::ginv(as.matrix(I_r + M)) # Ensure MASS is loaded or handle differently
  })
  
  # L' = L_base / sqrt(sigma2_k)
  L_prime = (1/sqrt(sigma2_k)) * L_base_factor
  
  # Precompute V = L' %*% Middle_inv
  V = L_prime %*% Middle_inv  # N_sum x r
  
  # Calculate diag(Sigma_k)_i = D_inv_i - D_inv_i^2 * sum_j V[i,j] * L_prime[i,j]
  # More efficiently: D_inv_i - D_inv_i^2 * rowSums(V * L_prime)
  Correction_term = D_inv_vec^2 * rowSums(V * L_prime)
  Sigma_k_diag_elems = D_inv_vec - Correction_term
  
  # Check for negative variances (shouldn't happen if stable)
  if(any(Sigma_k_diag_elems <= 0)){
    warning("Non-positive diagonal elements calculated for Sigma_k in exp_WeightedVec_lowrank. Clamping.")
    Sigma_k_diag_elems = pmax(Sigma_k_diag_elems, .Machine$double.eps)
  }
  
  Sigma_nn_diag = Sigma_k_diag_elems[1:n]
  return (list(diag = D_inv_vec,
               W = Sigma_k_list$W,
               L = L_base_factor,
               M0 = M0, 
               M = M, 
               Middle_inv = Middle_inv, 
               Sigma_nn_diag = Sigma_nn_diag,
               Sigma_k_diag_elems = Sigma_k_diag_elems))
}




exp_WeightedVec_lowrank = function(alpha_weight, mu, Sigma, beta_c, X_tilde){    # Current sigma^2 value
  n = nrow(X_tilde)
  N_sum = length(mu)
  mu_n = mu[1:n]
  Sigma_n = Sigma[1:n]
  
  if (!is.matrix(beta_c)) beta_c = matrix(beta_c, ncol=1)
  exp_arg = X_tilde %*% beta_c + mu_n + 0.5 * Sigma_n
  
  return(as.vector(alpha_weight) * exp(exp_arg))
}

# Calculates trace(Gamma_base * Sigma_k) using the low-rank representation of Sigma_k
# Gamma_base = G_diag + L L^T
# Sigma_k = D_inv - W L'^T D_inv
# trace(Gamma_base Sigma_k) = trace(G_diag D_inv) + trace(L L^T D_inv)
#                           - trace(G_diag W L'^T D_inv) - trace(L L^T W L'^T D_inv)
# Uses cyclic property: trace(A B) = trace(B A)
trace_Gamma_Sigma_lowrank <- function(Gamma_base_diag_vec, Sigma) {
  
  return(sum(Gamma_base_diag_vec  * Sigma))
}

# Log determinant of Sigma_k
# Based on det(Sigma) = det(D_inv) / det(I+M)
# where M = (1/sigma2) * L_base^T * D_inv * L_base
log_det_Sigma_lowrank <- function(Sigma) {
  
  return(sum(log(Sigma)))
}

F_ELBO_lowrank_mu = function(X_tilde, alpha_weight, mu, Sigma_k_list,beta_c, sigma2_k,
                             Gamma_base_precision # Base precision matrix (Needed for mu^T G mu)
                             ) {               # IG prior parameters for sigma2
  n = nrow(X_tilde)
  N_sum = ncol(X_tilde)
  N = N_sum - n
  mu_N = mu[(n+1):N_sum]
  
  # Likelihood related terms (using updated exp_WeightedVec_lowrank)
  # Pass L_base_factor and sigma2_k to the function
  exp_w_vec = exp_WeightedVec_lowrank(alpha_weight, mu, Sigma_k_list, beta_c, X_tilde)
  # Check for non-finite values from exp_WeightedVec_lowrank
  if (any(!is.finite(exp_w_vec))) {
    warning("Non-finite values in exp_w_vec within ELBO calculation.", call.=FALSE)
    exp_w_vec[!is.finite(exp_w_vec)] <- 0 # Or handle more robustly
  }
  log_likelihood_term = sum(mu_N) - sum(exp_w_vec)
  
  # Prior terms for mu (scaled by sigma2)
  # E_q[ log p(mu | sigma2) ] ~ -1/(2*sigma2) * E_q[mu^T Gamma_base mu] - 1/2 log det(sigma2 * inv(Gamma_base))
  # E_q[mu^T Gamma_base mu] = mu_k^T Gamma_base mu_k + tr(Gamma_base Sigma_k)
  muT_G_mu = as.numeric(Matrix::crossprod(mu, Gamma_base_precision %*% mu))
  
  log_prior_mu_energy_term = -1 / (2 * sigma2_k) * muT_G_mu
  
  
  ELBO = log_likelihood_term + log_prior_mu_energy_term
  
  if (!is.finite(ELBO)) {
    warning("Non-finite ELBO computed.", call.=FALSE)
    # You might want to inspect the components here if debugging
    # print(paste("log_lik:", log_likelihood_term))
    # print(paste("log_prior_E:", log_prior_mu_energy_term))
    # print(paste("log_prior_D:", log_prior_mu_det_term))
    # print(paste("entropy:", entropy_term))
    # print(paste("log_p_sig2:", log_prior_sigma2_term))
    return(NA)
  }
  
  return(as.numeric(ELBO))
}


# ELBO Calculation based on formula (12) interpretation and standard VI
# Assumes prior precision is sigma^{-2} * Gamma_base
# Assumes prior on sigma2 is IG(a, b)
F_ELBO_lowrank = function(X, X_tilde, alpha_weight, mu, Sigma, beta_c_k, sigma2_k, Gamma_base_diag_vec,
                                 log_det_gamma_base, a, b) {               # IG prior parameters for sigma2
  n = nrow(X_tilde)
  N = nrow(X)
  N_sum = n + N
  mu_N = mu[(n+1):N_sum]
  
  # Likelihood related terms (using updated exp_WeightedVec_lowrank)
  # Pass L_base_factor and sigma2_k to the function
  exp_w_vec = exp_WeightedVec_lowrank(alpha_weight, mu, Sigma, beta_c_k, X_tilde)
  # Check for non-finite values from exp_WeightedVec_lowrank
  if (any(!is.finite(exp_w_vec))) {
    warning("Non-finite values in exp_w_vec within ELBO calculation.", call.=FALSE)
    exp_w_vec[!is.finite(exp_w_vec)] <- 0 # Or handle more robustly
  }
  log_likelihood_term = sum(X %*% beta_c_k) + sum(mu_N) - sum(exp_w_vec)
  
  # Prior terms for mu (scaled by sigma2)
  # E_q[ log p(mu | sigma2) ] ~ -1/(2*sigma2) * E_q[mu^T Gamma_base mu] - 1/2 log det(sigma2 * inv(Gamma_base))
  # E_q[mu^T Gamma_base mu] = mu_k^T Gamma_base mu_k + tr(Gamma_base Sigma_k)
  muT_G_mu = as.numeric(Matrix::crossprod(mu, Gamma_base_precision %*% mu))
  tr_G_Sigma = trace_Gamma_Sigma_lowrank(Gamma_base_diag_vec,Sigma)
  
  # Check for NA from trace calculation
  if (is.na(tr_G_Sigma)) {
    warning("trace_Gamma_Sigma calculation failed in ELBO.", call.=FALSE)
    return(NA)
  }
  
  prior_energy_expect = muT_G_mu + tr_G_Sigma
  log_prior_mu_energy_term = -1 / (2 * sigma2_k) * prior_energy_expect
  
  # Log determinant of prior precision term: - 1/2 log det(sigma^2 * Gamma_base^{-1})
  # = - 1/2 * (N_sum * log(sigma2_k) - log_det_gamma_base)
  log_prior_mu_det_term = -0.5 * (N_sum * log(sigma2_k) - log_det_gamma_base)
  
  # Entropy term H(q) = 0.5 * log det(Sigma_k) + 0.5 * N_sum * (1 + log(2*pi))
  # We only need the log det part for optimization: 0.5 * log det(Sigma_k)
  log_det_sigma_val = log_det_Sigma_lowrank(Sigma)
  
  if (is.na(log_det_sigma_val)) {
    warning("Log_det computation failed in ELBO.", call.=FALSE)
    return(NA)
  }
  # Standard entropy term
  entropy_term = 0.5 * log_det_sigma_val
  
  # Prior terms for sigma2 ~ IG(a, b) -> log p(sigma2)
  # log p(sigma2) = - (a + 1) log sigma2 - b / sigma2 + a log b - lgamma(a)
  log_prior_sigma2_term = - (a + 1) * log(sigma2_k) - b / sigma2_k
  # Constants from IG prior (often dropped, but include for completeness)
  # ig_const = a * log(b) - lgamma(a) # Use b from input
  
  # Constant term n/2 from Eq 12 (origin unclear, but include per formula)
  const_term_n_half = n / 2.0
  
  # Assemble ELBO based on F = E_q[log p(Y, mu, sigma2)] - E_q[log q(mu)]
  # F = E_q[log p(Y|mu,beta)] + E_q[log p(mu|sigma2)] + log p(sigma2) - E_q[log q(mu)]
  # F = log_likelihood_term + (log_prior_mu_energy_term + log_prior_mu_det_term)
  #     + log_prior_sigma2_term + ig_const + entropy_term + N_sum/2*(1+log(2pi))
  # Drop constants N_sum/2*(1+log(2pi)) and ig_const if desired. Add n/2 per Eq 12.
  ELBO = log_likelihood_term + log_prior_mu_energy_term + log_prior_mu_det_term +
    entropy_term + log_prior_sigma2_term + const_term_n_half # + ig_const
  
  if (!is.finite(ELBO)) {
    warning("Non-finite ELBO computed.", call.=FALSE)
    # You might want to inspect the components here if debugging
    # print(paste("log_lik:", log_likelihood_term))
    # print(paste("log_prior_E:", log_prior_mu_energy_term))
    # print(paste("log_prior_D:", log_prior_mu_det_term))
    # print(paste("entropy:", entropy_term))
    # print(paste("log_p_sig2:", log_prior_sigma2_term))
    return(NA)
  }
  
  return(as.numeric(ELBO))
}


# Helper function for L-BFGS: Calculate Negative ELBO w.r.t. mu (ignoring constants w.r.t mu)
NegELBO_mu_lowrank <- function(mu, Sigma_k_list, beta_c_k, sigma2_k, Gamma_base_precision,
                               X_tilde, alpha_weight, n, N, tA_eN_vec, L_base_factor) {
  
  mu = as.matrix(mu)
  exp_w_vec = exp_WeightedVec_lowrank(alpha_weight, mu, Sigma_k_list, beta_c_k, X_tilde)
  mu_N = mu[(n + 1):(n + N), 1, drop=FALSE]
  
  likelihood_term1 = sum(mu_N)
  likelihood_term3 = -sum(exp_w_vec)
  # Prior term depends on sigma2_k
  quardratic_pen = -1 / (2 * sigma2_k) * as.numeric(Matrix::crossprod(mu, Gamma_base_precision %*% mu))
  
  neg_elbo = -(likelihood_term1 + likelihood_term3 + quardratic_pen)
  return(neg_elbo)
}

# Helper function for L-BFGS: Calculate Gradient of Negative ELBO w.r.t. mu
grad_NegELBO_mu_lowrank <- function(mu, Sigma_k_list, beta_c_k, sigma2_k, Gamma_base_precision,
                                    X_tilde, alpha_weight, n, N, tA_eN_vec, L_base_factor) {
  mu = as.matrix(mu)
  exp_w_vec = exp_WeightedVec_lowrank(alpha_weight, mu, Sigma_k_list, beta_c_k, X_tilde) # Pass args
  
  # Calculate f_mu (gradient of positive ELBO)
  f_mu_lik_term2 = matrix(c(exp_w_vec, rep(0, N)), ncol=1) # N_sum x 1
  f_mu_lik_term3 = matrix(tA_eN_vec, ncol=1) # N_sum x 1 (precomputed t(A)%*%e_N)
  
  # Prior gradient term (gradient of -1/(2*sigma2)*mu^T*G*mu is -1/sigma2*G*mu)
  f_mu_prior_term = -(1 / sigma2_k) * (Gamma_base_precision %*% mu)
  
  f_mu = f_mu_lik_term3 - f_mu_lik_term2 + f_mu_prior_term # grad(Lik) + grad(log p(mu))
  
  # Gradient for optim (minimization) is -f_mu (gradient of negative ELBO)
  grad = -f_mu
  return(as.vector(grad))
}

# Helper function for CG: Calculate Hessian-vector product (Hessian_of_NegELBO %*% v)
# Hessian of NegELBO = - Hessian of ELBO
# Hessian of ELBO = H_lik + H_prior
# H_lik = block_diag(exp_w_vec, 0)
# H_prior = -1/sigma2 * Gamma_base_precision
# HessianVecProd calculates (-H_lik - H_prior) %*% v = -(H_lik + H_prior) %*% v
HessianVecProd_mu_lowrank <- function(v, mu_k, Sigma_k_list, beta_c_k, sigma2_k, Gamma_base_precision,
                                      X_tilde, alpha_weight, n, N, L_base_factor) {
  v = as.matrix(v)
  N_sum = n + N
  if(nrow(v) != N_sum) stop("Vector v has incorrect dimension for HessianVecProd")
  
  exp_w_vec_k = exp_WeightedVec_lowrank(alpha_weight, mu_k, Sigma_k_list, beta_c_k, X_tilde) 
  
  # Term from prior Hessian: (-H_prior) %*% v = (1/sigma2 * Gamma_base_precision) %*% v
  term1 = (1 / sigma2_k) * (Gamma_base_precision %*% v)
  
  # Term from likelihood Hessian: (-H_lik) %*% v = -(block_diag(...) %*% v)
  v_n = v[1:n, 1, drop=FALSE]
  prod_part = as.vector(exp_w_vec_k) * v_n
  term2 = matrix(c(prod_part, rep(0, N)), ncol=1) # This is H_lik %*% v
  
  # Result is term1 + term2 because H_negELBO = H_prior_scaled + H_lik
  result = term1 + term2
  return(result) # Return as matrix
}

build_Hessian_sparse <- function(mu_k, Sigma_k_list, beta_c_k, sigma2_k,
                                 Gamma_base_precision,
                                 X_tilde, alpha_weight,
                                 n, N, L_base_factor) {
  
  # Total dimension of the parameter vector
  N_sum <- n + N
  
  # ------------------------------------------------------------
  # 1. Prior Hessian term
  #    H_prior = (1 / sigma2_k) * Gamma_base_precision
  #    Gamma_base_precision is assumed to be a sparse matrix
  # ------------------------------------------------------------
  H_prior <- (1 / sigma2_k) * Gamma_base_precision
  
  # ------------------------------------------------------------
  # 2. Likelihood Hessian term (diagonal structure)
  #    Only the first n components have non-zero curvature
  # ------------------------------------------------------------
  exp_w_vec_k <- exp_WeightedVec_lowrank(
    alpha_weight, mu_k, Sigma_k_list, beta_c_k, X_tilde
  )
  
  # Construct full diagonal vector of length (n + N)
  # First n entries correspond to likelihood curvature
  # Remaining N entries are zero
  diag_full <- c(as.vector(exp_w_vec_k), rep(0, N))
  
  # Sparse diagonal matrix
  H_lik <- Diagonal(x = diag_full)
  
  # ------------------------------------------------------------
  # 3. Combine prior and likelihood terms
  # ------------------------------------------------------------
  H <- H_prior + H_lik
  
  return(H)
}

# ============================================================
# Solve A x = v where A = Gamma_precision = X_star' X_star
# using sparse triangular solves instead of generic CG
#
# Since X_star is triangular:
#   A^{-1} v = (X_star' X_star)^{-1} v
#            = X_star^{-1} (X_star')^{-1} v
#
# More efficiently:
#   1. Solve X_star' y = v
#   2. Solve X_star  x = y
#
# This is exact, highly efficient, and exploits sparsity.
# ============================================================

Solve_GammaPrecision <- function(v, X_star,X_star_t) {
  
  y <- solve(X_star_t, v)
  x <- solve(X_star, y)
  
  return(x)
}

Mul_GammaPrecision <- function(v, grad.f_mu) {
  
  return(grad.f_mu %*% v)
}


# Basic Conjugate Gradient Solver (copied from thought process, seems reasonable)
cg_solver <- function(Ax_func, b, x0, tol = 1e-5, maxiter = 100, ...) {
  x <- as.matrix(x0); b <- as.matrix(b)
  r <- b - Ax_func(x, ...); p <- r
  rsold <- sum(r^2)
  if (sqrt(rsold) < tol) return(list(x = x, iterations = 0, converged = TRUE))
  
  iter = 0; converged = FALSE
  for (i in 1:maxiter) {
    iter = i; Ap <- Ax_func(p, ...); pAp = sum(p * Ap)
    if (!is.finite(pAp) || pAp == 0) {
      warning("CG: pAp is zero or non-finite. Stopping."); break }

    if (!is.finite(pAp) || pAp <= tol * norm(p, "2") * norm(Ap, "2")) { # 添加�?查，tol 可以�? sqrt(.Machine$double.eps)
      warning("CG: Matrix may not be positive definite or pAp is numerically zero/negative/non-finite. Stopping.")
      converged = FALSE # 标记未收�?
      break # �?出循�?
    }
    # --- 结束�?�? ---
    alpha <- rsold / pAp
    x <- x + alpha * p; r <- r - alpha * Ap
    rsnew <- sum(r^2)
    if (sqrt(rsnew) < tol) { converged = TRUE; break }
    if (!is.finite(rsold) || rsold == 0) {
      warning("CG: rsold is zero or non-finite. Stopping."); break }
    p <- r + (rsnew / rsold) * p; rsold <- rsnew
  }
  if (i == maxiter && !converged) warning(paste("CG did not converge within", maxiter, "iterations."))
  return(list(x = x, iterations = iter, converged = converged))
}

pcg_solver <- function(H_operator,b,X_star,X_star_t,x0,M_inv_func = NULL,tol = 1e-6,maxiter = 1000,relative = FALSE,verbose = FALSE) {
  # Ensure column vectors
  b <- as.numeric(b)
  #n <- length(b)
  x <- as.numeric(x0)
  # Initial residual
  r <- b - H_operator(x)
  z <- Solve_GammaPrecision(r,X_star,X_star_t)
  p <- z
  rz_old <- sum(r * z)
  # Initial norm
  b_norm <- sqrt(sum(b^2))
  if (b_norm == 0) b_norm <- 1
  residual_norm <- sqrt(sum(r^2))
  # Stopping threshold
  threshold <- if (relative) tol * b_norm else tol
  if (residual_norm < threshold) {
    return(list(
      x = x,
      iterations = 0,
      converged = TRUE,
      residual_norm = residual_norm
    ))
  }
  converged <- FALSE
  iter <- 0
  for (i in 1:maxiter) {
    iter <- i
    Ap <- H_operator(p)
    pAp <- sum(p * Ap)
  # SPD / numerical stability check
    if (pAp <= .Machine$double.eps) {
      warning("PCG stopped: non-positive or numerically unstable p^T A p.")
      break
    }
    alpha <- rz_old / pAp
    # if (!is.finite(alpha)) {
    #   warning("PCG stopped: non-finite step size.")
    #   break
    # }
    # Update solution
    x <- x + alpha * p
    # Update residual
    r <- r - alpha * Ap
    residual_norm <- sqrt(sum(r^2))
    if (verbose) {
      cat(sprintf("Iter %d: residual = %.6e\n", i, residual_norm))
    }
    # Convergence check
    if (residual_norm < threshold) {
      converged <- TRUE
      break
    }
    # Apply preconditioner
    z <- Solve_GammaPrecision(r,X_star,X_star_t)
    rz_new <- sum(r * z)
    # if (!is.finite(rz_new) || abs(rz_old) <= .Machine$double.eps) {
    #   warning("PCG stopped: numerical breakdown in residual inner product.")
    #   break
    # }
    beta <- rz_new / rz_old
    # if (!is.finite(beta)) {
    #   warning("PCG stopped: non-finite beta.")
    #   break
    # }
    # Update search direction
    p <- z + beta * p
    
    rz_old <- rz_new
  }
  if (!converged && iter == maxiter) {
    warning(sprintf("PCG did not converge within %d iterations.", maxiter))
  }
  return(list(
    x = x,
    iterations = iter,
    converged = converged,
    residual_norm = residual_norm
  ))
}

# Main algorithm incorporating sigma2 with inner loops
VIGGNP = function(A, A_tilde, # Note: A, A_tilde seem unused
                            X, X_tilde, alpha_weight,
                            mu0, Sigma0_list, beta_c0, sigma2_0, # Initial values
                            Gamma_base_precision, # Base precision (unscaled)
                            Gamma_base_diag_vec, L_base_factor, # Decomposed base prior parts
                            a, b, # IG prior parameters for sigma2
                            maxiter, tol, sigma2_tol = 1e-4, # Outer loop tolerances
                            inner_maxiter = 10, # Max iterations for inner loops
                            inner_tol_mu = 1e-5,  # Inner tolerance for mu update
                            inner_tol_sigma = 1e-5, # Inner tolerance for Sigma update
                            inner_tol_beta = 1e-5, # Inner tolerance for beta update
                            mu_update_method = "newton", # "newton", "lbfgs", "cg"
                            cg_tol = 1e-5, cg_maxiter = 100, cg_step, 
                            lbfgs_maxit = 100, lbfgs_factr = 1e7, fix_sigma = TRUE
){
  
  if (!mu_update_method %in% c("newton", "lbfgs", "cg","cg1","woodbury")) {
    stop("mu_update_method must be 'newton', 'lbfgs', 'cg', 'cg1' or 'woodbury'")
  }
  if (sigma2_0 <= 0) stop("Initial sigma2_0 must be positive.")
  if (a <= 0 || b <= 0) warning("IG parameters a, b should be positive.")
  if (!requireNamespace("MASS", quietly = TRUE) && mu_update_method == "newton") {
    warning("Package 'MASS' needed for ginv fallback in matrix solves. Please install it.", call. = FALSE)
  }
  
  
  pb = txtProgressBar(style = 3)
  
  # Initialize state
  mu_k = mu0 - mean(mu0)
  Sigma_k_list = Sigma0_list # Includes diag, W, L (should be L_base_factor)
  Sigma = Sigma0_list$diag
  beta_c_k = beta_c0
  sigma2_k = sigma2_0
  
  # # Precompute log determinant of Gamma_base_precision
  # log_det_gamma_base = tryCatch({
  #   Matrix::determinant(Gamma_base_precision, logarithm = TRUE)$modulus[1]
  # }, error = function(e){
  #   warning("Could not compute log determinant of Gamma_base_precision. Setting to 0.", call.=FALSE)
  #   0 # Handle error appropriately
  # })
  
  # Approximate log-determinant of the base precision matrix Gamma.
  # Instead of computing log|Gamma| exactly (which is computationally expensive),
  # we use the sum of log-diagonal entries as a proxy.
  # This ignores the off-diagonal (low-rank) structure and is therefore not exact.
  # However, since this term is constant with respect to all optimization variables
  # (beta, mu, Sigma, sigma2), the approximation does not affect the optimization
  # trajectory or the final solution, and only shifts the ELBO by a constant.
  log_det_gamma_base = sum(log(Gamma_base_diag_vec))
  
  # History tracking (for outer loop changes)
  dmu_norm_outer = c(); dSigma_norm_approx = c(); dbeta_norm_outer = c(); F_ELBO_vec = c()
  sigma2_hist = c(sigma2_k)
  beta_c_list = list(); beta_c_list[[1]] = beta_c0
  
  # Get dimensions and prior rank
  N = dim(X)[1]; n = dim(X_tilde)[1]; N_sum = n + N
  m = ncol(X) # Number of covariates
  r_prior = ncol(L_base_factor)
  I_r = Matrix::Diagonal(r_prior)
  X_star <- as(X_star, "dtCMatrix")
  X_star_t = t(X_star)
  
  # Precompute constants
  tA_eN_vec = c(rep(0, n), rep(1, N)) # Maps to mu_N part
  f_mu_lik_term3 = matrix(tA_eN_vec, ncol=1)
  tX_eN = Matrix::crossprod(X, matrix(1, nrow=N, ncol=1))
  
  # Calculate initial ELBO
  F_ELBO_vec[1] = F_ELBO_lowrank(X, X_tilde, alpha_weight, mu_k, Sigma, beta_c_k, sigma2_k, Gamma_base_diag_vec,
                                 log_det_gamma_base, a, b)
  
  ELBO_tol = tol[1] # Outer loop ELBO tolerance
  beta_tol = tol[2] # Outer loop beta tolerance
  star_time = Sys.time()
  
  # --- Outer Loop ---
  for (i in 1:maxiter) {
    # Store state at the START of the outer iteration for outer convergence check
    mu_outer_old = mu_k
    Sigma_k_diag_outer_old = Sigma
    beta_c_outer_old = beta_c_k
    sigma2_outer_old = sigma2_k
    
    for (inner_iter_beta in 1:inner_maxiter) {
      beta_c_inner_old = beta_c_k
      # Beta update depends on LATEST mu_k, Sigma_k_list, and current sigma2_k
      exp_w_vec_sigma_updated = exp_WeightedVec_lowrank(alpha_weight, mu_k, Sigma, beta_c_k, X_tilde)
      if (any(!is.finite(exp_w_vec_sigma_updated))) {
        warning(paste("Iter", i, "InnerBeta", inner_iter_beta, ": Non-finite exp_w_vec_sigma_updated. Clamping."), call.=FALSE); exp_w_vec_sigma_updated[!is.finite(exp_w_vec_sigma_updated)] <- .Machine$double.xmax }
      
      
      weights_vec = as.numeric(exp_w_vec_sigma_updated)
      if (length(weights_vec) != nrow(X_tilde)) stop("Dim mismatch weights_vec")
      
      # --- Newton Step for Beta ---
      if (m == 1) {
        grad.f_beta_raw = matrix(sum(weights_vec * X_tilde^2), 1, 1)
        # Gradient of positive ELBO w.r.t beta = t(X) 1_N - t(X_tilde) weights
        f_beta = tX_eN - matrix(sum(weights_vec * X_tilde), 1, 1)
      } else {
        Y_weighted = weights_vec * X_tilde
        grad.f_beta_raw = Matrix::crossprod(X_tilde, Y_weighted) # Hessian of NegELBO
        f_beta = tX_eN - Matrix::crossprod(X_tilde, weights_vec) # Gradient of positive ELBO
      }
      
      lambda_reg = 1e-6
      if (is(grad.f_beta_raw, "sparseMatrix")) {
        grad.f_beta_reg = grad.f_beta_raw + lambda_reg * Matrix::Diagonal(m)
      } else {
        grad.f_beta_reg = as.matrix(grad.f_beta_raw); diag(grad.f_beta_reg) = diag(grad.f_beta_reg) + lambda_reg
      }
      
      # Solve Hess(NegELBO) * dbeta = -grad(NegELBO) => grad.f_beta_reg * dbeta = f_beta
      # Newton step dbeta = Hess(NegELBO)^{-1} * f_beta
      dbeta_c = tryCatch({
        # Note: grad.f_beta_reg is Hessian of NegELBO, f_beta is grad(ELBO)
        # We want dbeta = - Hess(ELBO)^{-1} grad(ELBO)
        # Hess(ELBO) = - grad.f_beta_reg (approx, ignoring reg)
        # Solve grad.f_beta_reg * dbeta = f_beta
        Matrix::solve(grad.f_beta_reg, f_beta, sparse = FALSE, tol = .Machine$double.eps)
      }, error = function(e) {
        warning(paste("Iter", i, "InnerBeta", inner_iter_beta, ": Beta solve failed:", e$message), call.=FALSE); matrix(0, nrow = m, ncol = 1)
      })
      #print(norm(f_beta,'F'))
      # --- Line Search Parameters ---
      t_ls <- 1
      alpha_ls <- 1e-4   # Armijo condition parameter (controls sufficient increase)
      beta_ls <- 0.5     # Step size shrink factor
      max_ls_iter <- 20  # Maximum number of backtracking steps
      # Compute current ELBO value
      ELBO_current <- F_ELBO_lowrank(
        X, X_tilde, alpha_weight, mu_k, Sigma,
        beta_c_k, sigma2_k,
        Gamma_base_diag_vec, log_det_gamma_base,
        a, b
      )
      # Directional derivative: gradient^T * direction
      dir_deriv <- as.numeric(crossprod(f_beta, dbeta_c))
      # Safeguard: ensure ascent direction
      if (!is.finite(dir_deriv) || dir_deriv < 0) {
        warning("Non-ascent direction detected in beta update, fallback to small step", call.=FALSE)
        dir_deriv <- abs(dir_deriv)
      }
      # --- Backtracking Line Search Loop ---
      for (ls_iter in 1:max_ls_iter) {
        # Candidate update
        beta_candidate <- beta_c_k + t_ls * dbeta_c
        # Evaluate ELBO at candidate point
        ELBO_candidate <- tryCatch({
          F_ELBO_lowrank(
            X, X_tilde, alpha_weight, mu_k, Sigma,
            beta_candidate, sigma2_k,
            Gamma_base_diag_vec, log_det_gamma_base,
            a, b
          )
        }, error = function(e) {
          -Inf  # Treat failure as unacceptable step
        })
        # Check numerical stability
        if (!is.finite(ELBO_candidate)) {
          t_ls <- t_ls * beta_ls
          next
        }
        if (ELBO_candidate >= ELBO_current + alpha_ls * t_ls * dir_deriv) {
          break
        }
        # Reduce step size
        t_ls <- t_ls * beta_ls
      }
      # Final update
      beta_c_new <- beta_c_k + t_ls * dbeta_c
      beta_step <- t_ls
      beta_c_new = beta_c_k + beta_step*dbeta_c # Newton update: beta_new = beta_old + step
      
      # Check inner convergence for beta
      dbeta_norm_inner = norm(dbeta_c, "F")
      beta_c_k = beta_c_new # Accept inner update
      #print(dbeta_norm_inner)
      if (dbeta_norm_inner < inner_tol_beta) {
        beta_c_list[[length(beta_c_list)+1]] = beta_c_k # Store converged beta for this outer iter
        break # Exit inner beta loop
      }
    } # End inner beta loop
    # If loop finished by maxiter, store last beta
    if(inner_iter_beta == inner_maxiter) {
      beta_c_list[[length(beta_c_list)+1]] = beta_c_k
    }
    
    # --- 1. Inner Loop: Update mu ---
    for (inner_iter_mu in 1:inner_maxiter) {
      mu_inner_old = mu_k # Store mu at start of inner mu iteration
      
      if (mu_update_method == "newton") {
        # --- Newton Step for Mu ---
        exp_w_vec_k = exp_WeightedVec_lowrank(alpha_weight, mu_k, Sigma, beta_c_k, X_tilde)
        if (any(!is.finite(exp_w_vec_k))) {
          warning(paste("Iter", i, "InnerMu", inner_iter_mu, ": Non-finite exp_w_vec_k. Clamping."), call.=FALSE); exp_w_vec_k[!is.finite(exp_w_vec_k)] <- .Machine$double.xmax }
        
        f_mu_lik_term2 = matrix(c(exp_w_vec_k, rep(0, N)), ncol=1)
        #f_mu_lik_term3 = matrix(tA_eN_vec, ncol=1)
        f_mu_prior_term = -(1 / sigma2_k) * (Gamma_base_precision %*% mu_k)
        f_mu = f_mu_lik_term3 - f_mu_lik_term2 + f_mu_prior_term # grad(ELBO)
        
        diag_update_mat = Matrix::Diagonal(n = n, x = as.vector(exp_w_vec_k))
        block_diag_update = Matrix::bdiag(diag_update_mat, Matrix::Matrix(0, nrow=N, ncol=N, sparse=TRUE))
        grad.f_mu_unreg = block_diag_update + (1 / sigma2_k) * Gamma_base_precision # Hessian(ELBO) = H_lik + H_prior
        
        lambda_reg = 1e-6
        N_sum_local = nrow(grad.f_mu_unreg)
        if (is(grad.f_mu_unreg, "sparseMatrix")) {
          grad.f_mu = grad.f_mu_unreg + lambda_reg * Matrix::Diagonal(N_sum_local)
        } else {
          grad.f_mu = grad.f_mu_unreg; diag(grad.f_mu) = diag(grad.f_mu) + lambda_reg }
        
        # dmu = tryCatch({ # H * d = -f_mu => Newton step d = -H^{-1} f_mu
        #   Matrix::solve(grad.f_mu, f_mu, sparse = is(grad.f_mu, "sparseMatrix"), tol = .Machine$double.eps)
        # }, error = function(e) {
        #   warning(paste("Iter", i, "InnerMu", inner_iter_mu, ": Mu Newton solve failed:", e$message), call.=FALSE); matrix(0, nrow=length(mu_k), ncol=1) })
        if (!exists("dmu") || is.null(dmu)) {
          dmu <- Solve_GammaPrecision(f_mu, X_star,X_star_t)
        }
        pcg_result <- pcg_solver(
          grad.f_mu = grad.f_mu,          # Computes A %*% v
          b = f_mu,                              # Right-hand side
          x0 = dmu,# Initial guess
          X_star = X_star,     
          X_star_t = X_star_t,
          M_inv_func = M_inv_func,     # Preconditioner solve M^{-1}v
          tol = cg_tol,
          maxiter = cg_maxiter,verbose = FALSE)
        dmu = pcg_result$x
        
        #print(norm(dmu,'F'))
        # --- Line Search Parameters ---
        t_mu <- 1.0
        alpha_ls_mu <- 1e-4   # Armijo parameter
        beta_ls_mu <- 0.5     # Step shrink factor
        max_ls_iter_mu <- 20
        # Current ELBO
        ELBO_current_mu <- F_ELBO_lowrank(
          X, X_tilde, alpha_weight, mu_k, Sigma,
          beta_c_k, sigma2_k,
          Gamma_base_diag_vec, log_det_gamma_base,
          a, b
        )
        # Directional derivative: grad^T * direction
        dir_deriv_mu <- as.numeric(crossprod(f_mu, dmu))
        # Safeguard: ensure ascent direction
        if (!is.finite(dir_deriv_mu) || dir_deriv_mu < 0) {
          warning(paste("Iter", i, "InnerMu", inner_iter_mu, 
                        ": Non-ascent direction detected, fallback to small step"), call.=FALSE)
          dir_deriv_mu <- abs(dir_deriv_mu)
        }
        # --- Backtracking Line Search Loop ---
        for (ls_iter_mu in 1:max_ls_iter_mu) {
          # Candidate update (NOTE: sign is "+", since dmu = -H^{-1} grad)
          mu_candidate <- mu_k + t_mu * dmu
          # Evaluate ELBO
          ELBO_candidate_mu <- tryCatch({
            F_ELBO_lowrank(
              X, X_tilde, alpha_weight, mu_candidate, Sigma,
              beta_c_k, sigma2_k,
              Gamma_base_diag_vec, log_det_gamma_base,
              a, b
            )
          }, error = function(e) {
            -Inf
          })
          # Numerical stability check
          if (!is.finite(ELBO_candidate_mu)) {
            t_mu <- t_mu * beta_ls_mu
            next
          }
          # Armijo condition
          if (ELBO_candidate_mu >= ELBO_current_mu + alpha_ls_mu * t_mu * dir_deriv_mu) {
            break
          }
          # Reduce step size
          t_mu <- t_mu * beta_ls_mu
        }
        # Final update
        mu_new <- mu_k + t_mu * dmu
        mu_step <- t_mu
        mu_inner_diff = norm(mu_new - mu_inner_old, "F")
        mu_k = mu_new # Accept inner update
        if (mu_inner_diff < inner_tol_mu) {
          break # Exit inner mu loop
        }
        #print(norm(f_mu,'F'))
      } else if (mu_update_method == "lbfgs") {
        # --- L-BFGS Step for Mu ---
        mu_k_vec = as.vector(mu_k)
        optim_result = tryCatch({
          stats::optim(par = mu_k_vec, fn = NegELBO_mu_lowrank, gr = grad_NegELBO_mu_lowrank,
                       method = "L-BFGS-B", control = list(maxit = 1, factr = lbfgs_factr), # NOTE: maxit=1 for single step inside inner loop
                       Sigma_k_list = Sigma_k_list, beta_c_k = beta_c_k, sigma2_k = sigma2_k,
                       Gamma_base_precision = Gamma_base_precision, X_tilde = X_tilde,
                       alpha_weight = alpha_weight, n = n, N = N, tA_eN_vec = tA_eN_vec,
                       L_base_factor = L_base_factor)
        }, error = function(e) {
          warning(paste("Iter", i, "InnerMu", inner_iter_mu, ": optim (L-BFGS) failed:", e$message), call.=FALSE); list(par = mu_k_vec, convergence = -1) })
        if (optim_result$convergence != 0) {
          warning(paste("Iter", i, "InnerMu", inner_iter_mu, ": L-BFGS step warning. Code:", optim_result$convergence), call.=FALSE) }
        mu_new = matrix(optim_result$par, ncol=1)
        
      } else if (mu_update_method == "cg") {
        # --- CG Step for Mu (solving H*d = f_mu for step d) ---
        exp_w_vec_k_cg = exp_WeightedVec_lowrank(alpha_weight, mu_k, Sigma_k_list, beta_c_k, X_tilde)
        f_mu_lik_term2_cg = matrix(c(exp_w_vec_k_cg, rep(0, N)), ncol=1)
        f_mu_lik_term3_cg = matrix(tA_eN_vec, ncol=1)
        f_mu_prior_term_cg = -(1 / sigma2_k) * (Gamma_base_precision %*% mu_k)
        f_mu_cg = f_mu_lik_term3_cg - f_mu_lik_term2_cg + f_mu_prior_term_cg # grad(ELBO)
        
        b_cg = f_mu_cg # Right hand side for H d = b
        dmu0 = matrix(0, nrow = N_sum, ncol = 1) # Initial guess for step d
        # Hessian = build_Hessian_sparse(mu_k = mu_k, Sigma_k_list = Sigma_k_list, beta_c_k = beta_c_k,
        #                                sigma2_k = sigma2_k, Gamma_base_precision = Gamma_base_precision,
        #                                X_tilde = X_tilde, alpha_weight = alpha_weight, n = n, N = N,
        #                                L_base_factor = L_base_factor)
        
        cg_result = cg_solver(Ax_func = HessianVecProd_mu_lowrank, # Ax calculates H*v
                              b = b_cg, x0 = dmu0, tol = cg_tol, maxiter = cg_maxiter, # NOTE: maxiter=1 for single step
                              mu_k = mu_k, Sigma_k_list = Sigma_k_list, beta_c_k = beta_c_k,
                              sigma2_k = sigma2_k, Gamma_base_precision = Gamma_base_precision,
                              X_tilde = X_tilde, alpha_weight = alpha_weight, n = n, N = N,
                              L_base_factor = L_base_factor)
        if (!cg_result$converged && cg_maxiter > 1) { # Only warn if multiple CG iters were intended per step
          warning(paste("Iter", i, "InnerMu", inner_iter_mu, ": CG step did not converge."), call.=FALSE) }
        dmu = cg_result$x # Solution d = H^{-1} f_mu (approx)
        mu_step = cg_step[2]
        mu_new = mu_k + mu_step*dmu # Update mu_new = mu_k + d
        # mu_new = mu_k + dmu
        # # --- Backtracking Line Search ---
        # alpha = 1.0  # Start with full step
        # rho = 0.5    # Backtracking factor
        # c1 = 1e-4    # Armijo condition parameter
        # 
        # # Calculate gradient of Negative ELBO at mu_k
        # g_neg = -grad_NegELBO_mu_lowrank(mu_k, Sigma_k_list, beta_c_k, sigma2_k, Gamma_base_precision,
        #                                  X_tilde, alpha_weight, n, N, tA_eN_vec, L_base_factor)
        # g_neg = as.matrix(g_neg) # Ensure matrix form
        # 
        # # Calculate directional derivative grad(NegELBO)^T * dmu
        # dir_deriv = sum(g_neg * dmu)
        # 
        # # Check if dmu is a descent direction
        # if (dir_deriv >= 0) {
        #   warning(paste("Iteration", i, ": CG direction is not a descent direction (dir_deriv =", dir_deriv,"). Skipping mu update or using gradient step?"))
        #   # Option 1: Skip update
        #   mu_new = mu_k
        #   # Option 2: Use a small gradient step instead (less common with Newton-CG)
        #   # alpha_grad = 1e-4 # Or find via line search
        #   # mu_new = mu_k - alpha_grad * g_neg
        # } else {
        #   # Calculate current NegELBO
        #   neg_elbo_k = NegELBO_mu_lowrank(mu_k, Sigma_k_list, beta_c_k, sigma2_k, Gamma_base_precision,
        #                                   X_tilde, alpha_weight, n, N, tA_eN_vec, L_base_factor)
        #   
        #   max_ls_iter = 10 # Limit line search iterations
        #   ls_iter = 0
        #   while (ls_iter < max_ls_iter) {
        #     mu_trial = mu_k + alpha * dmu
        #     neg_elbo_trial = NegELBO_mu_lowrank(mu_trial, Sigma_k_list, beta_c_k, sigma2_k, Gamma_base_precision,
        #                                         X_tilde, alpha_weight, n, N, tA_eN_vec, L_base_factor)
        #     
        #     # Check Armijo condition
        #     if (neg_elbo_trial <= neg_elbo_k + c1 * alpha * dir_deriv) {
        #       break # Found suitable alpha
        #     }
        #     
        #     alpha = alpha * rho # Backtrack
        #     ls_iter = ls_iter + 1
        #     if (alpha < 1e-9){ # Avoid tiny steps
        #       warning(paste("Iteration", i, ": Line search alpha became too small."))
        #       break
        #     }
        #   }
        #   if(ls_iter == max_ls_iter) warning(paste("Iteration", i, ": Line search reached max iterations."))
        #   
        #   mu_new = mu_k + alpha * dmu # Update mu_new = mu_k + alpha * d
        # }
        # # --- End Line Search ---
      }
      else if (mu_update_method == "woodbury"){
        exp_w_vec_k_cg = exp_WeightedVec_lowrank(alpha_weight, mu_k, Sigma_k_list, beta_c_k, X_tilde)
        f_mu_lik_term2_cg = matrix(c(exp_w_vec_k_cg, rep(0, N)), ncol=1)
        f_mu_lik_term3_cg = matrix(tA_eN_vec, ncol=1)
        #f_mu_prior_term_cg = -(1 / sigma2_k) * (Gamma_base_precision %*% mu_k)
        f_mu_prior_term_cg = -(1 / sigma2_k) * (Gamma_base_diag_vec * mu_k+ L_base_factor %*% (t(L_base_factor) %*% mu_k))
        f_mu_cg = f_mu_lik_term3_cg - f_mu_lik_term2_cg + f_mu_prior_term_cg # grad(ELBO)
        
        b_cg = f_mu_cg # Right hand side for H d = b
        dmu = (Sigma_k_list$diag * as.matrix(b_cg)) - (1/sigma2_k) * Sigma_k_list$diag * (Sigma_k_list$L) %*% (Sigma_k_list$Middle_inv %*% (t(Sigma_k_list$L) %*% (Sigma_k_list$diag * as.matrix(b_cg)))) # Solution d = H^{-1} f_mu (approx)
        mu_step = cg_step[2]
        mu_new = mu_k + mu_step*dmu # Update mu_new
    } else if (mu_update_method == "cg1") {
      # --- Newton Step for Mu ---
      exp_w_vec_k = as.numeric(exp_WeightedVec_lowrank(alpha_weight, mu_k, Sigma, beta_c_k, X_tilde))
      if (any(!is.finite(exp_w_vec_k))) {
        warning(paste("Iter", i, "InnerMu", inner_iter_mu, ": Non-finite exp_w_vec_k. Clamping."), call.=FALSE); exp_w_vec_k[!is.finite(exp_w_vec_k)] <- .Machine$double.xmax }
      
      f_mu_lik_term2 = c(exp_w_vec_k, rep(0, N))
      #f_mu_lik_term3 = matrix(tA_eN_vec, ncol=1)
      f_mu_prior_term = -(1 / sigma2_k) * as.numeric(Gamma_base_precision %*% mu_k)
      f_mu = tA_eN_vec - f_mu_lik_term2 + f_mu_prior_term # grad(ELBO)
      lambda_reg = 1e-6
      
      # diag_update_mat = Matrix::Diagonal(n = n, x = as.vector(exp_w_vec_k))
      # block_diag_update = Matrix::bdiag(diag_update_mat, Matrix::Matrix(0, nrow=N, ncol=N, sparse=TRUE))
      # grad.f_mu_unreg = block_diag_update + (1 / sigma2_k) * Gamma_base_precision # Hessian(ELBO) = H_lik + H_prior
      # 
      
      # N_sum_local = nrow(grad.f_mu_unreg)
      # if (is(grad.f_mu_unreg, "sparseMatrix")) {
      #   grad.f_mu = grad.f_mu_unreg + lambda_reg * Matrix::Diagonal(N_sum_local)
      # } else {
      #   grad.f_mu = grad.f_mu_unreg; diag(grad.f_mu) = diag(grad.f_mu) + lambda_reg }
      H_operator <- function(v) { 
        out <- (1 / sigma2_k) * as.numeric(Gamma_base_precision %*% v) 
        out[1:n] <- out[1:n] + exp_w_vec_k * v[1:n] 
        out <- out + lambda_reg * v 
        return(out) }
      
      # dmu = tryCatch({ # H * d = -f_mu => Newton step d = -H^{-1} f_mu
      #   Matrix::solve(grad.f_mu, f_mu, sparse = is(grad.f_mu, "sparseMatrix"), tol = .Machine$double.eps)
      # }, error = function(e) {
      #   warning(paste("Iter", i, "InnerMu", inner_iter_mu, ": Mu Newton solve failed:", e$message), call.=FALSE); matrix(0, nrow=length(mu_k), ncol=1) })
      if (!exists("dmu_unconstrained") || is.null(dmu_unconstrained)) {
        dmu_unconstrained <- Solve_GammaPrecision(f_mu, X_star,X_star_t)
        #dmu_unconstrained = Matrix::solve(grad.f_mu, f_mu, sparse = is(grad.f_mu, "sparseMatrix"), tol = .Machine$double.eps)
      }
      pcg_unconstrained <- pcg_solver(
        H_operator = H_operator,          # Computes A %*% v
        b = f_mu,                              # Right-hand side
        x0 = dmu_unconstrained,# Initial guess
        X_star = X_star,       
        X_star_t = X_star_t,
        M_inv_func = M_inv_func,     # Preconditioner solve M^{-1}v
        tol = cg_tol,
        maxiter = cg_maxiter,verbose = FALSE)
      dmu_unconstrained = pcg_unconstrained$x
      
      all_one_vec = rep(1,nN)
      if (!exists("dmu_all_one") || is.null(dmu_all_one)) {
        dmu_all_one <- Solve_GammaPrecision(all_one_vec, X_star,X_star_t)
        #dmu_all_one = Matrix::solve(grad.f_mu, all_one_vec, sparse = is(grad.f_mu, "sparseMatrix"), tol = .Machine$double.eps)
      }
      pcg_all_one <- pcg_solver(
        H_operator = H_operator,          # Computes A %*% v
        b = all_one_vec,                              # Right-hand side
        x0 = dmu_all_one,# Initial guess
        X_star = X_star,       
        X_star_t = X_star_t,
        M_inv_func = M_inv_func,     # Preconditioner solve M^{-1}v
        tol = cg_tol,
        maxiter = cg_maxiter,verbose = FALSE)
      dmu_all_one = pcg_all_one$x
      schur_component = sum(dmu_all_one)
      if (schur_component < 1e-8){
        dmu = dmu_unconstrained - mean(dmu_unconstrained)
      } else {
        dmu = dmu_unconstrained - sum(dmu_unconstrained)/schur_component * dmu_all_one
      }
      
      
      #print(norm(dmu,'F'))
      # --- Line Search Parameters ---
      t_mu <- 1.0
      alpha_ls_mu <- 1e-4   # Armijo parameter
      beta_ls_mu <- 0.5     # Step shrink factor
      max_ls_iter_mu <- 20
      # Current ELBO
      ELBO_current_mu <- F_ELBO_lowrank(
        X, X_tilde, alpha_weight, mu_k, Sigma,
        beta_c_k, sigma2_k,
        Gamma_base_diag_vec, log_det_gamma_base,
        a, b
      )
      # Directional derivative: grad^T * direction
      dir_deriv_mu <- as.numeric(crossprod(f_mu, dmu))
      # Safeguard: ensure ascent direction
      if (!is.finite(dir_deriv_mu) || dir_deriv_mu < 0) {
        warning(paste("Iter", i, "InnerMu", inner_iter_mu, 
                      ": Non-ascent direction detected, fallback to small step"), call.=FALSE)
        dir_deriv_mu <- abs(dir_deriv_mu)
      }
      # --- Backtracking Line Search Loop ---
      for (ls_iter_mu in 1:max_ls_iter_mu) {
        # Candidate update (NOTE: sign is "+", since dmu = -H^{-1} grad)
        mu_candidate <- mu_k + t_mu * dmu
        # Evaluate ELBO
        ELBO_candidate_mu <- tryCatch({
          F_ELBO_lowrank(
            X, X_tilde, alpha_weight, mu_candidate, Sigma,
            beta_c_k, sigma2_k,
            Gamma_base_diag_vec, log_det_gamma_base,
            a, b
          )
        }, error = function(e) {
          -Inf
        })
        # Numerical stability check
        if (!is.finite(ELBO_candidate_mu)) {
          t_mu <- t_mu * beta_ls_mu
          next
        }
        # Armijo condition
        if (ELBO_candidate_mu >= ELBO_current_mu + alpha_ls_mu * t_mu * dir_deriv_mu) {
          break
        }
        # Reduce step size
        t_mu <- t_mu * beta_ls_mu
      }
      # Final update
      mu_new <- mu_k + t_mu * dmu
      mu_step <- t_mu
      mu_inner_diff = norm(mu_new - mu_inner_old, "F")
      mu_k = mu_new # Accept inner update
      if (mu_inner_diff < inner_tol_mu) {
        break # Exit inner mu loop
      }
      #print(norm(f_mu,'F'))
    }
    
    #Check inner convergence for mu
    mu_inner_diff = norm(mu_new - mu_inner_old, "F")
    mu_k = mu_new # Accept inner update
    
    if (mu_inner_diff < inner_tol_mu) {
      break # Exit inner mu loop
    }
    # End inner mu loop
    }
    
    # --- 2. Update Sigma (using updated mu_k, current sigma2_k) ---
    # Uses Woodbury update based on Sigma^{-1} = H = H_lik + H_prior
    # H = block_diag(exp_w_vec, 0) + (1/sigma2_k) * Gamma_base
    # H = block_diag(...) + (1/sigma2_k) * (Diag(G_diag) + L L^T)
    # H = Diag( H_diag' ) + (1/sigma2_k) * L L^T
    # where H_diag' = block_diag(...)_diag + (1/sigma2_k)*G_diag
    for (inner_iter_Sigma in 1:inner_maxiter){
      Sigma_k_diag_old = Sigma
      exp_w_vec_mu_updated = exp_WeightedVec_lowrank(alpha_weight, mu_k, Sigma, beta_c_k, X_tilde) # Pass args
      
      diag_lik_update_vec = c(as.vector(exp_w_vec_mu_updated), rep(0, N))
      H_tilde_diag_vec = diag_lik_update_vec + (1/sigma2_k) * Gamma_base_diag_vec
      H_tilde_diag_vec = pmax(H_tilde_diag_vec, .Machine$double.eps) # Ensure positivity
      
      D_inv_vec = 1 / H_tilde_diag_vec
      Sigma = D_inv_vec
      dSigma_diag_norm = norm(Sigma - Sigma_k_diag_old, "2")
      #print(dSigma_diag_norm)
      if(dSigma_diag_norm < 1e-6){
        break
      }
    }
    
    
    dSigma_diag_norm = norm(Sigma - Sigma_k_diag_outer_old, "2")
    dSigma_W_norm = dSigma_diag_norm
    dSigma_norm_approx_val = sqrt(dSigma_diag_norm^2 + dSigma_W_norm^2)
    dSigma_norm_approx = append(dSigma_norm_approx, dSigma_norm_approx_val, length(dSigma_norm_approx))
    
    
    
    if(fix_sigma == FALSE){
    # --- 4. Update sigma^2 (Once per outer iteration, using converged mu & Sigma) ---
    muT_G_mu_k = as.numeric(Matrix::crossprod(mu_k, Gamma_base_precision %*% mu_k))
    tr_G_Sigma_k = trace_Gamma_Sigma_lowrank(Gamma_base_diag_vec, Sigma)
    
    if (is.na(tr_G_Sigma_k)) {
      warning(paste("Iter", i, ": trace_Gamma_Sigma failed during sigma2 update. Skipping sigma2 update."), call.=FALSE)
      sigma2_new = sigma2_k # Keep old value
    } else {
      numerator_s2 = n + N + 2*a + 2
      denominator_s2 = muT_G_mu_k + tr_G_Sigma_k - 2*b
      
      if(denominator_s2 <= 0){
        warning(paste("Iter", i, ": Denominator for sigma2 update non-positive (", denominator_s2, "). Clamping sigma2."), call.=FALSE)
        # Clamp or keep old:
        # sigma2_new = sigma2_k
        sigma2_new = max(sigma2_k * 0.1, .Machine$double.eps) # Cautious decrease
      } else {
        sigma2_new = denominator_s2 / numerator_s2
      }
      sigma2_new = max(.Machine$double.eps, sigma2_new) # Ensure positivity
    }
      sigma2_k = sigma2_new # Accept update
    }
    
    sigma2_hist = append(sigma2_hist, sigma2_k, length(sigma2_hist))
    
    
    # --- Calculate Outer Loop Changes and ELBO ---
    dmu_norm_outer = append(dmu_norm_outer, norm(mu_k - mu_outer_old, "F"), length(dmu_norm_outer))
    
    dSigma_diag_norm_outer = norm(Sigma - Sigma_k_diag_outer_old, "2")
    dSigma_W_norm_outer = dSigma_diag_norm_outer
    
    
    dbeta_norm_outer = append(dbeta_norm_outer, norm(beta_c_k - beta_c_outer_old, "F"), length(dbeta_norm_outer)) # Overall change in beta
    
    F_ELBO_current = F_ELBO_lowrank(X, X_tilde, alpha_weight, mu_k, Sigma, beta_c_k, sigma2_k, Gamma_base_diag_vec,
                                    log_det_gamma_base, a, b)
    F_ELBO_vec = append(F_ELBO_vec, F_ELBO_current, length(F_ELBO_vec))
    
    # --- Outer Convergence Check ---
    # Use ELBO change from previous *outer* iteration's ELBO
    elbo_diff = abs(F_ELBO_current - F_ELBO_vec[length(F_ELBO_vec) - 1])
    # Use overall change in parameters from start to end of outer iteration
    beta_diff_outer = norm(beta_c_k - beta_c_outer_old, "F") # Or max(abs(...))
    sigma2_diff_outer = abs(sigma2_k - sigma2_outer_old)
    
    if(!is.na(elbo_diff) && elbo_diff <= ELBO_tol &&
       beta_diff_outer <= beta_tol &&
       sigma2_diff_outer <= sigma2_tol &&
       dSigma_diag_norm_outer <= 1e-5) {
      message("\nOuter loop convergence criteria met.")
      break
    }
    
    setTxtProgressBar(pb, i/maxiter)
  } # --- End of Outer loop ---
  
  end_time = Sys.time()
  close(pb)
  run_time = as.numeric(difftime(end_time, star_time, units = "secs"))
  
  # Return results including sigma2
  result = list(mu_k, dmu_norm_outer, Sigma_k_list, dSigma_norm_approx, beta_c_list, dbeta_norm_outer,
                sigma2_k, sigma2_hist, # Add sigma2 results
                F_ELBO_vec, run_time, mu_update_method)
  names(result) = c("mu", "dmu_norm", "Sigma_lowrank", "dSigma_norm_approx", "beta_hist", "dbeta_norm",
                    "sigma2", "sigma2_hist", # Add sigma2 names
                    "ELBO", "running_time", "mu_update_method")
  return(result)
}


# Example Usage (requires setup of inputs like decompose_gamma)

#' Decompose a sparse precision matrix into Diagonal + Low Rank approximation
#' (Copied from thought process - adjust inputs/outputs for base precision)
decompose_gamma <- function(Gamma_base_precision, rank_r, tol_eigen = sqrt(.Machine$double.eps)) {
  # ... (implementation as in thought process, ensuring output names are
  #      Gamma_base_diag_vec and L_base_factor) ...
  if (!is(Gamma_base_precision, "symmetricMatrix")) {
    if (!Matrix::isSymmetric(Gamma_base_precision)) stop("Input Gamma_base_precision must be symmetric.")
    Gamma_base_precision <- Matrix::forceSymmetric(Gamma_base_precision)
  }
  N_sum <- nrow(Gamma_base_precision)
  if(N_sum == 0) stop("Gamma_base_precision is empty.")
  
  message("Step 1: Extracting diagonal...")
  Gamma_base_diag_vec <- 1 * Matrix::diag(Gamma_base_precision)
  if (any(Gamma_base_diag_vec <= tol_eigen)) {
    warning("Some diagonal elements <= tol_eigen. Clamping."); Gamma_base_diag_vec[Gamma_base_diag_vec <= tol_eigen] <- tol_eigen }
  Gamma_diag_mat <- Matrix::Diagonal(n = N_sum, x = Gamma_base_diag_vec)
  
  message("Step 2: Calculating off-diagonal remainder matrix R...")
  R_offdiag <- Gamma_base_precision - Gamma_diag_mat
  if (!is(R_offdiag, "symmetricMatrix")) R_offdiag <- Matrix::forceSymmetric(R_offdiag)
  R_offdiag_for_eigs <- as(R_offdiag, "dgCMatrix") # Ensure dgCMatrix for RSpectra
  
  message(paste("Step 3: Computing partial eigen decomp (k =", rank_r, ")..."))
  tryCatch({
    eigen_decomp <- RSpectra::eigs_sym(R_offdiag_for_eigs, k = rank_r, which = "LM", opts = list(retvec = TRUE, ncv = min(N_sum, max(2 * rank_r + 1, 20))))
  }, error = function(e){ stop(paste("RSpectra::eigs_sym failed:", e$message)) })
  
  message("Step 4: Constructing L_base_factor...")
  eigenvalues <- eigen_decomp$values; eigenvectors <- eigen_decomp$vectors
  positive_idx <- which(eigenvalues > tol_eigen)
  
  if (length(positive_idx) == 0) {
    warning("No positive eigenvalues found. L_base_factor will be empty."); L_base_factor <- matrix(0.0, nrow = N_sum, ncol = 0); r_actual <- 0
  } else {
    pos_eigenvalues <- eigenvalues[positive_idx]; pos_eigenvectors <- eigenvectors[, positive_idx, drop = FALSE]
    order_idx <- order(pos_eigenvalues, decreasing = TRUE)
    pos_eigenvalues <- pos_eigenvalues[order_idx]; pos_eigenvectors <- pos_eigenvectors[, order_idx, drop = FALSE]
    L_base_factor <- sweep(pos_eigenvectors, MARGIN = 2, STATS = sqrt(pos_eigenvalues), FUN = `*`)
    r_actual <- ncol(L_base_factor); message(paste("  Actual rank r_actual =", r_actual))
  }
  message("Decomposition complete.")
  return(list(
    Gamma_base_diag_vec = Gamma_base_diag_vec,
    L_base_factor = L_base_factor,
    r_actual = r_actual
  ))
}
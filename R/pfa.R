#' Fit PFA (Perturbed Factor Analysis; script backend)
#'
#' PFA in the tutorial is provided as scripts (R + C++ via Rcpp). This wrapper expects
#' the PFA scripts to be bundled under `inst/extdata/pfa/` in this package source.
#' At runtime, the backend is loaded via [load_backend()], which sources all `.R` files
#' in that directory and compiles `PFA.cpp` (if present) via `Rcpp::sourceCpp()`.
#'
#' @param Y_list List of study matrices.
#' @param ... Passed to `PFA()` in the backend script. User can set hyper-parameters and numbers of iterations here.
#' @param k Positive integer number of shared factors.
#' @param center Logical; whether to center each variable within each study.
#' @param scale Logical; whether to scale each variable within each study.
#' @param nrun Number of MCMC iterations.
#' @param burn Number of burn-in iterations.
#' @return A `bifa_fit` object.
#' @export
fit_pfa <- function(
    Y_list = NULL,
    k = NULL,
    center = TRUE,
    scale = FALSE,
    nrun = 3000,
    burn = 2000,
    ...
) {
  if (is.null(Y_list)) stop("`Y_list` is required.", call. = FALSE)
  Y_list <- .as_matrix_list(Y_list)
  
  # Require k
  if (is.null(k) || length(k) != 1 || !is.finite(k) || k <= 0) {
    stop("PFA requires a positive integer number of factors: supply `k`.", call. = FALSE)
  }
  k <- as.integer(k)
  
  # Load PFA backend from inst/extdata/pfa/
  env <- load_backend("pfa")
  if (!exists("PFA", envir = env, inherits = FALSE)) {
    stop(
      "Could not find function `PFA` in the PFA backend script.\n",
      "Please ensure the PFA scripts are present under inst/extdata/pfa/ in the package source\n",
      "(e.g., FBPFA-PFA.R plus any helper .R files, and optionally PFA.cpp).",
      call. = FALSE
    )
  }
  
  # Scale each study, then stack
  Y_list_scaled <- lapply(Y_list, function(x) {
    scale(as.matrix(x), center = center, scale = scale)
  })
  X <-  as.matrix(do.call(rbind, Y_list_scaled))
  
  # Study id vector
  id <- rep(seq_along(Y_list_scaled), vapply(Y_list_scaled, nrow, integer(1)))
  
  pfa_fun <- get("PFA", envir = env, inherits = FALSE)
  fit <- pfa_fun(Y=t(X),
                 latentdim = k,
                 grpind = id,
                 Thin = 5,
                 Cutoff = 0.001,
                 Total_itr = nrun, burn = burn)
  
  .new_bifa_fit(
    "pfa",
    fit,
    meta = list(
      S = length(Y_list),
      P = ncol(Y_list[[1]]),
      k = k,
      center = center,
      scale = scale
    )
  )
}


#' Post-process PFA output (nutrition tutorial)
#'
#' This closely follows the code shown in the nutrition case study.
#'
#' @param fit A `bifa_fit` from `fit_pfa()` or a raw PFA fit.
#' @return A list with common/study-specific loadings and covariances.
#' @export
postprocess_pfa <- function(fit) {
  .require_pkg("MSFA", "PFA post-processing uses MSFA::sp_OP")
  fit <- .unwrap_fit(fit)
  if (is.null(fit$Loading) || is.null(fit$Latentsigma) || is.null(fit$Errorsigma) || is.null(fit$Pertmat)) {
    stop("PFA fit must contain $Loading, $Latentsigma, $Errorsigma, $Pertmat.", call. = FALSE)
  }
  # Determine posterior dimension (number of factors per sample)
  k_vec <- sapply(fit$Loading, ncol)
  mode_k <- as.numeric(names(sort(table(k_vec), decreasing = TRUE)[1]))
  
  # Filter posterior samples to those with mode_k
  keep_idx <- which(k_vec == mode_k)
  fit$Loading <- fit$Loading[keep_idx]
  fit$Latentsigma <- fit$Latentsigma[keep_idx]
  fit$Errorsigma <- fit$Errorsigma[keep_idx]
  fit$Pertmat <- fit$Pertmat[keep_idx]
  
  npost <- length(fit$Loading)
  p <- nrow(fit$Loading[[1]])
  k <- mode_k
  S <- dim(fit$Pertmat[[1]])[2]
  
  posteriorPhis <- array(0, dim = c(p, k, npost))
  posteriorLams <- vector("list", S)
  
  for(s in 1:S){
    posteriorLams[[s]] <- array(0, dim = c(p, k, npost))
    for(i in 1:npost){
      posteriorPhis[,,i] <- fit$Loading[[i]] %*% diag(fit$Latentsigma[[i]])
      posteriorLams[[s]][,,i] <- (solve(matrix(fit$Pertmat[[i]][, s], p, p)) - diag(p)) %*% posteriorPhis[,,i]
    }
  }
  
  # Varimax rotation
  est_Phi <- MSFA::sp_OP(posteriorPhis, itermax = 10, trace = FALSE)$Phi
  # For study-specific loadings: 
  # Study 1 should be the same as the shared loading
  # Only apply rotation to studies 2:S
  est_speLoad <- vector("list", S)
  est_speLoad[[1]] <- est_Phi
  
  if(S > 1) {
    for(s in 2:S) {
      est_speLoad[[s]] <- MSFA::sp_OP(posteriorLams[[s]], itermax = 10, trace = FALSE)$Phi
    }
  }
  # Estimated covariance components
  sharevar <- list()
  est_SigmaLambdaList <- vector("list", S)
  est_SigmaMarginal <- vector("list", S)
  est_Psi_list <- list()
  
  for(s in 1:S){
    post_SigmaLambda_s <- vector("list", npost)
    post_SigmaMarginal_s <- vector("list", npost)
    Psi <- vector("list", npost)
    
    for(i in 1:npost){
      sharevar[[i]] <- fit$Loading[[i]] %*% diag(fit$Latentsigma[[i]]^2) %*% t(fit$Loading[[i]]) + 
        diag(fit$Errorsigma[[i]]^2)
      Q_temp_inv <- solve(matrix(fit$Pertmat[[i]][, s], p, p))
      if (s==1) {
        post_SigmaLambda_s[[i]] = sharevar[[i]]
        post_SigmaMarginal_s[[i]] <- sharevar[[i]]
      } else {
        post_SigmaMarginal_s[[i]] <- Q_temp_inv %*% sharevar[[i]] %*% t(Q_temp_inv)
        post_SigmaLambda_s[[i]] <- post_SigmaMarginal_s[[i]] - sharevar[[i]]
      }
      Psi[[i]] <- diag(fit$Errorsigma[[i]]^2)
    }
    
    est_SigmaMarginal[[s]] <- Reduce('+', post_SigmaMarginal_s) / npost
    est_SigmaLambdaList[[s]] <- Reduce('+', post_SigmaLambda_s) / npost
    est_Psi_list[[s]] <- Reduce('+', Psi) / npost
  }
  
  est_Psi <- Reduce('+', est_Psi_list) / S
  est_SigmaPhi <- Reduce('+', sharevar) / npost
  est_Q <- Reduce('+', fit$Pertmat) / npost
  est_Q_list <- lapply(1:S, function(s) matrix(est_Q[, s], p, p))
  
  list(
    Phi = est_Phi,
    SigmaPhi = est_SigmaPhi,
    Psi = est_Psi,
    Q = est_Q_list,
    LambdaList = est_speLoad,
    SigmaLambdaList = est_SigmaLambdaList,
    SigmaMarginal = est_SigmaMarginal,
    mode_k = mode_k,
    kept_samples = length(keep_idx)
  )
}

#' Select K for PFA via nontrivial-column rule
#'
#' In the nutrition case study, a factor is considered "kept" if at least one
#' loading exceeds a small threshold in absolute value.
#'
#' @param post Output of `postprocess_pfa()`.
#' @param threshold Threshold for considering a column nontrivial (default 1e-3).
#' @return Integer K.
#' @export
select_k_pfa <- function(post, threshold = 1e-3) {
  Phi <- post$Phi
  keep_col <- apply(abs(Phi) > threshold, 2, any)
  sum(keep_col)
}

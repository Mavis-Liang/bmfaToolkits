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
  K_set <- vapply(fit$Loading, ncol, integer(1))
  mode_k <- as.numeric(names(which.max(table(K_set))))
  kept <- which(K_set == mode_k)

  Loading <- fit$Loading[kept]
  Latentsigma <- fit$Latentsigma[kept]
  Errorsigma <- fit$Errorsigma[kept]
  Pertmat <- fit$Pertmat[kept]
  npost <- length(kept)
  p <- nrow(Loading[[1]])
  S <- ncol(Pertmat[[1]]) / p

  posteriorPhis <- array(NA, dim = c(p, mode_k, npost))
  posteriorLams <- lapply(seq_len(S), function(s) array(NA, dim = c(p, mode_k, npost)))

  sharevar <- vector("list", npost)
  post_SigmaLambda <- lapply(seq_len(S), function(s) vector("list", npost))
  post_SigmaMarginal <- lapply(seq_len(S), function(s) vector("list", npost))
  Psi <- vector("list", npost)

  for (i in seq_len(npost)) {
    posteriorPhis[,,i] <- Loading[[i]] %*% diag(Latentsigma[[i]])
    sharevar[[i]] <- Loading[[i]] %*% diag(Latentsigma[[i]]^2) %*% t(Loading[[i]]) +
      diag(Errorsigma[[i]]^2)

    for (s in seq_len(S)) {
      Q_temp_inv <- solve(matrix(Pertmat[[i]][, s], p, p))
      post_SigmaMarginal[[s]][[i]] <- Q_temp_inv %*% sharevar[[i]] %*% t(Q_temp_inv)
      post_SigmaLambda[[s]][[i]] <- post_SigmaMarginal[[s]][[i]] - sharevar[[i]]
    }
    Psi[[i]] <- diag(Errorsigma[[i]]^2)
  }

  est_Phi <- MSFA::sp_OP(posteriorPhis, itermax = 10, trace = FALSE)$Phi
  est_speLoad <- lapply(posteriorLams, function(x) MSFA::sp_OP(x, itermax = 10, trace = FALSE)$Phi)

  est_SigmaPhi <- Reduce("+", sharevar) / npost
  est_SigmaLambdaList <- lapply(seq_len(S), function(s) Reduce("+", post_SigmaLambda[[s]]) / npost)
  est_SigmaMarginal <- lapply(seq_len(S), function(s) Reduce("+", post_SigmaMarginal[[s]]) / npost)

  est_Psi_list <- lapply(seq_len(S), function(s) diag(diag(est_SigmaMarginal[[s]] - est_SigmaPhi)))
  est_Psi <- Reduce("+", est_Psi_list) / S

  est_Q <- Reduce("+", Pertmat) / npost
  est_Q_list <- lapply(seq_len(S), function(s) matrix(est_Q[, s], p, p))

  list(Phi = est_Phi,
       SigmaPhi = est_SigmaPhi,
       LambdaList = est_speLoad,
       SigmaLambdaList = est_SigmaLambdaList,
       PsiList = est_Psi_list,
       Psi = est_Psi,
       QList = est_Q_list,
       SigmaMarginal = est_SigmaMarginal,
       mode_k = mode_k)
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

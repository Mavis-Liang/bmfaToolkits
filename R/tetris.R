#' Fit Tetris (single run)
#'
#' Runs one Tetris fit. Optionally uses a fixed A matrix (selected by choose.A()).
#'
#' @param Y_list List of study matrices.
#' @param alpha Either numeric alpha, or "auto" to use ceiling(1.25 * S).
#' @param beta Beta hyperparameter for Tetris.
#' @param fixed_bigT Logical; whether to run with fixed big T.
#' @param bigT Fixed T matrix used when fixed=TRUE.
#' @param nrun Number of iterations.
#' @param burn Burn-in iterations.
#' @param nprint Print frequency.
#' @param ... Passed to backend tetris().
#' @return A `bifa_fit` object.
#' @export
fit_tetris <- function(
    Y_list,
    alpha = c("auto"),
    beta = 1,
    fixed_bigT = FALSE,
    bigT = NULL,
    nrun = 3000,
    burn = 2000,
    nprint = 200,
    ...
) {
  Y_list <- .as_matrix_list(Y_list)
  env <- load_backend("tetris")
  
  if (!exists("tetris", envir = env, inherits = FALSE)) {
    stop("Could not find function `tetris` in the Tetris backend script. ",
         "Ensure Tetris scripts are under inst/extdata/tetris/.", call. = FALSE)
  }
  
  S <- length(Y_list)
  
  # alpha handling
  if (is.character(alpha)) {
    if (length(alpha) != 1 || alpha != "auto") {
      stop("`alpha` must be numeric or 'auto'.", call. = FALSE)
    }
    alpha_val <- ceiling(1.25 * S)
  } else {
    if (length(alpha) != 1 || !is.finite(alpha) || alpha <= 0) {
      stop("`alpha` must be a positive number (or 'auto').", call. = FALSE)
    }
    alpha_val <- alpha
  }
  
  # fixed handling
  if (isTRUE(fixed_bigT) && is.null(bigT)) {
    stop("If fixed_bigT=TRUE, you must provide `bigT`.", call. = FALSE)
  }
  
  dots <- list(...)
  # set defaults unless user overrides in ...
  if (is.null(dots$alpha))  dots$alpha  <- alpha_val
  if (is.null(dots$beta))   dots$beta   <- beta
  if (is.null(dots$nrun))   dots$nrun   <- nrun
  if (is.null(dots$burn))   dots$burn   <- burn
  if (is.null(dots$nprint)) dots$nprint <- nprint
  if (is.null(dots$fixed_bigT))  dots$fixed  <- fixed_bigT
  if (isTRUE(fixed_bigT) && is.null(dots$A_fixed)) dots$A_fixed <- bigT
  
  fit <- do.call(get("tetris", envir = env), c(list(Y_list), dots))
  
  .new_bifa_fit(
    "tetris",
    fit,
    meta = list(S = S, P = ncol(Y_list[[1]]), alpha = alpha_val, beta = beta,
                nrun = nrun, burn = burn, nprint = nprint)
  )
}

#' Select big_T (A matrix) for Tetris
#'
#' Uses choose.A() on a Tetris fit to select the A matrix.
#'
#' @param fit A `bifa_fit` returned by [fit_tetris()] (typically the non-fixed run).
#' @param alpha_IBP Alpha used by choose.A(). If NULL, uses fit metadata alpha.
#' @param S Number of studies. If NULL, uses fit metadata S.
#' @param ... Passed to choose.A().
#' @return The selected A matrix (big_T).
#' @export
select_T_tetris <- function(fit, alpha_IBP = NULL, S = NULL, ...) {
  if (!inherits(fit, "bifa_fit") || fit$backend != "tetris") {
    stop("`fit` must be a 'bifa_fit' with model == 'tetris'.", call. = FALSE)
  }
  
  env <- load_backend("tetris")
  if (!exists("choose.A", envir = env, inherits = FALSE)) {
    stop("Could not find function `choose.A` in the Tetris backend script.", call. = FALSE)
  }
  
  if (is.null(alpha_IBP)) alpha_IBP <- fit$meta$alpha
  if (is.null(S)) S <- fit$meta$S
  
  if (is.null(alpha_IBP) || !is.finite(alpha_IBP)) stop("Need `alpha_IBP`.", call. = FALSE)
  if (is.null(S) || !is.finite(S)) stop("Need `S`.", call. = FALSE)
  
  do.call(get("choose.A", envir = env), c(list(fit$fit), list(alpha_IBP = alpha_IBP, S = S), list(...)))
}


#' Post-process Tetris output
#'
#' The tutorial extracts study-specific loadings via `getLambda(fit)` from the Tetris
#' script and builds marginal covariances as Lambda Lambda^T + Psi.
#'
#' @param fit A `bifa_fit` from `fit_tetris()` or a raw Tetris fit.
#' @param num_samps Number of posterior samples to use. If NULL, uses min(2000, available).
#' @return A list with `Phi`, `SigmaPhi`, `LambdaList`, `SigmaLambdaList`, `PsiList`, `SigmaMarginal`.
#' @export
postprocess_tetris <- function(fit, num_samps = NULL) {
  fit <- .unwrap_fit(fit)
  
  # choose a fixed A matrix
  A <- fit$A
  if (is.list(A)) A <- A[[1]]
  if (is.null(A) || !is.matrix(A)) stop("Tetris fit has no usable A.", call. = FALSE)
  
  env <- load_backend("tetris")
  if (!exists("getLambda", envir = env, inherits = FALSE)) {
    stop("`getLambda` not found in Tetris backend.", call. = FALSE)
  }
  getLambda_fun <- get("getLambda", envir = env, inherits = FALSE)
  
  # ---- key: pick num_samps within available range
  n_avail <- if (!is.null(fit$Lambda)) length(fit$Lambda) else 0L
  if (n_avail < 1) stop("Tetris fit has no Lambda draws (fit$Lambda is empty).", call. = FALSE)
  
  if (is.null(num_samps)) num_samps <- min(2000L, n_avail)
  num_samps <- min(as.integer(num_samps), n_avail)
  
  # ---- do NOT swallow the real error; surface it
  Lambda <- tryCatch(
    getLambda_fun(fit, A, num_samps = num_samps),
    error = function(e) stop("getLambda() failed: ", e$message, call. = FALSE)
  )
  Lambda <- as.matrix(Lambda)
  
  S <- nrow(A)
  
  # shared
  est_Phi <- Lambda[, colSums(A) == S, drop = FALSE]
  est_SigmaPhi <- tcrossprod(est_Phi)
  
  # study-specific
  P_shared <- diag(as.integer(colSums(A) == S))
  T_s <- vector("list", S)
  est_LambdaList <- vector("list", S)
  for (s in 1:S) {
    T_s[[s]] <- diag(A[s, ])
    Lambda_s <- Lambda %*% (T_s[[s]] - P_shared)
    keep <- colSums(abs(Lambda_s) > 0) > 0
    est_LambdaList[[s]] <- Lambda_s[, keep, drop = FALSE]
  }
  est_SigmaLambdaList <- lapply(est_LambdaList, tcrossprod)
  
  # marginal
  Psi <- vector("list", S)
  est_SigmaMarginal <- vector("list", S)
  for (s in 1:S) {
    Psi[[s]] <- diag(Reduce("+", fit$Psi[[s]]) / length(fit$Psi[[s]]))
    est_SigmaMarginal[[s]] <- Lambda %*% T_s[[s]] %*% t(Lambda) + Psi[[s]]
  }
  
  list(
    Phi = est_Phi, SigmaPhi = est_SigmaPhi,
    LambdaList = est_LambdaList, SigmaLambdaList = est_SigmaLambdaList,
    Psi = Psi, T_s = T_s,
    SigmaMarginal = est_SigmaMarginal
  )
}



#' Run the full Tetris pipeline: initial fit -> choose T -> fixed refit
#'
#' Mirrors the common workflow:
#' 1) tetris(..., fixed_bigT=FALSE)
#' 2) big_T <- choose.A(fit, ...)
#' 3) tetris(..., fixed_bigT=TRUE, A_fixed=big_T)
#'
#' @param Y_list List of study matrices.
#' @param alpha Either numeric alpha, or "auto" to use ceiling(1.25 * S).
#' @param beta Beta hyperparameter.
#' @param nrun,burn,nprint MCMC controls for both runs.
#' @param chooseA_args List of extra arguments passed to choose.A().
#' @param ... Extra args passed to tetris() in both runs.
#' @return A `bifa_fit` (fixed refit) with `stage1` and `big_T` attached.
#' @export
fit_tetris_2step <- function(
    Y_list,
    alpha = "auto",
    beta = 1,
    nrun = 3000,
    burn = 2000,
    nprint = 200,
    chooseA_args = list(),
    ...
) {
  # stage 1: initial fit (not fixed)
  fit1 <- fit_tetris(
    Y_list = Y_list,
    alpha = alpha,
    beta = beta,
    fixed_bigT = FALSE,
    nrun = nrun,
    burn = burn,
    nprint = nprint,
    ...
  )
  
  # choose T
  big_T <- do.call(
    select_T_tetris,
    c(list(fit = fit1), chooseA_args)
  )
  
  # stage 2: fixed refit
  fit2 <- fit_tetris(
    Y_list = Y_list,
    alpha = fit1$meta$alpha,   # ensure numeric alpha carried through
    beta = beta,
    fixed_bigT = TRUE,
    bigT = big_T,
    nrun = nrun,
    burn = burn,
    nprint = nprint,
    ...
  )
  postprocess_tetris(fit2)
}


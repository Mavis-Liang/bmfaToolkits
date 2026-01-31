#' Fit VI-MSFA CAVI (variational MSFA)
#'
#' Thin wrapper around `VIMSFA::cavi_msfa()` (from blhansen/VI-MSFA).
#'
#' Conventions in this package:
#' - Input is `Y_list`: a list of S matrices (n_s x P) with common P (same variables and column order).
#' - We optionally center/scale each study matrix *before* calling CAVI.
#' - The returned object is wrapped as a `bifa_fit` with backend `"cavi"`.
#'
#' @param Y_list List of study matrices (each `n_s x P`).
#' @param k Integer; number of common factors (`K` in VI-MSFA).
#' @param j_s Integer; study-specific factor counts (`J_s` in VI-MSFA). Either length 1 or length S.
#' @param scaling Logical; whether to scale each study matrix (default FALSE).
#' @param centering Logical; whether to center each study matrix (default TRUE).
#' @param center,scale Deprecated aliases for `centering` / `scaling` (kept for backward compatibility).
#' @param ... Passed to `VIMSFA::cavi_msfa()`.
#' @return A `bifa_fit` object.
#' @export
fit_cavi <- function(Y_list, k, j_s,
                     scaling = FALSE, centering = TRUE,
                     center = NULL, scale = NULL,
                     ...) {
  .require_pkg("VIMSFA", "CAVI backend (VI-MSFA)")

  # Back-compat: if the caller used the old argument names, honor them
  if (missing(centering) && !is.null(center)) centering <- isTRUE(center)
  if (missing(scaling)   && !is.null(scale))  scaling   <- isTRUE(scale)

  Y_list <- .as_matrix_list(Y_list)
  S <- length(Y_list)
  if (length(j_s) == 1) j_s <- rep(as.integer(j_s), S)
  if (length(j_s) != S) stop("j_s must be length 1 or length S.", call. = FALSE)

  Y_list_scaled <- lapply(Y_list, function(x) scale(as.matrix(x), center = centering, scale = scaling))

  dots <- list(...)
  # VI-MSFA has an internal `scale` option; default FALSE here because we've preprocessed.
  if (is.null(dots$scale)) dots$scale <- FALSE

  fit <- tryCatch(
    do.call(VIMSFA::cavi_msfa, c(list(X_s = Y_list_scaled, K = as.integer(k), J_s = as.integer(j_s)), dots)),
    error = function(e) {
      stop(
        "Calling VIMSFA::cavi_msfa failed.\n",
        "Tip: install the backend with bmfaToolkits::install_backend('cavi').\n\n",
        "Original error: ", e$message,
        call. = FALSE
      )
    }
  )

  .new_bifa_fit("cavi", fit, meta = list(S = S, P = ncol(Y_list[[1]]), k = as.integer(k), j_s = as.integer(j_s)))
}

#' Post-process CAVI output
#'
#' Extracts posterior means from `VIMSFA::cavi_msfa()` output and builds:
#' - `Phi` and `SigmaPhi = Phi Phi^T`
#' - `LambdaList` and `SigmaLambdaList[[s]] = Lambda_s Lambda_s^T`
#' - `PsiList` (study-specific residual variances)
#' - `SigmaMarginal[[s]] = SigmaPhi + SigmaLambda_s + diag(Psi_s)`
#'
#' This standardizes CAVI output so that downstream code (metrics/plots) can treat methods uniformly.
#'
#' @param fit A `bifa_fit` from `fit_cavi()` or a raw `VIMSFA::cavi_msfa()` fit object.
#' @return A list with `Phi`, `SigmaPhi`, `LambdaList`, `SigmaLambdaList`, `PsiList`, `SigmaMarginal`.
#' @export
post_CAVI <- function(fit) {
  fit <- .unwrap_fit(fit)
  stopifnot(is.list(fit))

  if (is.null(fit$mean_phi) || is.null(fit$mean_lambda_s) || is.null(fit$mean_psi_s)) {
    stop("CAVI fit object missing mean_phi / mean_lambda_s / mean_psi_s.", call. = FALSE)
  }

  Phi <- as.matrix(fit$mean_phi)                 # P x K
  SigmaPhi <- tcrossprod(Phi)

  LambdaList <- fit$mean_lambda_s                # list of (P x J_s)
  S <- length(LambdaList)

  SigmaLambdaList <- lapply(LambdaList, function(L) {
    if (is.null(L)) return(NULL)
    tcrossprod(as.matrix(L))
  })

  psi_list <- fit$mean_psi_s                     # list of length-P vectors
  P <- nrow(Phi)
  PsiList <- lapply(psi_list, function(v) as.numeric(v))

  SigmaMarginal <- lapply(seq_len(S), function(s) {
    SigmaPhi + SigmaLambdaList[[s]] + diag(PsiList[[s]], nrow = P)
  })

  list(
    Phi = Phi,
    SigmaPhi = SigmaPhi,
    LambdaList = LambdaList,
    SigmaLambdaList = SigmaLambdaList,
    PsiList = PsiList,
    SigmaMarginal = SigmaMarginal
  )
}

#' Post-process CAVI output (standard name)
#'
#' @param fit A `bifa_fit` from `fit_cavi()` or a raw CAVI fit.
#' @param ... Ignored (kept for interface consistency).
#' @return Same as `post_CAVI()`.
#' @export
postprocess_cavi <- function(fit, ...) {
  post_CAVI(fit)
}

#' Select K for CAVI from the common covariance
#'
#' Uses the eigen-proportion rule on `SigmaPhi`.
#'
#' @param post Output of `postprocess_cavi()` / `post_CAVI()`.
#' @param cutoff Eigen proportion cutoff.
#' @export
select_k_cavi <- function(post, cutoff = 0.05) {
  select_k_from_sigma(post$SigmaPhi, cutoff = cutoff)
}

#' Select J_s for CAVI from study-specific covariances
#'
#' Uses the eigen-proportion rule on each `SigmaLambdaList[[s]]`.
#'
#' @param post Output of `postprocess_cavi()` / `post_CAVI()`.
#' @param cutoff Eigen proportion cutoff.
#' @export
select_js_cavi <- function(post, cutoff = 0.05) {
  lapply(post$SigmaLambdaList, select_k_from_sigma, cutoff = cutoff)
}

#' Refit CAVI after selecting the number of factors
#'
#' This mirrors the BMSFA \"2-step\" workflow:
#' 1) Fit CAVI with an initial (potentially over-specified) `(k, j_s)`.\n
#' 2) Post-process and select `K` and `J_s` from the estimated covariance components.\n
#' 3) Refit CAVI with the selected factor counts.\n
#'
#' @param Y_list Data list.
#' @param fit0 Optional initial fit (from `fit_cavi()`).
#' @param post_fit0 Optional post-processed initial fit (from `postprocess_cavi()`). If both are supplied,
#'   `post_fit0` is used.
#' @param cutoff Eigen proportion cutoff.
#' @param ... Passed to `fit_cavi()` when refitting.
#' @return Post-processed output of the refit.
#' @export
refit_cavi <- function(Y_list, fit0 = NULL, post_fit0 = NULL, cutoff = 0.05, ...) {
  if (is.null(post_fit0)) {
    if (is.null(fit0)) stop("Provide `fit0` or `post_fit0` for refit_cavi().", call. = FALSE)
    post_fit0 <- postprocess_cavi(fit0)
  }

  K  <- select_k_cavi(post_fit0, cutoff = cutoff)
  js <- unlist(select_js_cavi(post_fit0, cutoff = cutoff))

  message(sprintf("Selected K = %d", K))
  message(sprintf("Selected j_s = c(%s)", paste(js, collapse = ", ")))

  fit2 <- fit_cavi(Y_list, k = K, j_s = js, ...)
  postprocess_cavi(fit2)
}

#' Fit CAVI with automatic factor selection (2-step)
#'
#' Convenience wrapper: fit an initial CAVI, select `K` and `J_s`, then refit.\n
#' If you already have an initial post-processed fit (e.g. from a fast run), pass it via `post_fit0`
#' to skip the initial fitting stage.
#'
#' @param Y_list Data list.
#' @param post_fit0 Optional post-processed initial fit (or NULL to compute internally).
#' @param k Initial number of common factors (required if `post_fit0` is NULL).
#' @param j_s Initial study-specific factors (required if `post_fit0` is NULL).
#' @param cutoff Eigen proportion cutoff.
#' @param ... Passed to `fit_cavi()` for both stages.
#' @return Post-processed output of the refit.
#' @export
fit_cavi_2step <- function(Y_list, post_fit0 = NULL, k = NULL, j_s = NULL, cutoff = 0.05, ...) {
  if (is.null(post_fit0)) {
    if (is.null(k) || is.null(j_s)) {
      stop("If `post_fit0` is NULL, you must provide initial `k` and `j_s`.", call. = FALSE)
    }
    fit0 <- fit_cavi(Y_list, k = k, j_s = j_s, ...)
    post_fit0 <- postprocess_cavi(fit0)
    return(refit_cavi(Y_list, post_fit0 = post_fit0, cutoff = cutoff, ...))
  }
  refit_cavi(Y_list, post_fit0 = post_fit0, cutoff = cutoff, ...)
}

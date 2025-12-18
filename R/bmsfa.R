#' Fit BMSFA via MSFA::sp_msfa
#'
#' @param Y_list List of study matrices.
#' @param k Number of common factors.
#' @param j_s Study-specific factors (length S).
#' @param outputlevel,scaling,centering,control Passed to `MSFA::sp_msfa()`.
#' @param ... Passed to `MSFA::sp_msfa()`.
#' @return A `bifa_fit` object.
#' @export
fit_bmsfa <- function(Y_list, k, j_s,
                      outputlevel = 1,
                      scaling = FALSE, centering = TRUE,
                      control = list(nrun = 3000, burn = 2000),
                      ...) {
  .require_pkg("MSFA", "BMSFA backend")
  Y_list <- .as_matrix_list(Y_list)
  S <- length(Y_list)
  if (length(j_s) == 1) j_s <- rep(as.integer(j_s), S)
  if (length(j_s) != S) stop("j_s must be length 1 or length S.", call. = FALSE)

  fit <- MSFA::sp_msfa(Y_list, k = k, j_s = j_s,
                       outputlevel = outputlevel,
                       scaling = scaling, centering = centering,
                       control = control, ...)
  .new_bifa_fit("bmsfa", fit, meta = list(S = S, P = ncol(Y_list[[1]]), k = k, j_s = j_s))
}

#' Post-process BMSFA output
#'
#' OP on Phi and each Lambda_s, form
#' SigmaPhi, SigmaLambdaList, estimate Psi_s by averaging `fit$psi[[s]]` over iterations,
#' and build marginal covariance as SigmaPhi + SigmaLambda_s + diag(Psi_s).
#'
#' @param fit A `bifa_fit` from `fit_bmsfa()` or a raw `MSFA::sp_msfa()` fit.
#' @return A list with `Phi`, `SigmaPhi`, `LambdaList`, `SigmaLambdaList`, `PsiList`, `SigmaMarginal`.
#' @export
postprocess_bmsfa <- function(fit) {
  .require_pkg("MSFA", "BMSFA post-processing")
  fit <- .unwrap_fit(fit)

  est_Phi <- MSFA::sp_OP(fit$Phi, trace = FALSE)$Phi
  est_SigmaPhi <- tcrossprod(est_Phi)

  est_LambdaList <- lapply(fit$Lambda, function(x) MSFA::sp_OP(x, trace = FALSE)$Phi)
  est_SigmaLambdaList <- lapply(est_LambdaList, function(x) tcrossprod(x))

  S <- length(est_SigmaLambdaList)
  est_PsiList <- lapply(seq_len(S), function(s) {
    apply(fit$psi[[s]], c(1, 2), mean)
  })

  est_margin_cov <- lapply(seq_len(S), function(s) {
    est_SigmaPhi + est_SigmaLambdaList[[s]] + diag(as.vector(est_PsiList[[s]]))
  })

  list(Phi = est_Phi,
       SigmaPhi = est_SigmaPhi,
       LambdaList = est_LambdaList,
       SigmaLambdaList = est_SigmaLambdaList,
       PsiList = est_PsiList,
       SigmaMarginal = est_margin_cov)
}

#' Select K for BMSFA from common covariance
#' @param post Output of `postprocess_bmsfa()`.
#' @param cutoff Eigen proportion cutoff.
#' @export
select_k_bmsfa <- function(post, cutoff = 0.05) {
  select_k_from_sigma(post$SigmaPhi, cutoff = cutoff)
}

#' Select J_s for BMSFA from study-specific covariances
#' @param post Output of `postprocess_bmsfa()`.
#' @param cutoff Eigen proportion cutoff.
#' @export
select_js_bmsfa <- function(post, cutoff = 0.05) {
  lapply(post$SigmaLambdaList, select_k_from_sigma, cutoff = cutoff)
}

#' Fit BMSFA all-at-once
#' If `post_fit0` is NULL, fit an initial BMSFA with given `k` and `j_s`, then automatically select
#' the number of factors and refit.
#' Otherwise, use the `post_fit0` then select the number of factors and refit.
#'
#' @param Y_list Data list.
#' @param post_fit0 Post-processed initial fit (or NULL).
#' @param k Initial number of common factors (if `post_fit0` is NULL).
#' @param j_s Initial study-specific factors (if `post_fit0` is NULL).
#' @param cutoff Eigen proportion cutoff.
#' @param ... Passed to `fit_bmsfa()`.
#' @export
fit_bmsfa_2step <- function(Y_list, post_fit0 = NULL, k = NULL, j_s = NULL, cutoff = 0.05, ...) {
  if (is.null(post_fit0)) {
    # quick initial fit; user should pass a better starting point if needed
    fit0 <- fit_bmsfa(Y_list, k = k, j_s = j_s, ...)
    post_fit0 <- postprocess_bmsfa(fit0)
  }
  K <- select_k_bmsfa(post_fit0, cutoff = cutoff)
  js <- unlist(select_js_bmsfa(post_fit0, cutoff = cutoff))
  message(sprintf("Selected K = %d", K))
  message(sprintf("Selected j_s = c(%s)", paste(js, collapse = ", ")))
  fit2 <- fit_bmsfa(Y_list, k = K, j_s = js, ...)
  postprocess_bmsfa(fit2)
}

#' Fit independent factor analysis (Ind FA) via MSFA::sp_fa (per study)
#'
#' @param Y_list List of study matrices.
#' @param j_s Either a single integer or a length-S integer vector of factors.
#' @param scaling,centering,control Passed to `MSFA::sp_fa()`.
#' @param ... Passed to `MSFA::sp_fa()`.
#' @return A `bifa_fit` object whose `$fit` is a list of fits (one per study).
#' @export
fit_ind_fa <- function(Y_list, j_s, scaling = FALSE, centering = TRUE,
                       control = list(nrun = 3000, burn = 2000), ...) {
  .require_pkg("MSFA", "Ind FA backend")
  Y_list <- .as_matrix_list(Y_list)
  S <- length(Y_list)
  if (length(j_s) == 1) j_s <- rep(as.integer(j_s), S)
  if (length(j_s) != S) stop("j_s must be length 1 or length S.", call. = FALSE)

  fit_list <- lapply(seq_len(S), function(s) {
    MSFA::sp_fa(Y_list[[s]], k = j_s[s], scaling = scaling, centering = centering,
                control = control, ...)
  })
  .new_bifa_fit("ind_fa", fit_list, meta = list(S = S, P = ncol(Y_list[[1]]), j_s = j_s))
}

#' Post-process Ind FA output (nutrition tutorial)
#'
#' @param fit A `bifa_fit` from `fit_ind_fa()` or a raw list of `MSFA::sp_fa()` fits.
#' @return A list with `LambdaList`, `SigmaLambdaList`, `Psi` (list), `SigmaMarginal` (list).
#' @export
postprocess_ind_fa <- function(fit) {
  .require_pkg("MSFA", "Ind FA post-processing")
  fit_list <- .unwrap_fit(fit)
  S <- length(fit_list)

  est_LambdaList <- lapply(seq_len(S), function(s) {
    MSFA::sp_OP(fit_list[[s]]$Lambda, trace = FALSE)$Phi
  })
  est_SigmaLambdaList <- lapply(est_LambdaList, function(x) tcrossprod(x))

  est_SigmaMarginal <- lapply(seq_len(S), function(s) {
    apply(fit_list[[s]]$Sigma, c(1, 2), mean)
  })

  Psi <- vector("list", S)
  for (s in seq_len(S)) {
    Psi_chain <- vector("list", dim(fit_list[[s]]$Sigma)[3])
    for (i in seq_len(dim(fit_list[[s]]$Sigma)[3])) {
      Psi_chain[[i]] <- fit_list[[s]]$Sigma[, , i] - tcrossprod(fit_list[[s]]$Lambda[, , i])
    }
    Psi[[s]] <- Reduce("+", Psi_chain) / length(Psi_chain)
  }

  list(LambdaList = est_LambdaList,
       SigmaLambdaList = est_SigmaLambdaList,
       Psi = Psi,
       SigmaMarginal = est_SigmaMarginal)
}

#' Select J_s for Ind FA (per-study) using eigen proportion rule
#'
#' @param post Output of `postprocess_ind_fa()`.
#' @param cutoff Proportion threshold.
#' @export
select_js_ind_fa <- function(post, cutoff = 0.05) {
  lapply(post$SigmaLambdaList, select_k_from_sigma, cutoff = cutoff)
}

#' Fit Ind FA all at once
#' If `post_fit0` is NULL, fit an initial IndFA with given `j_s`, then automatically select
#' the number of factors and refit.
#'
#' @param Y_list Data list.
#' @param post_fit0 Initial post-processing result (or NULL to compute internally).
#' @param j_s Initial number of factors (or NULL to use all).
#' @param cutoff Eigen proportion cutoff.
#' @param ... Passed to `fit_ind_fa()`.
#' @export
fit_ind_fa_2step <- function(Y_list, post_fit0 = NULL, j_s=NULL, cutoff = 0.05, ...) {
  if (is.null(post_fit0)) {
    fit0 <- fit_ind_fa(Y_list, j_s = j_s, ...)
    post_fit0 <- postprocess_ind_fa(fit0)
  }
  js <- unlist(select_js_ind_fa(post_fit0, cutoff = cutoff))
  message(sprintf("Selected j_s = c(%s)\n", paste(js, collapse = ", ")))
  fit2 <- fit_ind_fa(Y_list, j_s = js, ...)
  postprocess_ind_fa(fit2)
}

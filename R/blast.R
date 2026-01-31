#' Fit BLAST (script backend)
#'
#' BLAST is a research-code repository (not a CRAN/Bioconductor package). In this package we
#' treat BLAST as a **script backend**, similar to PFA / Tetris:
#'
#' - When developing from source: put the BLAST scripts under `inst/extdata/blast/`.
#' - After installation: they will be available at `system.file("extdata","blast", package="bmfaToolkits")`.
#' - At runtime, `fit_blast()` loads the backend via [load_backend()] and calls the backend's `fit_blast()`.
#'
#' @param Y_list List of study matrices (each n_s x P).
#' @param k Number of common factors.
#' @param q_s Study-specific factor counts (length 1 or length S).
#' @param n_MC Number of BLAST Monte Carlo iterations.
#' @param center Whether to center each study matrix (default TRUE).
#' @param scale Whether to scale each study matrix (default FALSE).
#' @param subsample_index Integer vector of variable indices used by BLAST when `sample_outer_product = FALSE`.
#'   If `NULL`, we use `1:min(100, P)` to match the BLAST defaults used in this project.
#' @param sample_outer_product Logical, forwarded to BLAST.
#' @param ... Passed to BLAST's `fit_blast()` implementation.
#' @return A `bifa_fit` object.
#' @export
fit_blast <- function(
    Y_list,
    k,
    q_s,
    n_MC = 10000,
    center = TRUE,
    scale = FALSE,
    subsample_index = NULL,
    sample_outer_product = FALSE,
    ...
) {
  Y_list <- .as_matrix_list(Y_list)
  S <- length(Y_list)
  if (length(q_s) == 1) q_s <- rep(as.integer(q_s), S)
  if (length(q_s) != S) stop("q_s must be length 1 or length S.", call. = FALSE)

  Y_list_scaled <- lapply(Y_list, function(x) scale(as.matrix(x), center = center, scale = scale))
  P <- ncol(Y_list_scaled[[1]])
  if (is.null(subsample_index)) subsample_index <- seq_len(min(100L, P))

  env <- load_backend("blast")
  backend_fun <- get("fit_blast", envir = env, inherits = TRUE)

  fit <- tryCatch(
    backend_fun(
      Y = Y_list_scaled,
      k = as.integer(k),
      q_s = as.integer(q_s),
      n_MC = as.integer(n_MC),
      subsample_index = subsample_index,
      sample_outer_product = isTRUE(sample_outer_product),
      ...
    ),
    error = function(e) {
      stop(
        "Calling BLAST's fit_blast() failed.\n",
        "Tip: ensure the BLAST scripts are present under inst/extdata/blast/ (source) or in the installed package extdata/blast/ directory.\n\n",
        "Original error: ", e$message,
        call. = FALSE
      )
    }
  )

  .new_bifa_fit("blast", fit, meta = list(S = S, P = P, k = as.integer(k), q_s = as.integer(q_s)))
}

#' Post-process BLAST output
#'
#' BLAST objects are not standardized, so this helper extracts posterior means and returns
#' a common structure used throughout this package:
#'
#' - `Phi` (common loadings) and `SigmaPhi`
#' - `LambdaList` (study-specific loadings) and `SigmaLambdaList`
#' - `PsiList` (residual variances; if unavailable, `NA` is returned)
#' - `SigmaMarginal` per study
#'
#' Behind the scenes, BLAST often provides covariance components directly as posterior means
#' of outer products (e.g., `Lambda_outer_mean` and `Gammas_outer_mean`). If these are present,
#' we use them; otherwise we compute `tcrossprod()` from the corresponding loadings means.
#'
#' @param fit A `bifa_fit` from `fit_blast()` or a raw BLAST fit list.
#' @return A list with `Phi`, `SigmaPhi`, `LambdaList`, `SigmaLambdaList`, `PsiList`, `SigmaMarginal`.
#' @export
post_BLAST <- function(fit) {
  fit <- .unwrap_fit(fit)
  stopifnot(is.list(fit))

  # --- infer S from fitted output ---
  S <- NA_integer_
  if (!is.null(fit$Ms) && is.list(fit$Ms)) {
    S <- length(fit$Ms)
  } else if (!is.null(fit$Fs) && is.list(fit$Fs)) {
    S <- length(fit$Fs)
  } else if (!is.null(fit$Gammas_outer_mean)) {
    S <- dim(fit$Gammas_outer_mean)[3]
  } else if (!is.null(fit$Gammas_mean)) {
    S <- dim(fit$Gammas_mean)[3]
  }
  if (is.na(S) || length(S) != 1L) stop("post_BLAST(): Cannot infer scalar number of studies S.", call. = FALSE)

  # --- common loadings and common covariance ---
  if (is.null(fit$Lambda_mean)) stop("post_BLAST(): fit$Lambda_mean not found", call. = FALSE)
  Phi <- as.matrix(fit$Lambda_mean)

  # BLAST may provide posterior mean of Lambda Lambda^T directly
  if (is.null(fit$Lambda_outer_mean)) {
    SigmaPhi <- tcrossprod(Phi)
  } else {
    SigmaPhi <- as.matrix(fit$Lambda_outer_mean)
  }

  P <- nrow(Phi)

  # --- infer q_s (study-specific factor counts) ---
  q_s <- rep(0L, S)
  if (!is.null(fit$Fs) && is.list(fit$Fs) && length(fit$Fs) == S) {
    q_s <- as.integer(sapply(fit$Fs, function(Fs) if (is.null(Fs)) 0L else ncol(as.matrix(Fs))))
  } else if (!is.null(fit$Gammas_mean)) {
    q_s <- rep(as.integer(dim(fit$Gammas_mean)[2]), S)
  }

  # --- study-specific loadings and covariance ---
  LambdaList <- vector("list", S)
  SigmaLambdaList <- vector("list", S)

  if (!is.null(fit$Gammas_mean)) {
    Gm <- fit$Gammas_mean  # P x max(q_s) x S
    for (s in seq_len(S)) {
      js <- q_s[s]
      if (js <= 0L) {
        LambdaList[[s]] <- matrix(0, nrow = P, ncol = 0)
      } else {
        LambdaList[[s]] <- matrix(Gm[, seq_len(js), s], nrow = P, ncol = js)
      }
    }
  } else {
    for (s in seq_len(S)) LambdaList[[s]] <- matrix(0, nrow = P, ncol = 0)
  }

  if (!is.null(fit$Gammas_outer_mean)) {
    for (s in seq_len(S)) SigmaLambdaList[[s]] <- as.matrix(fit$Gammas_outer_mean[, , s, drop = FALSE][,,1])
  } else {
    SigmaLambdaList <- lapply(LambdaList, function(L) tcrossprod(L))
  }

  # --- residual variances ---
  psi_hat <- rep(NA_real_, P)
  if (!is.null(fit$Sigma_2s_samples)) {
    Sig <- as.matrix(fit$Sigma_2s_samples)  # N_mc x P
    if (ncol(Sig) == P) psi_hat <- colMeans(Sig)
  }
  PsiList <- replicate(S, psi_hat, simplify = FALSE)

  # --- marginal covariance per study ---
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

#' Post-process BLAST output (standard name)
#'
#' @param fit A `bifa_fit` from `fit_blast()` or a raw BLAST fit.
#' @param ... Ignored (kept for interface consistency).
#' @return Same as `post_BLAST()`.
#' @export
postprocess_blast <- function(fit, ...) {
  post_BLAST(fit)
}

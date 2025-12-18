#' Fit SUFA
#'
#' Thin wrapper around `SUFA::fit_SUFA()`
#'
#' @param Y_list List of study matrices.
#' @param k Positive integer for maximum number of factors shared (`qmax` in SUFA).
#' @param nrun Number of MCMC iterations.
#' @param center Whether to center each study matrix (default TRUE).
#' @param scale Whether to scale each study matrix (default FALSE).
#' @param ... Passed to `SUFA::fit_SUFA()`.
#' @return A `bifa_fit` object.
#' @export
fit_sufa <- function(
    Y_list,
    k = NULL,
    nrun = 3000,
    center = TRUE,
    scale = FALSE,
    ...
) {
  .require_pkg("SUFA", "SUFA backend")
  if (is.null(k) || length(k) != 1 || !is.finite(k) || k <= 0) {
    stop("SUFA requires a positive integer `k` (or alias `k`).", call. = FALSE)
  }

  qmax <- as.integer(k)
  
  Y_list <- .as_matrix_list(Y_list)
  
  # Scale each study: default center=TRUE, scale=FALSE
  Y_list_scaled <- lapply(Y_list, function(x) {
    scale(as.matrix(x), center = center, scale = scale)
  })
  
  dots <- list(...)
  
  # If user didn't provide, set defaults that match your workflow
  if (is.null(dots$qmax)) dots$qmax <- k
  if (is.null(dots$nrun)) dots$nrun <- nrun
  
  # Call SUFA::fit_SUFA with first arg as the scaled Y_list
  fit <- tryCatch(
    do.call(SUFA::fit_SUFA, c(list(Y_list_scaled), dots)),
    error = function(e) {
      stop(
        "Calling SUFA::fit_SUFA failed.\n",
        "Tip: try specifying `qmax` and `nrun` (and check SUFA installation/system deps).\n\n",
        "Original error: ", e$message,
        call. = FALSE
      )
    }
  )
  
  .new_bifa_fit(
    "sufa",
    fit,
    meta = list(
      S = length(Y_list_scaled),
      P = ncol(Y_list_scaled[[1]]),
      k = k,
      nrun = nrun,
      center = center,
      scale = scale
    )
  )
}


#' Post-process SUFA output
#'
#' @param fit A `bifa_fit` from `fit_sufa()` or a raw SUFA fit.
#' @return A list with `Phi`, `SigmaPhi`, `SigmaMarginal`, `PsiList`.
#' @export
postprocess_sufa <- function(fit) {
  .require_pkg("SUFA", "SUFA post-processing")
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("SUFA post-processing requires 'doParallel'. Install it with install.packages('doParallel').",
         call. = FALSE)
  }
  # SUFA calls library(doParallel) internally; make sure it's attached
  if (!"package:doParallel" %in% search()) {
    suppressPackageStartupMessages(library(doParallel))
  }
  
  fit <- .unwrap_fit(fit)

  all <- dim(fit$Lambda)[3]
  burnin <- floor(all * 0.8) # We will use the last 20% samples
  # shared and study-specific loading matrices
  loadings <- SUFA::lam.est.all(fit, burn = burnin)
  # Obtain common covariance matrix and loading from fitting
  est_Phi <- loadings$Shared
  est_SigmaPhi <- SUFA::SUFA_shared_covmat(fit, burn = burnin)
  est_Psi <- diag(colMeans(fit$residuals))
  # Study-specific loadings
  est_LambdaList <- loadings$Study_specific
  
  # Obtain study-specific covariance matrices
  S <- length(fit$A)
  marginal_cov <- SUFA::sufa_marginal_covs(fit, burn = burnin)
  est_SigmaLambdaList <- list()
  for (s in 1:S) {
    est_SigmaLambdaList[[s]] <- marginal_cov[,,s] - est_SigmaPhi
  }
  
  return(list(SigmaPhi = est_SigmaPhi, Phi = est_Phi, 
              SigmaLambdaList = est_SigmaLambdaList,
              LambdaList = est_LambdaList, 
              Psi = est_Psi,
              SigmaMarginal = lapply(1:S, function(s) marginal_cov[,,s])
  ))
}

#' Fit MOM-SS (MOM-SS / BFR.BE backend)
#'
#' This is a thin wrapper around `BFR.BE::BFR.BE.EM.CV()` (or similar entry points).
#' Because the upstream API may evolve, we keep this wrapper flexible:
#' - If you pass `...` with `X` and `M` (the membership matrix), they are forwarded.
#' - If you pass `Y_list`, we will stack it and construct `X` and `M` accordingly.
#'
#' @param Y_list Optional list of study matrices. If provided, we create a stacked matrix.
#' @param k Positive integer number of shared factors.
#' @param X Optional covariate matrix or list of covariate matrices aligned to `Y_list`.
#' @param scaling Logical; whether to scale the data.
#' @param center_within_study Logical; whether to center each study's data before fitting. Recommended unless you have strong reasons not to.
#' @param folds Number of folds for cross-validation within MOM-SS.
#' @param ... Parameters that can be pass directly to `BFR.BE::BFR.BE.EM.CV()`.
#' @return A `bifa_fit` object.
#' @export
fit_momss <- function(
    Y_list,
    k,
    X = NULL,
    scaling = FALSE,
    center_within_study = TRUE,
    folds = 10,
    ...
) {
  .require_pkg("BFR.BE", "MOM-SS backend (BFR.BE)")
  
  Y_list <- .as_matrix_list(Y_list)
  if (center_within_study) {
    Y_list <- lapply(Y_list, function(x) scale(x, center = TRUE, scale = FALSE))
  }
  
  x <- as.matrix(do.call(rbind, Y_list))
  storage.mode(x) <- "double"
  
  N_s <- vapply(Y_list, nrow, integer(1))
  b <- as.matrix(Matrix::bdiag(lapply(N_s, function(n) matrix(1, nrow = n, ncol = 1))))
  storage.mode(b) <- "double"
  
  v <- NULL
  if (!is.null(X)) {
    v <- if (is.list(X)) as.matrix(do.call(rbind, lapply(X, as.matrix))) else as.matrix(X)
    storage.mode(v) <- "double"
    if (nrow(v) != nrow(x)) stop("`X` must match nrow(rbind(Y_list)).", call. = FALSE)
  }
  
  fit <- .with_safe_mvtnorm_logdet(function() {
    .with_attached_pkgs(c("plyr","MASS","mombf","matrixcalc","Matrix","BFR.BE"), function() {
      BFR.BE::BFR.BE.EM.CV(
        x = x, v = v, b = b, q = k,
        scaling = scaling, folds = folds, ...
      )
    })
  })
  
  .new_bifa_fit("momss", fit, meta = list(S = length(Y_list), P = ncol(Y_list[[1]]), k = k))
}



#' Post-process MOM-SS output
#'
#' The tutorial uses `fit$M` as the loading matrix and `fit$Sigma` as the list of
#' study-specific marginal covariances.
#'
#' @param fit A `bifa_fit` from `fit_momss()` or a raw MOM-SS fit.
#' @return A list with `Phi`, `SigmaPhi`, `SigmaMarginal`, `PsiList`.
#' @export
postprocess_momss <- function(fit) {
  fit <- .unwrap_fit(fit)

  # Loadings
  est_Phi <- fit$Mpost
  if (is.null(est_Phi)) stop("Cannot find loadings matrix in MOM-SS fit (expected $Mpost).", call. = FALSE)

  est_SigmaPhi <- tcrossprod(est_Phi)
  
  # Marginal covariance
  S <- dim(fit$sigma)[2]
  est_PsiList <- est_SigmaMarginal <-  list()
  for(s in 1:S){
    est_PsiList[[s]] <- fit$sigma[,s]
    est_SigmaMarginal[[s]] <- est_SigmaPhi + diag(fit$sigma[,s])
  }
  # last S columns of fit$Theta are the study-specific intercepts
  est_alphas <- fit$Theta[, (dim(fit$Theta)[2]-S+1):dim(fit$Theta)[2]]
  # The rest are coeficients for the known covariates
  est_B <- fit$Theta[, 1:(dim(fit$Theta)[2]-S)]
  
  return(list(Phi = est_Phi, SigmaPhi = est_SigmaPhi, Psi = est_PsiList, alpha = est_alphas, B = est_B,
              SigmaMarginal = est_SigmaMarginal))
}

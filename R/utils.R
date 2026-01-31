#' Internal helpers (not exported)
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

.require_pkg <- function(pkg, why = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0(
      "Required package '", pkg, "' is not installed.\n",
      if (!is.null(why)) paste0("Needed for: ", why, "\n") else "",
      "Try: bmfaToolkits::install_backend(\"", tolower(pkg), "\") or install it manually."
    )
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}

# For MOMSS that needs several packages to be attached. (Because it does not use 
# namespaces properly)
.with_attached_pkgs <- function(pkgs, thunk) {
  stopifnot(is.function(thunk))
  pkgs <- unique(pkgs)
  newly <- pkgs[!paste0("package:", pkgs) %in% search()]
  
  for (p in newly) {
    suppressPackageStartupMessages(
      library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
    )
  }
  
  on.exit({
    for (p in rev(newly)) {
      tag <- paste0("package:", p)
      if (tag %in% search()) {
        # If another package requires it (e.g., glmnet -> Matrix), detach() will error.
        # We MUST NOT fail during cleanup.
        try(detach(tag, character.only = TRUE, unload = FALSE), silent = TRUE)
      }
    }
  }, add = TRUE)
  
  thunk()
}


# MOMSS's logdet function from mvtnorm fails on ordinary matrices. We temporarily 
# replace it with a version that works for both.
.with_safe_mvtnorm_logdet <- function(thunk) {
  stopifnot(is.function(thunk))
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Need package 'mvtnorm' installed.", call. = FALSE)
  }
  
  ns <- asNamespace("mvtnorm")
  
  # Save old binding
  old_logdet <- get("logdet", envir = ns, inherits = FALSE)
  
  # Replacement that works for ordinary matrices too
  new_logdet <- function(x, ...) {
    if (is.matrix(x) || inherits(x, "Matrix")) {
      x <- as.matrix(x)
      storage.mode(x) <- "double"
      
      # fast path for diagonal matrices
      if (nrow(x) == ncol(x) &&
          all(x[upper.tri(x)] == 0) &&
          all(x[lower.tri(x)] == 0)) {
        return(sum(log(diag(x))))
      }
      
      return(as.numeric(determinant(x, logarithm = TRUE)$modulus))
    }
    old_logdet(x, ...)
  }
  
  # Temporarily swap it in mvtnorm namespace (binding exists, so this is allowed)
  unlockBinding("logdet", ns)
  assign("logdet", new_logdet, envir = ns)
  lockBinding("logdet", ns)
  
  on.exit({
    unlockBinding("logdet", ns)
    assign("logdet", old_logdet, envir = ns)
    lockBinding("logdet", ns)
  }, add = TRUE)
  
  thunk()
}

# For MSFA: if it's DeVito's version or Mavis's version that supports centering/scaling
.msfa_capabilities <- function() {
  f <- getFromNamespace("sp_msfa", "MSFA")
  args <- names(formals(f))
  list(
    has_center = "center" %in% args,
    has_scale  = "scale"  %in% args,
    version    = as.character(utils::packageVersion("MSFA"))
  )
}


.find_pkg_root <- function(start = getwd()) {
  p <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in 1:15) {
    if (file.exists(file.path(p, "DESCRIPTION"))) return(p)
    parent <- dirname(p)
    if (identical(parent, p)) break
    p <- parent
  }
  ""
}


.as_matrix_list <- function(Y) {
  if (is.matrix(Y)) return(list(Y))
  if (is.data.frame(Y)) return(list(as.matrix(Y)))
  if (is.list(Y)) {
    out <- lapply(Y, function(x) {
      if (is.data.frame(x)) x <- as.matrix(x)
      if (!is.matrix(x)) stop("Each element of Y_list must be a matrix/data.frame.", call. = FALSE)
      x
    })
    return(out)
  }
  stop("Input must be a matrix, data.frame, or list of matrices/data.frames.", call. = FALSE)
}

.check_same_p <- function(Y_list) {
  p <- ncol(Y_list[[1]])
  if (!all(vapply(Y_list, ncol, integer(1)) == p)) {
    stop("All studies must have the same number of variables (same number of columns).", call. = FALSE)
  }
  invisible(p)
}

.new_bifa_fit <- function(backend, fit, meta = list()) {
  structure(list(backend = backend, fit = fit, meta = meta), class = c("bifa_fit", paste0("bifa_fit_", backend)))
}

.unwrap_fit <- function(x) {
  if (inherits(x, "bifa_fit")) return(x$fit)
  x
}

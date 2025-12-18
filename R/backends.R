#' Check which backends are installed and usable
#'
#' Many methods in the tutorial are implemented in external packages (often on GitHub
#' or Bioconductor). This helper reports availability and provides install hints.
#'
#' @return A data.frame with backend name and availability.
#' @export
available_backends <- function() {
  backends <- data.frame(
    backend = c("bmsfa", "stack_fa", "ind_fa", "momss", "sufa", "pfa", "tetris", "curated_ovarian"),
    requires = c("MSFA", "MSFA", "MSFA", "BFR.BE", "SUFA", "scripts", "scripts", "curatedOvarianData"),
    stringsAsFactors = FALSE
  )
  backends$available <- vapply(backends$backend, function(b) {
    if (b %in% c("bmsfa","stack_fa","ind_fa")) return(requireNamespace("MSFA", quietly = TRUE))
    if (b == "momss") return(requireNamespace("BFR.BE", quietly = TRUE))
    if (b == "sufa") return(requireNamespace("SUFA", quietly = TRUE))
    if (b == "curated_ovarian") return(requireNamespace("curatedOvarianData", quietly = TRUE))
    if (b %in% c("pfa","tetris")) return(nzchar(backend_script_dir(b)))
    FALSE
  }, logical(1))
  backends$hint <- vapply(backends$backend, function(b) {
    if (b %in% c("pfa","tetris")) return("Bundled in package under inst/extdata/")
    if (b == "curated_ovarian") return("BiocManager::install(\"curatedOvarianData\")")
    paste0("bifaToolkits::install_backend(\"", b, "\")")
  }, character(1))
  backends
}

#' Install optional backends used in the tutorial
#'
#' This is a convenience wrapper for non-CRAN dependencies.
#'
#' - For GitHub packages, we call `remotes::install_github()`.
#' - For Bioconductor packages, we call `BiocManager::install()`.
#' - For script-based methods (PFA / Tetris), scripts are bundled under `inst/extdata/`.
#'
#' @param backend One of: "bmsfa", "momss", "sufa", "curated_ovarian".
#' @param ... Passed to `remotes::install_github()` or `BiocManager::install()` when applicable.
#' @export
install_backend <- function(backend, ...) {
  backend <- tolower(backend)
  
  .require_pkg("remotes", "installing backends from GitHub")
  
  # safer download defaults (esp. inside RStudio)
  op <- options(timeout = 60, download.file.method = "libcurl")
  on.exit(options(op), add = TRUE)
  
  pat <- Sys.getenv("GITHUB_PAT", unset = NA_character_)
  
  if (backend %in% c("bmsfa", "stack_fa", "ind_fa")) {
    .require_pkg("remotes", "installing MSFA from GitHub")
    
    ok <- FALSE
    
    # Try GitHub API first
    ok <- tryCatch({
      remotes::install_github("Mavis-Liang/MSFA", ...)
      TRUE
    }, error = function(e) FALSE)
    
    # Fallback 1: git clone (no api.github.com)
    if (!ok) {
      ok <- tryCatch({
        remotes::install_git("https://github.com/Mavis-Liang/MSFA.git", ...)
        TRUE
      }, error = function(e) FALSE)
    }
    
    # Fallback 2: direct archive (no api.github.com, no git)
    if (!ok) {
      remotes::install_url("https://github.com/Mavis-Liang/MSFA/archive/refs/heads/master.tar.gz", ...)
    }
    
    return(invisible(TRUE))
  }
  
  if (backend == "momss") {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("sparseMatrixStats", ask = FALSE, update = FALSE)
    if (!requireNamespace("mombf", quietly = TRUE)) install.packages("mombf")
    
    remotes::install_github(
      "AleAviP/BFR.BE",
      auth_token = if (!is.na(pat)) pat else NULL,
      ...
    )
    return(invisible(TRUE))
  }
  
  if (backend == "sufa") {
    remotes::install_github(
      "noirritchandra/SUFA",
      build_vignettes = FALSE,
      auth_token = if (!is.na(pat)) pat else NULL,
      ...
    )
    return(invisible(TRUE))
  }
  
  if (backend == "curated_ovarian") {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("curatedOvarianData", ask = FALSE, update = FALSE)
    return(invisible(TRUE))
  }
  
  stop("Unknown backend: ", backend, call. = FALSE)
}

#' Get the installed directory for a script-based backend
#'
#' Script backends (PFA/Tetris) are bundled inside this package under `inst/extdata/`.
#' Use this helper to locate them after installation.
#'
#' @param backend "pfa" or "tetris".
#' @return Full path to the backend directory inside the installed package.
#' @export
backend_script_dir <- function(backend) {
  backend <- tolower(backend)
  
  # 1) installed package
  dir1 <- system.file("extdata", backend, package = "bifaToolkits")
  if (nzchar(dir1) && dir.exists(dir1)) {
    return(normalizePath(dir1, winslash = "/", mustWork = TRUE))
  }
  
  # 2) source tree (e.g., knitting vignettes from repo)
  root <- .find_pkg_root()
  if (nzchar(root)) {
    dir2 <- file.path(root, "inst", "extdata", backend)
    if (dir.exists(dir2)) {
      return(normalizePath(dir2, winslash = "/", mustWork = TRUE))
    }
  }
  
  ""
}


#' Load a bundled script backend into a private environment
#'
#' For PFA/Tetris, this sources all `.R` files under `inst/extdata/<backend>/`.
#' For PFA, it also compiles `PFA.cpp` (if present) via `Rcpp::sourceCpp()`.
#'
#' @param backend "pfa" or "tetris".
#' @return An environment containing the sourced functions.
#' @export
load_backend <- function(backend) {
  dir <- backend_script_dir(backend)
  if (!nzchar(dir) || !dir.exists(dir)) {
    stop(
      "Bundled backend directory not found for: ", backend, "\n",
      "Looked for installed extdata and for source inst/extdata/", backend, "/.\n",
      "Tip: if knitting the vignette from source, run from the package repo so DESCRIPTION is findable.",
      call. = FALSE
    )
  }
  
  env <- new.env(parent = globalenv())
  env$rgamma <- stats::rgamma
  env$rnorm  <- stats::rnorm
  env$runif  <- stats::runif
  env$rexp   <- stats::rexp
  
  r_files <- list.files(dir, pattern = "\\.R$", full.names = TRUE)
  if (length(r_files) == 0) {
    stop("No .R files found for backend '", backend, "' under: ", dir, call. = FALSE)
  }
  
  # ---- Source in stable order ----
  r_files <- normalizePath(r_files, winslash = "/", mustWork = TRUE)
  
  if (backend == "tetris") {
    if (!requireNamespace("clue", quietly = TRUE)) {
      stop("Tetris backend requires 'clue'. Install it with install.packages('clue').",
           call. = FALSE)
    }
    
    if (!requireNamespace("R.utils", quietly = TRUE)) {
      stop("Tetris backend requires 'R.utils'. Install it with install.packages('R.utils').",
           call. = FALSE)
    }
    
    if (!requireNamespace("combinat", quietly = TRUE)) {
      stop("Tetris backend requires 'combinat'. Install it with install.packages('combinat').",
           call. = FALSE)
    }
    
    for (f in sort(r_files)) source(f, local = env, chdir = TRUE)
    
  } else if (backend == "pfa") {
    if (!requireNamespace("expm", quietly = TRUE)) {
      stop("Tetris backend requires 'R.utils'. Install it with install.packages('R.utils').",
           call. = FALSE)
    }
    first  <- r_files[grepl("FBPFA-PFA\\.R$", r_files)]
    second <- r_files[grepl("FBPFA-PFA with fixed latent dim\\.R$", r_files)]
    
    for (f in c(first, second)) source(f, local = env, chdir = TRUE)
    
  } else {
    stop("Unknown script backend: ", backend, call. = FALSE)
  }
  
  env
}

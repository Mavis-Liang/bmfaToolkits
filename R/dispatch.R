#' Fit an integrative factor model with a chosen method within a general interface
#'
#' @param method One of: "bmsfa", "momss", "sufa", "pfa", "tetris", "stack_fa", "ind_fa", "cavi", "blast".
#' @param ... Passed to the corresponding `fit_*()` function.
#' @export
fit_integrative_fa <- function(method, ...) {
  method <- tolower(method)
  switch(
    method,
    bmsfa = fit_bmsfa(...),
    momss = fit_momss(...),
    sufa = fit_sufa(...),
    pfa = fit_pfa(...),
    tetris = fit_tetris(...),
    stack_fa = fit_stack_fa(...),
    ind_fa = fit_ind_fa(...),
    cavi = fit_cavi(...),
    blast = fit_blast(...),
    stop("Unknown method: ", method, call. = FALSE)
  )
}

#' Post-process an integrative factor model fit
#'
#' @param method Backend name.
#' @param fit A fit object from the corresponding `fit_*()` function.
#' @param ... Passed to the corresponding `postprocess_*()` function.
#' @export
postprocess_integrative_fa <- function(method, fit, ...) {
  method <- tolower(method)
  switch(
    method,
    bmsfa = postprocess_bmsfa(fit, ...),
    momss = postprocess_momss(fit, ...),
    sufa = postprocess_sufa(fit, ...),
    pfa = postprocess_pfa(fit, ...),
    tetris = postprocess_tetris(fit, ...),
    stack_fa = postprocess_stack_fa(fit, ...),
    ind_fa = postprocess_ind_fa(fit, ...),
    cavi = postprocess_cavi(fit, ...),
    blast = postprocess_blast(fit, ...),
    stop("Unknown method: ", method, call. = FALSE)
  )
}

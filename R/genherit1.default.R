## Method genherit1.default
##
##' Internal lmmtools function
##'
##' Default method for calculating the generalised heritability for a single site or trait.
##' The only available method is 'asreml'.
##' @title genherit1 for method 'default'
##' @param asr a model object
##' @param \ldots Additional parameter
##' @export
##'
genherit1.default <- function(asr, ...) {
    stop("Currently the only supported method is \"asreml\"")
}

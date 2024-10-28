##' Generic function for Generalised Heritability for a single trial
##'
##' Estimates the generalized heritability of one trial
##' @title generic genherit1
##' @param asr A model object.
##' @param id A character string giving the name of the genetic effect.
##' @param only A character vector specifying which variables form the prediction model. The default value is \code{id}.
##' @param Gmat A genetic variance or inverse variance matrix.  It must have a
##' logical attribute called \code{"INVERSE"}.  If the attribute is \code{TRUE},
##' \code{Gmat} is an inverse variance matrix.  If the attribute is \code{FALSE},
##' \code{Gmat} is a variance matrix.  The default for \code{Gmat} is \code{NULL}
##' in which case the variance matrix is assumed to be the scaled identity.
##' @param ... Other arguments to \code{asreml} can be given.
##' @return The estimation of generalised heritability (a number between 0 and 1)
##' @author Ari Verbyla \email{averbyla at avdataanalytics.com.au}
##'
##' @references
##' Oakey, H., Verbyla, A. P., Cullis, B. R., Wei, X., & Pitchford, W. S. (2007). Joint
##' modeling of additive and non-additive genetic line effects in multi-environment trials.
##' Theoretical and Applied Genetics, 114(8), 1319-1332.
##' @export
##'
genherit1 <- function(asr, ...) {
    UseMethod("genherit1")
}

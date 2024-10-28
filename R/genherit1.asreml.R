##' Generalised Heritability for a single trial
##'
##' Estimates the generalized heritability of one trial
##' @title genherit1 for method 'asreml'
##' @param asr An \code{asreml} model object.
##' @param id A character string giving the name of the genetic effect.
##' @param only A character vector specifying which variables form the prediction model. The default value is \code{id}.
##' @param Gmat A genetic variance or inverse variance matrix.  It must have a
##' logical attribute called \code{"INVERSE"}.  If the attribute is \code{TRUE},
##' \code{Gmat} is an inverse variance matrix.  If the attribute is \code{FALSE},
##' \code{Gmat} is a variance matrix.  The default for \code{Gmat} is \code{NULL}
##' in which case the variance matrix is assumed to be the scaled identity.
##' @param ... Other arguments to \code{asreml} can be given.
##' @export
##' @return The estimation of generalised heritability (a number between 0 and 1)
##' @author Ari Verbyla \email{averbyla at avdataanalytics.com.au}
##'
##' @references
##' Oakey, H., Verbyla, A. P., Cullis, B. R., Wei, X., & Pitchford, W. S. (2007). Joint
##' modeling of additive and non-additive genetic line effects in multi-environment trials.
##' Theoretical and Applied Genetics, 114(8), 1319-1332.
##' @export
##'
genherit1.asreml <- function(asr, id = 'Genotype', only = NULL, Gmat = NULL, ...) {
    response <- as.character(asr$call[[2]][2])
    if(is.null(only)) {
        if(!is.null(Gmat)) stop("\n Error:  argument only is NULL, but Gmat is not NULL\n")
        only <- id
    }
    if(packageVersion("asreml") == "3.0") {
        mypred <- predict(asr, classify=id, maxiter=1, only=only, vcov=TRUE, ...)$predictions
        which.vc <- grep(id, names(asr$gammas))
    } else {
        mypred <- predict(asr, classify=id, maxit=1, only=only, vcov=TRUE, ...)
        which.vc <- grep(id, names(asr$vparameters), fixed=TRUE)
    }
    my.vcov <- mypred$vcov
    if(is.null(Gmat)) {
##################  version 3 might require summary(asr)$varcomp[which.vc,2]
        h2 <- 1- sum(diag(my.vcov))/(summary(asr)$varcomp[which.vc,1]*nrow(my.vcov))
    }  else {
        if(is.null(attr(Gmat, "INVERSE"))) stop('\n Gmat must have an "INVERSE" attribute\n')
        if(!attr(Gmat, "INVERSE")) Gmat <- ginv(Gmat)
        h2 <- 1- sum(diag(Gmat%*%my.vcov))/nrow(my.vcov)
    }
    names(h2) <- response
    h2
}

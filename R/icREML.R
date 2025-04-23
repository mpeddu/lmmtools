## Function icREML
##
##' icREML
##'

##' This function calculates the AIC and BIC for each fitted asreml model in a
##' list.  Options are set and each model is updated.  The elements required
##' for the AIC and BIC are then calculated.
##'
##' @title Find the AIC and BIC for a \code{list} of models fitted using
##' \code{asreml}
##' @param fm A \code{list} of asreml fitted model objects
##' @param scale A scalar to scale the variance matrix of the estimated
##' fixed effects (to ensure numerical stability of a log-determinant).
##' Default value is 1.
##' @param logdet The log determinant that modifies the residual
##' log-likelihood to form the full log-likelihood is to be included in the
##' table.
##' @return A data frame.  The data frame has the following components
##' \itemize{
##' \item \code{model} : the names of the models
##' \item \code{res.loglik} : the residual log-likelihood for each model
##' \item \code{full.loglik} : the full log-likelihood for each model
##' \item \code{p} :  the number of fixed effects parameters for each model
##' \item \code{q} : the number of (non-zero) variance parameters for each model.
##' \item \code{b} : the number of variance parameters that are fixed or on the
##' boundary.  These parameters are not counted in the AIC or BIC.
##' \item \code{AIC} : the AIC for each model
##' \item \code{BIC} : the BIC for each model
##' \item \code{logdet} : the log-determinant used in adjusting the residual
##' log-likelihood for each model
##' }
##'
##' @importFrom stringr str_extract
##'
##' @export
##'
icREML <- function(fm, scale=1, logdet=FALSE) {
    stopifnot(requireNamespace("asreml"))

    if(length(fm) == 0) stop(" The model list is empty\n")
    if(!is.list(fm)) stop(" Models need to be in a list\n")
    if(is.null(names(fm))) namesfm <- paste0("fm", 1:length(fm))
    else namesfm <- names(fm)
    fm <- lapply(1:length(fm), function(el, fm) {
        if(is.null(fm[[el]]$Cfixed)) {
            cat(" Updating model ", names(fm)[el], " for likelihood calculation\n")
            myfm <- fm[[el]]
            if(packageVersion("asreml") >= "4.0") {
                asr.opt <- asreml::asreml.options()
                asreml::asreml.options(Cfixed = TRUE, gammaPar=FALSE)
                out <- asreml::update(myfm, maxit=1)
                asreml::asreml.options(asr.opt)
            }
            else out <- asreml::update(myfm, maxit=1, Cfixed=TRUE)
        }
        else {
##            print(el)
##            print(names(fm)[el])
            out <- fm[[el]]
        }
        out}, fm=fm)
    logl <- lapply(fm, function(el) el$loglik)
    summ <- lapply(fm, function(el) asreml::summary.asreml(el, coef=TRUE)$coef.fixed)
    which.X0 <- lapply(summ, function(el) !is.na(el[, "z.ratio"]))
    p.0 <- lapply(which.X0, function(el) sum(el))
    Cfixed <- lapply(fm, function(el) {
        if(!is.null(dim(el))) {
#            ord <- rownames(summary(el, coef=TRUE)$coef.fixed)
#            el$Cfixed <- el$Cfixed[ord, ord, drop=FALSE]
            summv <- asreml::summary.asreml(el)$varcomp
            uR <- grep("units", dimnames(summv)[[1]])
            if(length(uR) != 0) {
                if (summv$bound[uR] == "F")
                    el$Cfixed  <- el$Cfixed/el$sigma2
            }
        }
        el$Cfixed
    })
    logdetC <- lapply(1:length(fm), function(el, Cfixed, which.X0, scale) {
        if(!is.null(dim(Cfixed[[el]]))) {
            mysvdd <- svd(scale*Cfixed[[el]][which.X0[[el]], which.X0[[el]]])$d
            mysvdd <- mysvdd[mysvdd > 0]
            ldc <- sum(log(mysvdd))
        }
        else ldc <- log(scale*Cfixed[[el]])
        ldc},
        Cfixed, which.X0, scale)
    vparam <- lapply(fm, function(el) asreml::summary.asreml(el)$varcomp)
    q.0 <- lapply(vparam, function(el) sum(!(el$bound == "F" | el$bound == "B" | el$bound == "C")) + sum(el$bound[!is.na(stringr::str_extract(dimnames(el)[[1]], "cor"))] == "B"))
    b.0 <- lapply(vparam, function(el) sum(el$bound == "F" | el$bound == "B") - sum(el$bound[!is.na(stringr::str_extract(dimnames(el)[[1]], "cor"))] == "B"))
    full.logl <- lapply(1:length(fm), function(el, logl, logdetC, p.0) {
        logl[[el]] - logdetC[[el]]/2}, logl, logdetC, p.0)
    aic <- unlist(lapply(1:length(fm), function(el, full.logl, p.0, q.0) {
        -2*full.logl[[el]] + 2*(p.0[[el]] + q.0[[el]])}, full.logl, p.0, q.0))
    bic <- unlist(lapply(1:length(fm), function(el, full.logl, p.0, q.0, fm) {
        -2*full.logl[[el]] + log(fm[[el]]$nedf+p.0[[el]])*(p.0[[el]] + q.0[[el]])},
        full.logl, p.0, q.0, fm))
    results <- data.frame(model=namesfm, res.loglik = unlist(logl), full.loglik = unlist(full.logl), p=unlist(p.0),
                          q=unlist(q.0), b = unlist(b.0), AIC = aic, BIC = bic)
    if(logdet) results$logdet <- unlist(logdetC)
    row.names(results) <- 1:dim(results)[1]
    results
}

## faSummary function
##
##' faSummary
##'
##' For a list of models from a diagonal form (always assumed first in the list),
##' to a sequential set of factor analysis (FA) models, tabulate various statistics
##' to provide guidance on model choice of an FA model.
##' @title faSummary
##' @param fm a list of fitted models, from a diagonal form to factor analytic models
##' of increasing order
##' @param data the data frame used in fitting the models
##' @param Trait a character string giving the names of the trait factor (for example site)
##' @param id a character string giving the name of the genotype/variety factor.  Default is
##' "id"
##' @return a data frame with various statistics and diagnostics.
##' @export
##'
faSummary <- function(fm, data, Trait, id = "id") {
###
###  fm is a list of an increasing sequence of fa models:  the first is assumed to be diag
###  data is the data frame of the data used to fit the models in fm
###  Trait is a character that is the Trait in the fa model
###
    n.mod <- length(fm)
    pp <- length(levels(data[, Trait]))
    n <- dim(data)[1]
    loglik <- unlist(lapply(fm, function(el) el$loglik))
    lrt <- c(NA, 2*diff(loglik))
    summ.bnd <- lapply(fm, function(el) asreml::summary.asreml(el)$varcomp[, "bound"])
    np <- unlist(lapply(summ.bnd, function(el) sum(el != "B" & el != "F")))
    lrt.df  <- c(NA, diff(np))
    lrt.pval <- c(NA, 1-pchisq(lrt[-1], df=lrt.df[-1]))
    aic <- round(-2*loglik + 2*np,2)
    bic <- round(-2*loglik + np*log(n),2)
    vm.lst <- list()
    perc.var <- c()
    fn.val <- fn.cor <- norm.var <- meyer <- rep(NA, n.mod)
    for (ii in 1:n.mod) {
        el <- fm[[ii]]
        grp <- grep(paste0(":", id), names(el$vparameters))
        faval <- el$vparameters[grp]
        n.fa <- (length(faval)/pp)-1
        psi <- faval[1:pp]
        varmat <- diag(psi)
        perc.var[ii] <- NA
        if(ii > 1) {
            Lambda <- matrix(faval[(pp+1):length(faval)], ncol=n.fa)
            LLt <- Lambda %*% t(Lambda)
            varmat <- LLt + diag(psi)
            perc.var[ii] <- sum(diag(LLt))/sum(diag(varmat))
        }
        vm.lst[[ii]] <- varmat
        meyer[ii] <- sum(diag(varmat))
        if(ii > 1) {
            sdii <- diag(1/sqrt(diag(vm.lst[[ii]])))
            sdii1 <- diag(1/sqrt(diag(vm.lst[[ii-1]])))
            corii <- sdii %*% vm.lst[[ii]] %*% sdii
            corii1 <- sdii1 %*% vm.lst[[ii-1]] %*% sdii1
            vm.diff <- vm.lst[[ii]]-vm.lst[[ii-1]]
            cor.diff <- corii-corii1
            fn.val[ii] <- Frobenius(vm.diff)/Frobenius(vm.lst[[ii-1]])
            fn.cor[ii] <- Frobenius(cor.diff)/Frobenius(corii1)
            norm.var[ii] <- sqrt(sum(diag(vm.diff)^2))/sqrt(sum(diag(vm.lst[[ii-1]])^2))
        }
    }
    names(loglik) <- names(lrt) <- names(lrt.df) <- names(lrt.pval) <- names(aic) <- names(bic) <- names(vm.lst) <- names(perc.var) <- names(fn.val) <- names(fn.cor) <- names(norm.var) <- names(fm)
    fa.df <- data.frame(model = names(fm), loglik=round(loglik,2), seq.lrt=round(lrt,2), df=lrt.df, pval = round(lrt.pval,4), AIC = aic, BIC = bic, perc.var=round(perc.var*100,2), fn.vm = round(fn.val*100,2), meyer = round(meyer,2), norm.var=round(norm.var*100,2), fn.cor=round(fn.cor*100,2))
    fa.df
}


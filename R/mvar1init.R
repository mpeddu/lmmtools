## mvar1init function
##
##' mvar1 initial value function: mvar1init
##'
##' This function sets up starting values for the parameters
##' of a multivariate autoregressive process of order 1
##' @title mvar1 initial values
##' @param gam a vparameters.table object from an \code{asreml}
##' call using \code{start.values=TRUE}
##' @param init.var initial variances for each level of the factor
##' given by the \code{Traits} argument with a names attribute that
##' matches \code{Traits}(see below)
##' @param spatial a list of character strings of the model for the spatial
##' effects
##' @param Traits a vector of character strings giving the name of the traits
##' that defines the mvar1 structure in addition to rows and columns
##' @param n.trait the number of levels of the trait factor that
##' should match the length of the initial variances specified
##' in the argument \code{init.var}
##' @return a vparamters.table object with initial values of the
##' mvar1 parameters set
##' @export
##'
mvar1init <-
function(gam, init.var, spatial=list("Run", "Range"), Traits=NULL, n.trait=length(init.var)) {
    if(is.null(Traits)) stop(" You need to supply the character vector of trait names for the MVAR1 terms\n")
    if(n.trait == 0) stop(" No initial values in init.var for fa terms\n")
    if(n.trait != length(Traits)) stop(" Length of initial values does not match the length of the Traits vector\n")
    cat("\n All fa1 specific variances set to 1e-6 and fixed as asreml can't handle zero variances\n\n")
    for (el in Traits) {
        tlist <- c(el, spatial)
        cat(" Term with ", unlist(tlist), "\n")
        which.var  <- lgrep(c(tlist, "var"), gam$Component)
        gam$Value[which.var] <- 1e-06
        gam$Constraint[which.var] <- "F"
        which.fa1 <- lgrep(c(tlist, "fa1"), gam$Component)
        names(which.fa1) <- names(init.var)
        gam$Value[which.fa1] <- sqrt(init.var[el])/100
        gam$Value[which.fa1[el]] <- sqrt(init.var[el])
        cat(" fa1 initial values ", gam$Value[which.fa1], "\n")
    }
    gam
}

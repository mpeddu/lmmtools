## Put together all variance parameters in stage 1 analysis
##
## covMat
##
##' Collate all variance parameters in a stage 1 analysis
##'
##' Stage 1 analysis variance parameters
##'
##' @title Stage 1 analysis collation of variance parameters
##' @param vc a list of varcomp summary objects from each of the stage 1
##' analyses using \code{asreml}
##' @param which.comb A matrix of combinations of levels of the factor
##' used to subset the stage 1 analyses
##' @param ltrait the levels of the trait indexing the subsetting in
##' the stage 1 analyses
##'
##' @importFrom stringr str_sort
##'
##' @export
##'
covMat  <- function(vc, which.comb, ltrait) {
    tmp <- unlist(lapply(vc, function(el) dimnames(el)[[1]]))
    tmp.unique <- unique(tmp)
    dim1nam <- stringr::str_sort(tmp.unique, numeric=TRUE)
    tmp <- matrix(ltrait[which.comb], dim(which.comb))
    dim2nam <- apply(tmp, 1, function(el) paste(as.character(el), collapse=","))
    d1 <- length(dim1nam)
    d2 <- length(dim2nam)
    covMat <- matrix(NA, d1, d2)
    dimnames(covMat) <- list(dim1nam, dim2nam)
    for (ii in 1:d1) {
        which.term <- dim1nam[ii]
        covMat[ii, ]  <- unlist(sapply(vc, function(el, which.term) {
            test <- dimnames(el)[[1]] == which.term
            if(any(test)) out <- el$component[test]
            else out <- NA
            out}, which.term=which.term))
    }
    covMat
}

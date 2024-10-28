## Function Frobenius
##
##' Frobenius norm
##'
##' This function calculates a distance measure using the Frobenius norm
##' @title Frobenius Norm
##' @param D a matrix
##' @return the Frobenius norm, a scalar
##' @author Ari Verbyla (averbyla at avdataanalytics.com.au)
##' @export
##'
Frobenius <- function(D) {
    DtD <- t(D) %*% D
    F2 <-  sum(diag(DtD))
    sqrt(F2)
}

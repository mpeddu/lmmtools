##
## Take a symmetric matrix and form a three column (matrix) sparse form
##
##' mat2df
##'
##' Convert a symmetric matrix to a three column (matrix) sparse form
##'
##' @title Conversion of a symmetric matrix to a three column (matrix) sparse form
##' @param mat A symmetric matrix.
##' @return A data frame with row, column and values from the lower triangle
##' of the matrix \code{mat}
##' @author Ari Verbyla (averbyla at avdataanalytics.com.au)
##' @export
##'
mat2df  <- function(mat) {
    if(diff(dim(mat)) != 0) stop(" Error: The matrix is not square\n")
    if(!all(mat==t(mat))) stop(" Error: The matrix is not symmetric\n")
    val  <- mat[lower.tri(mat, diag=TRUE)]
    kk <- 1:dim(mat)[2]
    mygrid <- expand.grid(kk,kk)
    mygrid <- as.matrix(mygrid[mygrid$Var1 >= mygrid$Var2,])
    spmat  <- cbind.data.frame(mygrid, val)
    names(spmat) <- c("Row", "Col", "Value")
    spmat
}



##' mat2sparse
##'
##' Convert a symmetric matrix to a three column (matrix) sparse form
##'
##' @param mat A symmetric matrix.
##' @return A matrix with row, column and values from the lower triangle
##' of the matrix \code{mat}, which is row-major order required by
##' \code{asreml}.
##' @export
##'
mat2sparse  <- function(mat) {
    if(diff(dim(mat)) != 0) stop(" Error: The matrix is not square\n")
    if(!all(mat==t(mat))) stop(" Error: The matrix is not symmetric\n")
    val  <- mat[lower.tri(mat, diag=TRUE)]
    kk <- 1:dim(mat)[2]
    mygrid <- expand.grid(kk,kk)
    mygrid <- as.matrix(mygrid[mygrid$Var1 <= mygrid$Var2,])
    spmat  <- cbind(mygrid, val)
    dimnames(spmat) <- list(1:dim(spmat)[1], c("Row", "Col", "Value"))
    spmat
}



## wheat data documentation
##
##' WheatData: variety field trial
##'
##' The field trial was almost a complete randomised block design,
##' consisted of 22 rows by 15 columns, with three blocks of 5 columns.
##' The 107 varieties were replicated three times (with some standards
##' replicated twice in each block).  The trait of interest was grain
##' yield (grams per hectare).  The aim of the analysis was variety selection.
##' The variables are as follows:
##'
##' \itemize{
##' \item yield. The grain yield in grams per hectare.
##' \item Variety.  The variety (as a factor) in a particular plot, as
##' defined by \code{Row} and \code{Col}, see below.
##' \item Block.  The Block in which the plot occurs.  A factor.
##' \item row.  A variate for the row a plot is in.
##' \item col.  A variate for the column a plot is in.
##' \item Rowcode.  A factor required for some analyses.
##' \item Colcode.  A factor required for some analyses.
##' \item Row.  A factor for the row a plot is in.
##' \item Col.  A factor for the column a plot is in.
##' }
##'
##' @references{Gilmour et al. 1997, Journal of Agricultural, Biological
##' and Environmental Statistics 2, 269-293.}
##' @name WheatData
##' @docType data
NULL

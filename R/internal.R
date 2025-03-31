##  merge formulae into a single formula
##
##' Merge two formulae into a single formula
##'
##' The function takes two formulae and combines them into a single formula
##'
##' @title Merge two formulae
##' @param form1 a formula
##' @param form2 a formula
##' @param LHS a logical scalar that indicates if there is a response in \code{form1}
##' and \code{form2}
##' @param ... additional arguments
##' @return A single model formula
##' @export
##'
merge_formula <- function(form1, form2, LHS=TRUE, ...){

                                        # get character strings of the names for the responses
                                        # (i.e. left hand sides, lhs)
    if(LHS) {
        lhs1 <- deparse(form1[[2]])
        lhs2 <- deparse(form2[[2]])
        if(lhs1 != lhs2) stop('both formulas must have the same response')
        lhs <- lhs1
        which.one <- 3
    }
    else {
        lhs <- NULL
        which.one <- 2
    }

                                        # get character strings of the right hand sides
    rhs1 <- strsplit(deparse(form1[[which.one]]), " \\+ ")[[1]]
    rhs2 <- strsplit(deparse(form2[[which.one]]), " \\+ ")[[1]]

                                        # create the merged rhs and lhs in character string form
    rhs <- c(rhs1, rhs2)

                                        # put the two sides together with the amazing
                                        # reformulate function
    out <- reformulate(rhs, lhs)

                                        # set the environment of the formula (i.e. where should
                                        # R look for variables when data aren't specified?)
    environment(out) <- parent.frame()

    return(out)
}


#' @rdname lgrep
#' @export
igrep <- function(patterns, x, ...){
    if(length(patterns) == 1)
        grep(patterns[[1]], x, ...)
    else {
        xt <- x
        for(i in 1:length(patterns)){
            ind <- grep(patterns[[i]], x, ...)
            x <- x[ind]
        }
        pmatch(x, xt)
    }
}

#' Search for patterns in lists
#'
#' @param patterns list of character vectors.
#' @param x object for which patterns are sought.
#' @param ... Additional parameters.
#' @return \code{lgrep}: The index of matches of the intersection of patterns in \code{x}.
#' \code{igrep}: The index of matches of the patterns in \code{x}.
#' @export
lgrep <- function(patterns, x, ...){
    if(length(patterns) == 1)
        grep(patterns[[1]], x, ...)
    else {
        ind <- list()
        for(ii in 1:length(patterns)){
            ind[[ii]] <- grep(patterns[[i]], x, ...)
        }
    Reduce(intersect, ind)
    }
}


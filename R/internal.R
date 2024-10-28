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
##' @param \code{...} additional arguments
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


## An addition operator that works for formulae
##
##' An addition operator that works for formulae
##'
##' @title Addition operator for formulae
##' @param e1 and e2 Two model formulae
##' @return The sum of the two model formulae
##' @export
##'
Ops.formula <- function(e1, e2){
	FUN <- .Generic
	if(FUN == '+'){
		out <- merge(e1, e2)
		environment(out) <- parent.frame()
		return(out)
	}
	else stop('can not yet subtract formula objects')
}

##'
##' igrep(patterns, x, ...)
##' lgrep(patterns, x, ...)
##' @title Internal splinetools function.
##' @param patterns list of character vectors.
##' @param x object for which patterns are sought.
##' @param \ldots Additional parameters.
##' @return The index of matches of patterns in x.
##' @author Julian Taylor <julian.taylor@adelaide.edu.au>
##' @export
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

##' lmmtools internal functions
##'
##' These functions are used internally in lmmtools.
##'
##' igrep(patterns, x, ...)
##' lgrep(patterns, x, ...)
##' @title Internal lmmtools function.
##' @param patterns: list of character vectors.
##' @param x: object for which patterns are sought.
##' @param ... Additional parameters.
##' @return The index of matches of the intersection of patterns in x.
##' @author Ari Verbyla (averbyla at avdataanalytics.com.au)
##' @export
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


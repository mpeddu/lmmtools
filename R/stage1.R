## Stage 1 analysis function
##
##' A function for stage 1 analyses of a two stage analysis
##'
##' A function that calls asreml for each component of a stage 1
##' analysis.
##' @title Stage 1 analysis function using asreml
##' @param fixed A model formula for the fixed effects passed to
##'     asreml
##' @param random A model formula for the random effects passed to
##'     asreml
##' @param data  The data frame or data table containing the data to
##'      be used in analysis.
##' @param residual A model formula for the residual effects passed to
##'     asreml
##' @param sparse A model formula for the sparse effects passed to
##'     asreml
##' @param family A family call passed to asreml
##' @param specific A list of terms to be added to the model that are
##'     specific to Trait values.  Each component of the list is
##'     itself a list with the type of effect to be included as a
##'     component (fixed or random), and a vector specifying which
##'     Trait values require the term as another component.  This
##'     component can have names (character vector) or be numeric.  In
##'     the latter case the values are assumed to correspond to the
##'     order of the Trait factor levels.
##' @param Genetic A character vector specifying the name of the
##'     genetic effect that is being fitted, The default is "id".
##' @param predGenetic A list specifying the \code{id} argument which
##'     gives the name of the genetic term to be fitted at stage 2
##'     (default, the same as the \code{Genetic} argument, the
##'     \code{classify} term for prediction of genetic effects
##'     (default, the same as the \code{Genetic} argument), and any
##'     terms to remove, called \code{rm} in the list, from the model
##'     for prediction of fixed genetic effects.  The default for the
##'     classify is set to the value of the \code{Genetic} argument by
##'     default, while \code{rm} is initially set to \code{NULL}
##' @param predict.options A list of options for the prediction step
##'     in the stage 1 analyses, that are applicable to all stage 1
##'     analyses.
##' @param type A character vector specifying the type of interaction
##'     between Trait and Genetic.  The default is "random".  The
##'     alternative is "fixed".
##' @param Trait A character string that defines the trait involved.
##' @param n.trait The number of levels of Trait to be used in the
##'     first stage analysis.  The default is 1.
##' @param keep.models A logical that flags whether all the models
##'     fitted at stage 1 should be saved as part of the output list.
##'     If \code{type} is random, the returned models have random
##'     Genetic effects, otherwise fixed Genetic effects.
##' @param mvar1.init.obj A list of size 3.  The component \code{val}
##'     is a vector of initial values (variances) for
##'     each of the levels of Trait.  The other two components are
##'     \code{Run} and \code{Range} which are characters that specify
##'     the factors for the field Run and Range.  If NULL, the
##'     multivariate ar1 is not being fitted and this is the default.
##'     This structure will only be used for repeated measures data.
##' @param s1trace A logical scalar indicating if trace printing is to
##'     be switched on for debugging
##' @param ... Further arguments to be used in the underlying
##'     \code{asreml} fit, for example \code{workspace}, \code{maxit}
##'     that are standard asreml arguments.
##' @return A list of named components.  \code{call} contains the call
##'     to this function. \code{Vt} contains a three column data frame
##'     of the variance matrix of the first stage effects, row,
##'     column, and value.  It has two attributes, namely
##'     \code{"INVERSE"} which is \code{FALSE}, and \code{"rowNames"}
##'     required for \code{vm} in the second stage analysis in
##'     \code{asreml}.  \code{bdVt} is a sparse block diagonal
##'     variance matrix for use in the second stage analysis using
##'     \code{vm} in package \code{asreml}. \code{V} and \code{W} are
##'     lists of the variance and weight matrices from the stage 1
##'     analyses. \code{pred.df} is a data frame with Trait by genetic
##'     effects for use in the second stage analysis.  \code{type}
##'     records the type of first stage analysis used (fixed or
##'     random).  \code{vc} is a list of estimated variance parameters
##'     for each stage 1 analysis.  \code{vcMat} is a list of matrices
##'     of variance parameters for all stage 1 models. \code{models}
##'     is a list of the fitted model objects from the stage 1 fits if
##'     \code{keep.models} has been set to \code{TRUE}.
##'
##' @importFrom methods as
##' @importFrom stats formula reformulate
##' @importFrom utils head str
##' @importFrom gtools combinations
##' @importFrom MASS ginv
##'
##' @export
##'
stage1 <- function(fixed, random, data, residual=NULL, sparse=NULL, family=NULL,
                   specific = list(), Genetic="id",
                   predGenetic = list(id = Genetic, classify=Genetic, rm=NULL),
                   predict.options = NULL, type="random", Trait="Trait", n.trait=1,
                   keep.models = FALSE, mvar1.init.obj =NULL, s1trace=FALSE, ...) {

    stopifnot(requireNamespace("asreml"))

    call <- match.call()
    asr.options <- asreml::asreml.options()
    ltrait <- levels(data[, Trait])
    l.Trait <- length(ltrait)
    asreml::asreml.options(Cfixed=TRUE, extra=3, ai.sing=TRUE)
    V <- W <- pred.lst <- vc <- mymodels <- list()
    sigma2 <- c()
    subnames <- character()
    which.comb <- gtools::combinations(l.Trait, n.trait)
    predGen  <- predGenetic$id
    fremove <- vector(mode="list", length=l.Trait)
    if(!is.null(predGenetic$rm)) {
        for (ii in 1:length(predGenetic$rm)) {
            rmii <- predGenetic$rm[[ii]]
            fremove[rmii$which] <- lapply(rmii$which,
                                          function(el, rmii, fremove) append(fremove[[el]], rmii$terms),
                                          rmii, fremove)
        }
    }
    if(s1trace) print(fremove)
    if(n.trait > 1) predGen <- paste(Trait, predGen, sep=":")
    ddd <- list(...)
    for (ii in 1:dim(which.comb)[1]) {
        which.ii <- which.comb[ii, ]
        which.l <- ltrait[which.ii]
        if(n.trait == 1) subnames[ii] <- as.character(which.l)
        else subnames[ii] <-  paste(as.character(which.l), collapse=",")
        cat(Trait, " value(s) ", which.l, "\n")
        data.l <- subset(data, data[, Trait] %in% which.l)
        data.l <- droplevels(data.l)
        if(s1trace) print(str(data.l))
        lev.gen <- levels(data.l[, Genetic])
        l.gen <- length(lev.gen)

######################  Fix

        tmp <- unlist(lapply(specific, function(el, which.l, ltrait) {
            if(!is.numeric(el$which)) ww <- el$which %in% which.l
            else ww <- ltrait[el$which] %in% which.l
            any(ww, na.rm=TRUE)}, which.l=which.l, ltrait=ltrait))


        afixed <- fixed
        arandom <- random


        if(any(tmp)) {
            this.spec <- specific[tmp]
            which.type <- unlist(lapply(this.spec, function(el) el$type))
            if(any(which.type == "fixed")) {
                which.fixed <- names(this.spec)[which.type == "fixed"]
                add.fixed <- formula(paste0(as.character(fixed[[2]]), " ~ ",
                                            paste(which.fixed, collapse = " + ")))
                afixed <- merge_formula(afixed, add.fixed)
            }
            if(any(which.type == "random")) {
                which.random <- names(this.spec)[which.type == "random"]
                add.random <- formula(paste0(" ~ ", paste(which.random, collapse = " + ")))
                arandom <- merge_formula(arandom, add.random, LHS=FALSE)
            }
        }
###################################

        if(s1trace) {
            print(l.gen)
            print(afixed)
            print(arandom)
        }
        mycall <- list()
        mycall[[1]] <- as.name("asreml")
        mycall$fixed <- afixed
        mycall$random <- arandom
        mycall$residual <- residual
        mycall$sparse <- sparse
        mycall$family <- family
        mycall$data <- data.l
        if(is.null(mycall$maxit)) mycall$maxit <- 100
        mycall <- c(mycall, ddd)
        mycall <- as(mycall, "call", strict=TRUE)
        if(n.trait > 1) {
            if(!is.null(mvar1.init.obj)) {
                Traits <- paste0(Trait, which.l)
                if(s1trace) print(Traits)
                mycall$data[, Traits] <- mycall$data[, Trait]
                mycall$data[, Traits] <-
                    lapply(mycall$data[, Traits],
                           function(el, which.l)
                               ordered(el, levels=which.l), which.l)
                str(mycall$data)
                mycall$start.values <- TRUE
                Run <- mvar1.init.obj$Run
                Range <- mvar1.init.obj$Range
                mvar1.trm <- unlist(lapply(Traits, function(el, Run, Range) paste0("fa(", el, ", 1):ar1(", Run, "):ar1(", Range, ")"), Run=Run, Range=Range))
                if(s1trace) print(mvar1.trm)
                rhs1 <- strsplit(deparse(arandom[[2]]), " \\+ ")[[1]]
                if(s1trace) print(rhs1)
                rhs <- c(rhs1, mvar1.trm)
                if(s1trace) print(rhs)
                aranmvar1 <- reformulate(rhs, NULL)
                if(s1trace) print(arandom)
                if(s1trace) print(aranmvar1)
                mycall$random <- aranmvar1
                if(s1trace) print(mycall$random)
                sv <- eval(mycall)
                gam <- sv$vparameters.table
                gam <- mvar1init(gam, mvar1.init,obj$val[which.ii], Traits = Traits, n.trait = which.l)
                if(s1trace) print(gam)
                mycall$G.param <- gam
                mycall$start.values <- FALSE
            }
        }
        asm <- eval(mycall)
        summ <- summary(asm)$varcomp
        mymodels[[ii]] <- asm
        if(s1trace) print(summ)
        frhs <- strsplit(deparse(afixed[[3]]), split=" \\+ ")[[1]]
        if(s1trace) print(frhs)
        if(!is.null(fremove[[ii]])) {
            frhs <- frhs[!(frhs %in% fremove[[ii]])]
        }
        lhs <- deparse(afixed[[2]])
        if(type == "random") frhs <- c(frhs, predGen)
        afixed <- reformulate(frhs, lhs)
        mycall$fixed <- afixed
##            ffixed <- formula(paste(" . ~ . + ", predGen))
        if(s1trace) print(afixed)
#####
#####  Need to get the random term with Genetic in it.  It might be a compound structure
#####
        if(type=="random") {
            mycall$start.values <- TRUE
            mycall$G.param <- asm$G.param
            mycall$R.param <- asm$R.param
            sv <- eval(mycall)
            gam <- sv$vparameters.table
            gam[, 3] <- "F"
            mycall$G.param <- gam
            mycall$R.param <- gam
            mycall$start.values <- FALSE
            ar.terms <- strsplit(deparse(arandom[[2]]), split=" \\+ ")[[1]]
            which.aterm <- grep(Genetic, ar.terms)
            ar.terms <- ar.terms[- which.aterm]
            arandom <- reformulate(ar.terms, NULL)
            mycall$random <- arandom
##            rm.terms <- paste(ar.terms[which.aterm], collapse = "-")
##            rrandom <- formula(paste(" ~ . - ", rm.terms))
            if(s1trace) print(arandom)
            cat("  De-regression\n")
            asm <- eval(mycall) ##update(asm, fixed. = ffixed, random. = rrandom, G.param=gam, R.param=gam)
        }
        sigma2[ii] <- asm$sigma2
        summ <- summary(asm)$varcomp
        if(s1trace) print(summ)
        uR <- grep("units!R", dimnames(summ)[[1]])
        if(length(uR) != 0) {
            if (summ$bound[uR] == "F")
                sigma2 <- 1
        }
        ## if(s1trace) print(asm$call)
##        pred <- predict(asm, classify = predGen,
##                        vcov=TRUE, maxit=1,  ...)
        mypred <- list()
        mypred[[1]] <- as.name("predict.asreml")
        mypred$object <- quote(asm)
        mypred$classify <- predGenetic$classify
        mypred$vcov <- TRUE
        mypred$maxit <- 1
        if(!is.null(predict.options)) {
            if(s1trace) print(predict.options)
            mypred <- c(mypred, predict.options)
        }
        if(s1trace) print(mypred)
        mypred <- as(mypred, "call", strict=TRUE)
        cat(" Prediction step\n")
        pred <- eval(mypred)
        if(s1trace) print(names(pred))
        if(s1trace) print(head(pred$pvals))
        pev <- pred$vcov
        if(s1trace) print(pev[1:5,1:5])
##        if(length(predGenetic) > 0) {
##            lev.predgen <- levels(pred$pvals[, predGen])
##            which.gen <- pred$pvals[, predGen] %in% lev.gen
##            pred$pvals <- pred$pvals[which.gen, ]
##            pred$pvals <- droplevels(pred$pvals)
##            pev <- pev[which.gen, which.gen]
##        }
        if(n.trait == 1) {
            nam <- as.character(levels(data.l[, predGen]))
            trt <- rep(ltrait[ii], length(levels(data.l[, predGen])))
        }
        else {
            nam <- paste(as.character(levels(data.l[, Trait])),
                         as.character(levels(data.l[, predGen])), sep="\\.")
        }
        pev2 <- asm$Cfixed
        if(s1trace) print(dim(pev))
        which.omit <- is.na(pred$pvals$predicted.value)
        pvals <- pred$pvals[!which.omit,]
        pev <- pev[!which.omit, !which.omit]
        if(s1trace) print(dim(pev))
        if(s1trace) print(pev[1:4,1:4])
        if(s1trace) print(pev2[1:4,1:4])
        pevi <- MASS::ginv(as.matrix(pev))
        if(s1trace) print(pevi[1:4,1:4])
        Vt <- as.matrix(pev)
        nam <- nam[!which.omit]
        dimnames(Vt) <- list(nam, nam)
        Vt.mat <- mat2sparse(Vt)
        if(ii == 1) Vtfull.mat <- Vt.mat
        else {
            curr.dim <- max(Vtfull.mat[, "Row"])
            Vt.mat[, "Row"] <- Vt.mat[, "Row"] + Vstart
            Vt.mat[, "Col"] <- Vt.mat[, "Col"] + Vstart
            Vtfull.mat <- rbind(Vtfull.mat, Vt.mat)
        }
        Vstart  <- max(Vtfull.mat[, "Row"])
        Wt <- pevi
        wt <- diag(Wt)
        if(any(wt < 0)) {
            how.many <- sum(wt < 0)
            cat(how.many, " Negative weight(s)\n")
            cat(" Setting the negative weights to a small positive number\n")
            diag(W11)[(1:dim(W11)[1])[diag(W11) < 0]] <- wt[wt < 0] <- min(abs(wt))
        }
        pvals$wt <- wt
        if(n.trait == 1) {
            pvals <- cbind.data.frame(factor(trt[!which.omit]), pvals)
            names(pvals)[1] <- Trait
        }
        pred.lst[[ii]] <- pvals
        V[[ii]] <- Vt
        W[[ii]] <- Wt
        vc[[ii]] <- summ
##        dimnames(W[[ii]]) <- list(nam[!which.omit],nam[!which.omit])
        ##        attr(W[[ii]], "INVERSE") <- TRUE
    }
    Vtfull.mat <- Vtfull.mat[, c(2,1,3)]
    dimnames(Vtfull.mat)[[2]][1:2]  <- dimnames(Vtfull.mat)[[2]][2:1]
    attr(Vtfull.mat, "INVERSE")  <- FALSE
    attr(Vtfull.mat, "rowNames")  <- as.character(1:Vstart)
    bdVt <- bdiag(V)
    dimnames(bdVt) <- list(1:dim(bdVt)[1], 1:dim(bdVt)[1])
    attr(bdVt, "INVERSE") <- FALSE
    bdWt <- bdiag(W)
    dimnames(bdWt) <- list(1:dim(bdWt)[1], 1:dim(bdWt)[1])
    attr(bdWt, "INVERSE") <- TRUE
    pred.df <- do.call("rbind.data.frame", pred.lst)
    names(pred.df)[grep("predicted.value", names(pred.df))] <-
        as.character(asm$call$fixed[[2]])
    pred.df$TG <- factor(1:dim(bdWt)[1])
    names(vc) <- names(mymodels) <- names(V)  <- names(W) <- subnames
    vcMat <- covMat(vc, which.comb, ltrait)
    asreml::asreml.options(asr.options)
    out.list <- list(call = call, Vt = Vtfull.mat, bdVt = bdVt, bdWt = bdWt,
                   V = V, W = W, pred.df = pred.df, type = type, vc = vc,
                   vcMat = vcMat)
    if(keep.models) out.list$models <- mymodels
    invisible(out.list)
}

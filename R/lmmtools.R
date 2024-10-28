# lmmtools documentation
#
#' lmmtools:  A package of tools for linear mixed models
#'
#' The lmmtools package provides a number of functions that
#' supplement the fitting of linear mixed models.  These functions are
#' \enumerate{
#' \item \emph{heritability}: functions for calculating heritability from
#' linear mixed models, including pedigree, genomic relationship matrices
#' and independent effects, and combinations of these.
#' \item \emph{information criteria}: a function to calculate the Akaike and
#' Bayesian Information Criteria from a fit of linear mixed models, in a list
#' of models, when the model are fitted using residual maximum likelihood.
#' \item \emph{initial value function}: for a multivariate autoregressive
#' process of order 1.
#' \item \emph{summary of fa model fits}: a function that provides tests,
#' information criteria, and diagnostics for a list of linear mixed models that
#' have factor analytic models for genetic effects.
#' \item \emph{stage 1 fitting function}: a function that provides for the first
#' stage of fitting in a two stage analysis of multi-site, multi-trait or
#' repeated measures analysis.  Produces the estimated response in a data frame that
#' includes diagonal weights, but also a proper weight matrix for the second stage
#' of analysis.  In addition, a de-regression function is available that implements
#' the original animal breeding approach.
#' }
#' The functions are currently only available as \code{asreml} methods.
#'
#' @docType package
#' @name lmmtools
NULL

#' @title BayesRep package
#' 
#' @description The BayesRep package provides various tools for Bayesian
#'     analysis of replication studies.
#'
#' \code{\link{repPosterior}} visualizes the posterior distribution of the
#'     effect size based on both studies. \code{\link{BFs}} computes the
#'     sceptical Bayes factor (Pawel and Held, 2022), \code{\link{BFr}} computes
#'     the replication Bayes factor (Verhagen and Wagenmakers, 2014), and
#'     \code{\link{BFe}} computes the equality of effect size Bayes factor
#'     (Bayarri and Mayorall, 2002).
#'
#' These functions take effect estimates and their standard errors from original
#'     and replication study as inputs. Throughout, original effect estimate and
#'     standard error are denoted by \code{to} and \code{so} and replication
#'     effect estimate and standard error are denoted \code{tr} and \code{sr}.
#'     It is assumed that each effect estimate is normally distributed around
#'     its true underlying effect size with variance equal to its squared
#'     standard error \deqn{\code{to} \, | \, \theta_o \sim \mathrm{N}(\theta_o,
#'     \code{so}^2) ~ \mathrm{and} ~ \code{tr} \, | \, \theta_r \sim
#'     \mathrm{N}(\theta_r, \code{sr}^2).}{to | theta_o ~ N(theta_o, so^2) and
#'     tr | theta_r ~ N(theta_r, sr^2).} These assumptions may be inadequate for
#'     studies with small sample size (there are special functions for data with
#'     continuous outcomes and standardized mean difference effect size,
#'     \code{\link{BFsSMD}} and \code{\link{BFrSMD}}, and binary outcomes with
#'     log odds ratio effects, \code{\link{BFslogOR}} and
#'     \code{\link{BFrlogOR}}, which are based on the exact distribution of the
#'     data). If not specified otherwise, it is assumed that the true effect
#'     sizes from both studies are the same (\eqn{\theta_o = \theta_r}).
#'
#'
#' @references Bayarri, M. and Mayorall, A. (2002). Bayesian Design of
#'     "Successful" Replications. The American Statistician, 56(3): 207-214.
#'     \doi{10.1198/000313002155}
#'
#' Verhagen, J. and Wagenmakers, E. J. (2014). Bayesian tests to quantify the result
#' of a replication attempt. Journal of Experimental Psychology: General,
#' 145:1457-1475. \doi{10.1037/a0036731}
#'
#' Pawel, S. and Held, L. (2022). The sceptical Bayes factor for the
#'     assessment of replication success. Journal of the Royal Statistical
#'     Society Series B: Statistical Methodology, 84(3): 879-911.
#'     \doi{10.1111/rssb.12491}
#'
#' 
#' @docType package
#' @name BayesRep

NULL

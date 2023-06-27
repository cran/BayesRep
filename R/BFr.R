## non-vectorized version of BFr
.BFr <- function(to, so, tr, sr, ss = 0, truncate = FALSE, log = FALSE,
                 zo = NULL, zr = NULL, c = NULL, g = 0) {
    ## check inputs
    stopifnot(
        length(truncate) == 1,
        is.logical(truncate),
        !is.na(truncate),

        length(log) == 1,
        is.logical(log),
        !is.na(log)
    )
    if (!missing(to) && !missing(so) && !missing(tr) && !missing(sr)) {
        ## parametrization 1
        stopifnot(
            length(to) == 1,
            is.numeric(to),
            is.finite(to),

            length(so) == 1,
            is.numeric(so),
            is.finite(so),
            0 < so,

            length(tr) == 1,
            is.numeric(tr),
            is.finite(tr),

            length(sr) == 1,
            is.numeric(sr),
            is.finite(sr),
            0 < sr,

            length(ss) == 1,
            is.numeric(ss),
            is.finite(ss),
            0 <= ss
        )
        ## compute relative prior variance
        zo <- to/so
        zr <- tr/sr
        c <- so^2/sr^2
        g <- ss^2/so^2
    } else {
        ## parametrization 2
        stopifnot(
            length(zo) == 1,
            is.numeric(zo),
            is.finite(zo),

            length(zr) == 1,
            is.numeric(zr),
            is.finite(zr),

            length(c) == 1,
            is.numeric(c),
            is.finite(c),
            0 < c,

            length(g) == 1,
            is.numeric(g),
            is.finite(g),
            0 <= g
        )
    }

    ## compute bf
    logbf <- stats::dnorm(x = zr, mean = 0, sd = sqrt(1 + c*g), log = TRUE) -
        stats::dnorm(x = zr, mean = zo*sqrt(c), sd = sqrt(1 + c), log = TRUE)
    if (truncate == TRUE) {
        ## add truncation correction factor
        logbf <- logbf + stats::pnorm(q = sign(zo)*zo, log.p = TRUE) -
            stats::pnorm(q = sign(zo)*(zo + zr*sqrt(c))/sqrt(1 + c), log.p = TRUE)
    }
    if (log == TRUE) return(logbf)
    else return(exp(logbf))
}

#' @title Generalized replication Bayes factor
#'
#' @description Computes the generalized replication Bayes factor
#' 
#' @details The generalized replication Bayes factor is the Bayes factor
#'     contrasting the sceptic's hypothesis that the effect size is about zero
#'     \deqn{H_{\mathrm{S}}: \theta \sim \mathrm{N}(0, \code{ss}^2)}{HS: theta ~
#'     N(0, ss^2)} to the advocate's hypothesis that the effect size is
#'     compatible with its posterior distribution based on the original study
#'     and a uniform prior \deqn{H_{\mathrm{A}}: \theta \sim f(\theta \, | \,
#'     \mathrm{original~study}).}{HA: theta ~ f(theta | original study).} The
#'     standard replication Bayes factor from Verhagen and Wagenmakers (2014) is
#'     obtained by specifying a point-null hypothesis \code{ss = 0} (the
#'     default).
#'
#' The function can be used with two input parametrizations, either on the
#' absolute effect scale (\code{to}, \code{so}, \code{tr}, \code{sr}, \code{ss})
#' or alternatively on the relative *z*-scale (\code{zo}, \code{zr}, \code{c},
#' \code{g}). If an argument on the effect scale is missing, the *z*-scale is
#' automatically used and the other non-missing arguments on the effect scale
#' ignored.
#' 
#' @param to Original effect estimate
#' @param so Standard error of the original effect estimate
#' @param tr Replication effect estimate
#' @param sr Standard error of the replication effect estimate
#' @param ss Standard devation of the sceptical prior under
#'     \eqn{H_\mathrm{S}}{HS}. Defaults to \code{0}
#' @param truncate Logical indicating whether advocacy prior should be truncated
#'     to direction of the original effect estimate (i.e., a one-sided test).
#'     Defaults to \code{FALSE}
#' @param log Logical indicating whether the natural logarithm of the Bayes
#'     factor should be returned. Defaults to \code{FALSE}
#' @param zo Original *z*-value \code{zo} = \code{to}/\code{so} (alternative
#'     parametrization for \code{to} and \code{so})
#' @param zr Replication *z*-value \code{zr} = \code{tr}/\code{sr} (alternative
#'     parametrization for \code{tr} and \code{sr})
#' @param c Relative variance \code{c = so^2/sr^2} (alternative parametrization
#'     for \code{so} and \code{sr})
#' @param g Relative prior variance \code{g = ss^2/so^2}. Defaults to \code{0}
#'     (alternative parametrization for \code{ss})
#'
#' @return The generalized replication Bayes factor
#'     \eqn{\mathrm{BF}_{\mathrm{SA}}}{BF_SA}. \eqn{\mathrm{BF}_{\mathrm{SA}} <
#'     1}{BF_SA < 1} indicates that the data favour the advocate's hypothesis
#'     \eqn{H_{\mathrm{A}}}{HA} (replication success), whereas
#'     \eqn{\mathrm{BF}_{\mathrm{SA}} > 1}{BF_SA > 1} indicates that the data
#'     favour the sceptic's hypothesis \eqn{H_{\mathrm{S}}}{HS} (replication
#'     failure).
#' 
#' @seealso \code{\link{BFrSMD}}, \code{\link{BFrlogOR}}
#' 
#' @author Samuel Pawel
#'
#' @references Verhagen, J. and Wagenmakers, E. J. (2014). Bayesian tests to
#'     quantify the result of a replication attempt. Journal of Experimental
#'     Psychology: General, 145:1457-1475. \doi{10.1037/a0036731}
#'
#' Ly, A., Etz, A., Marsman, M., Wagenmakers, E. J. (2019). Replication Bayes
#' factors from evidence updating. Behavior Research Methods, 51(6):2498-2508.
#' \doi{10.3758/s13428-018-1092-x}
#'
#' Pawel, S. and Held, L. (2022). The sceptical Bayes factor for the
#'     assessment of replication success. Journal of the Royal Statistical
#'     Society Series B: Statistical Methodology, 84(3): 879-911.
#'     \doi{10.1111/rssb.12491}
#' 
#' @examples
#' to <- 2
#' tr <- 2.5
#' so <- 1
#' sr <- 1
#' BFr(to = to, so = so, tr = tr, sr = sr)
#' BFr(zo = to/so, zr = tr/sr, c = so^2/sr^2)
#'
#'  
#' @export 
BFr <- Vectorize(.BFr)


## non-vectorized version of BFrlogOR
.BFrlogOR <- function(ao, bo, nTo = ao + bo, co, do, nCo = co + do,
                      ar, br, nTr = ar + br, cr, dr, nCr = cr + dr, ss,
                      method = c("integration", "hypergeo")) {

    ## check inputs
    ## check inputs
    if (!missing(bo)) {
        ## parametrization 1
        stopifnot(
            length(bo) == 1,
            is.numeric(bo),
            is.finite(bo),
            0 <= bo
        )
    }
    if (!missing(do)) {
        ## parametrization 1
        stopifnot(
            length(do) == 1,
            is.numeric(do),
            is.finite(do),
            0 <= do
        )
    }
    if (!missing(br)) {
        ## parametrization 1
        stopifnot(
            length(br) == 1,
            is.numeric(br),
            is.finite(br),
            0 <= br
            )
    }
    if (!missing(dr)) {
        ## parametrization 1
        stopifnot(
            length(dr) == 1,
            is.numeric(dr),
            is.finite(dr),
            0 <= dr
            )
    }
    stopifnot(
        length(ao) == 1,
        is.numeric(ao),
        is.finite(ao),
        0 <= ao,

        length(nTo) == 1,
        is.numeric(nTo),
        is.finite(nTo),
        0 <= nTo,

        length(co) == 1,
        is.numeric(co),
        is.finite(co),
        0 <= co,

        length(nCo) == 1,
        is.numeric(nCo),
        is.finite(nCo),
        0 <= nCo,

        length(ar) == 1,
        is.numeric(ar),
        is.finite(ar),
        0 <= ar,

        length(nTr) == 1,
        is.numeric(nTr),
        is.finite(nTr),
        0 <= nTr,

        length(cr) == 1,
        is.numeric(cr),
        is.finite(cr),
        0 <= cr,

        length(nCr) == 1,
        is.numeric(nCr),
        is.finite(nCr),
        0 <= nCr,

        length(ss) == 1,
        is.numeric(ss),
        is.finite(ss),
        0 <= ss,

        !is.null(method)
    )
    method <- match.arg(method)

    ## exact likelihood function
    likExact <- function(a, nT, c, nC, logOR) {
        ## integrate out the proportion in the control group from the likelihood
        ## using the translation-invariant Jeffreys' prior (p0*(1 - po))^(-0.5)
        intFun <- function(p0) {
            stats::dbinom(x = a, size = nT,
                          prob = 1/(1 + exp(-(log(p0/(1 - p0)) + logOR)))) *
                stats::dbinom(x = c, size = nC, prob = p0) *
                stats::dbeta(x = p0, shape1 = 0.5, shape2 = 0.5)
        }
        lik <- try(stats::integrate(f = intFun, lower = 0, upper = 1)$value,
                   silent = TRUE)
        if (inherits(lik, "try-error")) return(NaN)
        else return(lik)
    }

    ## compute marginal likelihood under sceptical prior H_S
    if (ss == 0) {
        ## point-null
        mS <- likExact(a = ar, nT = nTr, c = cr, nC = nCr, logOR = 0)
    } else {
        ## composite-null
        mS <- try(stats::integrate(f = function(logOR) {
            vapply(X = logOR, FUN = function(logOR) {
                likExact(a = ar, nT = nTr, c = cr, nC = nCr, logOR = logOR) *
                    stats::dnorm(x = logOR, mean = 0, sd = ss)
            }, FUN.VALUE = 1)
        }, lower = -Inf, upper = Inf)$value, silent = TRUE)
        if (inherits(mS, "try-error")) return(NaN)
    }

    ## compute marginal likelihood under advocacy prior H_A:
    ## 1) posterior density of logOR based on original study and Jeffreys priors
    if (method == "hypergeo") {
        ## 1a) analytical solution using hypergeometric function
        e <- ao + 0.5
        f <- nTo - ao + 0.5
        g <- co + 0.5
        h <- nCo - co + 0.5
        C <- beta(a = e + g, b = f + h)/(beta(a = e, b = f) * beta(a = g, b = h))
        postExact <- function(logOR) {
            suppressWarnings({
                postdens <-  ifelse(logOR < 0,
                                    hypergeo::hypergeo(A = e + f, B = e + g,
                                                       C = e + f + g + h,
                                                       z = 1 - exp(logOR))*C*exp(e*logOR),
                                    hypergeo::hypergeo(A = e + f, B = f + h,
                                                       C = e + f + g + h,
                                                       z = 1 - exp(-logOR))*C*exp(-f*logOR))
            })
            return(abs(postdens))
        }
    } else {
        ## 1b) using analytical posterior for p0 and p1, then CoV and integrating out p0
        postExact <- function(logOR) {
            intFun <- function(p0) {
                stats::dbeta(x = stats::plogis(q = logOR + log(p0/(1 - p0))),
                             shape1 = 0.5 + ao, shape2 = 0.5 + nTo - ao) *
                    exp(logOR) * (1 - p0) * p0 / (p0 * exp(logOR) + (1 - p0))^2 * ## CoV constant
                    stats::dbeta(x = p0, shape1 = 0.5 + co, shape2 = 0.5 + nCo - co)
            }
            postdens <- try(stats::integrate(f = intFun, lower = 0, upper = 1)$value,
                            silent = TRUE)
            if (inherits(postdens, "try-error")) return(NaN)
            else return(postdens)
        }
    }
    ## 2) integrate replication likelihood wrt to posterior of logOR
    mA <- try(stats::integrate(f = function(logOR) {
        vapply(X = logOR, FUN = function(logOR) {
            likExact(a = ar, nT = nTr, c = cr, nC = nCr, logOR = logOR) *
                postExact(logOR = logOR)
        }, FUN.VALUE = 1)
    }, lower = -Inf, upper = Inf)$value, silent = TRUE)
    if (inherits(mA, "try-error")) return(NaN)

    ## compute BF
    bfSA <- mS/mA
    return(bfSA)
}

#' @title Generalized replication Bayes factor for logOR effect sizes
#'
#' @description Computes the generalized replication Bayes factor for log odds
#'     ratio (logOR) effect sizes
#'
#' @details This function computes the generalized replication Bayes factor for
#'     log odds ratio (logOR) effect sizes using an exact binomial likelihood
#'     for the data instead of the normal approximation used in
#'     \code{\link{BFr}} (for details, see Section 4 in Pawel and Held, 2022).
#'
#'
#' @param ao Number of cases in original study treatment group
#' @param bo Number of non-cases in original study treatment group
#' @param nTo Number of participants in original study treatment group (specify
#'     alternatively to \code{b})
#' @param co Number of cases in original study control group
#' @param do Number of non-cases in original study control group
#' @param nCo Number of participants in original study control group (specify
#'     alternatively to \code{d})
#' @param ar Number of cases in replication study treatment group
#' @param br Number of non-cases in replication study treatment group
#' @param nTr Number of participants in replication study treatment group
#'     (specify alternatively to \code{b})
#' @param cr Number of cases in replication study control group
#' @param dr Number of non-cases in replication study control group
#' @param nCr Number of participants in replication study control group (specify
#'     alternatively to \code{d})
#' @param ss Standard deviation of the sceptical prior under
#'     \eqn{H_\mathrm{S}}{HS}. Defaults to \code{0}
#' @param method Method to compute posterior density. Either
#'     \code{"integration"} (default) or \code{"hypergeo"}
#'
#' @return The generalized replication Bayes factor
#'     \eqn{\mathrm{BF}_{\mathrm{SA}}}{BF_SA}. \eqn{\mathrm{BF}_{\mathrm{SA}} <
#'     1}{BF_SA < 1} indicates that the data favour the advocate's hypothesis
#'     \eqn{H_{\mathrm{A}}}{HA} (replication success), whereas
#'     \eqn{\mathrm{BF}_{\mathrm{SA}} > 1}{BF_SA > 1} indicates that the data
#'     favour the sceptic's hypothesis \eqn{H_{\mathrm{S}}}{HS} (replication
#'     failure).
#'
#' @references Verhagen, J. and Wagenmakers, E. J. (2014). Bayesian tests to
#'     quantify the result of a replication attempt. Journal of Experimental
#'     Psychology: General, 145:1457-1475. \doi{10.1037/a0036731}
#'
#' Pawel, S. and Held, L. (2022). The sceptical Bayes factor for the assessment
#'     of replication success. Journal of the Royal Statistical Society Series
#'     B: Statistical Methodology, 84(3): 879-911. \doi{10.1111/rssb.12491}
#'
#' @author Samuel Pawel
#'
#' @examples
#' data("SSRPexact")
#' balafoutas2012 <- subset(SSRPexact, study == "Balafoutas and Sutter (2012), Science")
#' with(balafoutas2012,
#'      BFrlogOR(ao = ao, bo = bo, co = co, do = do, ar = ar, br = br, cr = cr, dr = dr,
#'               ss = 0))
#'
#' @export
BFrlogOR <- Vectorize(.BFrlogOR)

## non-vectorized version of BFrSMD
.BFrSMD <- function(to, no, n1o = no, n2o = no, tr, nr, n1r = nr, n2r = nr,
                    ss, type = c("two.sample", "one.sample", "paired")) {

    ## input checks
    stopifnot(
        length(to) == 1,
        is.numeric(to),
        is.finite(to),

        length(tr) == 1,
        is.numeric(tr),
        is.finite(tr),

        length(n1o) == 1,
        is.numeric(n1o),
        is.finite(n1o),
        0 < n1o,

        length(n2o) == 1,
        is.numeric(n2o),
        is.finite(n2o),
        0 < n2o,

        length(n1r) == 1,
        is.numeric(n1r),
        is.finite(n1r),
        0 < n1r,

        length(n2r) == 1,
        is.numeric(n2r),
        is.finite(n2r),
        0 < n2r,

        length(ss) == 1,
        is.numeric(ss),
        is.finite(ss),
        0 <= ss,

        !is.null(type)
    )
    type <- match.arg(type)
    if (type != "two.sample") {
        if (n1o != n2o) {
            warning(paste0('different n1o and n2o supplied but type set to "', type,
                           '", using no = n1o'))
        }
        if (n1r != n2r) {
            warning(paste0('different n1r and n2r supplied but type set to "', type,
                           '", using nr = n1r'))
        }
    }

    ## compute df and effective sample size depending on test type
    if (type == "two.sample") {
        df <- n1r + n2r - 2
        nstar <- 1/(1/n1r + 1/n2r)
        dfo <- n1o + n2o - 2
        nstaro <- 1/(1/n1o + 1/n2o)
    } else {
        df <- n1r - 1
        nstar <- n1r
        dfo <- n1o -1
        nstaro <- n1o
    }

    ## compute marginal likelihood under sceptical prior H_S
    if (ss == 0) {
        ## point-null
        mS <- stats::dt(x = tr, df = df, ncp = 0)
    } else {
        ## composite-null
        mS <- try(stats::integrate(f = function(SMD) {
            suppressWarnings({
                stats::dt(x = tr, df = df, ncp = SMD*sqrt(nstar)) *
                    stats::dnorm(x = SMD, mean = 0, sd = ss)
            })
        }, lower = -Inf, upper = Inf)$value, silent = TRUE)
        if (inherits(mS, "try-error")) return(NaN)
    }

    ## compute marginal likelihood under advocacy prior H_A:
    ## 1) posterior density of SMD based on original study and flat prior
    postExact <- function(SMD) {
        postdens <- try(stats::integrate(f = function(prec) {
            ## stats::dnorm(x = SMD, mean = to*prec/sqrt(nstaro), sd = sqrt(1/nstaro)) *
            ##     stats::dgamma(x = prec, shape = dfo + 1, rate = dfo)
            stats::dnorm(x = SMD, mean = to*sqrt(prec/nstaro), sd = sqrt(1/nstaro)) *
                stats::dgamma(x = prec, shape = (dfo + 1)/2, rate = dfo/2)
        }, lower = 0, upper = Inf)$value, silent = TRUE)
        if (inherits(postdens, "try-error")) return(NaN)
        else return(postdens)
    }
    ## 2) integrate tr likelihood wrt to posterior of SMD
    mA <- try(stats::integrate(f = function(SMD) {
        vapply(X = SMD, FUN = function(SMD) {
            suppressWarnings({
                stats::dt(x = tr, df = df, ncp = SMD*sqrt(nstar)) * postExact(SMD)
            })
        }, FUN.VALUE = 1)
    }, lower = -Inf, upper = Inf)$value, silent = TRUE)
    if (inherits(mA, "try-error")) return(NaN)

    ## compute BF
    bfSA <- mS/mA
    return(bfSA)
}


#' @title Generalized replication Bayes factor for SMD effect sizes
#'
#' @description Computes the generalized replication Bayes factor for
#'     standardized mean difference (SMD) effect sizes
#'
#' @details This function computes the generalized replication Bayes factor for
#'     standardized mean difference (SMD) effect sizes using an exact
#'     *t*-likelihood for the data instead of the normal approximation used in
#'     \code{\link{BFr}} (for details, see Section 4 in Pawel and Held, 2022).
#'     Data from both studies are summarized by \eqn{t}-statistics and sample
#'     sizes. The following types of \eqn{t}-tests are accepted:
#'
#' - Two-sample \eqn{t}-test where the SMD represents the standardized
#' mean difference between two group means (assuming equal variances in
#' both groups).
#' - One-sample \eqn{t}-test where the SMD represents the standardized
#' mean difference to the null value.
#' - Paired \eqn{t}-test where the SMD represents the standardized mean
#' difference score.
#'
#' @param to \eqn{t}-statistic from the original study
#' @param no Sample size of the original study (per group)
#' @param n1o Sample size in group 1 of the original study (only required for
#'     two-sample \eqn{t}-test with unequal group sizes)
#' @param n2o Sample size in group 2 of the original study (only specify if
#'     unequal group sizes)
#' @param tr \eqn{t}-statistic from the replication study
#' @param nr Sample size of the replication study (per group)
#' @param n1r Sample size in group 1 of the replication study (only required for
#'     two-sample \eqn{t}-test with unequal group sizes)
#' @param n2r Sample size in group 2 of the replication study (only required for
#'     two-sample \eqn{t}-test with unequal group sizes)
#' @param ss Standard devation of the sceptical prior under
#'     \eqn{H_\mathrm{S}}{HS}. Defaults to \code{0}
#' @param type Type of \eqn{t}-test associated with \eqn{t}-statistic. Can be
#'     `"two.sample"`, `"one.sample"`, `"paired"`. Defaults to `"two.sample"`
#'
#' @return The generalized replication Bayes factor
#'     \eqn{\mathrm{BF}_{\mathrm{SA}}}{BF_SA}. \eqn{\mathrm{BF}_{\mathrm{SA}} <
#'     1}{BF_SA < 1} indicates that the data favour the advocate's hypothesis
#'     \eqn{H_{\mathrm{A}}}{HA} (replication success), whereas
#'     \eqn{\mathrm{BF}_{\mathrm{SA}} > 1}{BF_SA > 1} indicates that the data
#'     favour the sceptic's hypothesis \eqn{H_{\mathrm{S}}}{HS} (replication
#'     failure).
#'
#' @seealso \code{\link{BFr}}, \code{\link{BFrlogOR}}
#'
#' @references Verhagen, J. and Wagenmakers, E. J. (2014). Bayesian tests to
#'     quantify the result of a replication attempt. Journal of Experimental
#'     Psychology: General, 145:1457-1475. \doi{10.1037/a0036731}
#'
#' Pawel, S. and Held, L. (2022). The sceptical Bayes factor for the assessment
#'     of replication success. Journal of the Royal Statistical Society Series
#'     B: Statistical Methodology, 84(3): 879-911. \doi{10.1111/rssb.12491}
#'
#' @author Samuel Pawel
#'
#' @examples
#' data("SSRPexact")
#' morewedge2010 <- subset(SSRPexact, study == "Morewedge et al. (2010), Science")
#' with(morewedge2010,
#'      BFrSMD(to = to, n1o = n1o, n2o = n2o, tr = tr, n1r = n1r, n2r = n2r, ss = 0))
#'
#' @export
BFrSMD <- Vectorize(.BFrSMD)

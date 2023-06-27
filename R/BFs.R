## non-vectorized version of BFs
.BFs <- function(to, so, tr, sr, truncate = FALSE, zo = NULL, zr = NULL,
                 c = NULL) {
    ## check inputs
    stopifnot(
        length(truncate) == 1,
        is.logical(truncate),
        !is.na(truncate)
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
            0 < sr
        )
        ## compute relative quantities
        zo <- to/so
        zr <- tr/sr
        c <- so^2/sr^2
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
            0 < c
        )
    }

    ## compute some quantities for initial checks
    gMinBFo <- max(c(zo^2 - 1, 0)) ## relative prior variance at minBFo
    minBFo <- BFo(zo = zo, g = gMinBFo) ## minBFo
    BFr0 <- BFr(zo = zo, zr = zr, c = c, g = 0,
                truncate = truncate) ## Replication BF (BF_SA(zr; g=0))
    BFrgMinBFo <- BFr(zo = zo, zr = zr, c = c, g = gMinBFo,
                      truncate = truncate) ## BF_SA(zr; g=gminBFo)

    ## when BF_SA(zr; gMinBFo) <= minBFo:
    ## BFs = minBFo
    if (BFrgMinBFo <= minBFo) {
        vss <- gMinBFo

        ## otherwise determine g* where BFo and BFr intersect:
        ## BFs = BFo(g*) = BFr(g*)
    } else {

        ## when c = 1: use Lambert W function to compute g*
        if (c == 1) {
            x <- 0.5*(zo^2 + zr^2)/sqrt(2)*exp(-0.5*(zo^2 + 0.5*(zr - zo)^2))
            if (truncate == TRUE) {
                x <- x*stats::pnorm(q = sign(zo)*(zo + zr)/sqrt(2))/stats::pnorm(q = abs(zo))
            }
            res <- -0.5*(zo^2 + zr^2)/lamW::lambertWm1(x = -x) - 1
            if (res < 0) {
                vss <- NaN
            } else {
                vss <- res
            }

            ## when c != 1: compute g* numerically
        } else {
            BFdiff <- function(g) {
                logBFo <- BFo(zo = zo, g = g, log = TRUE)
                logBFr <- BFr(zo = zo, zr = zr, c = c, g = g, log = TRUE,
                              truncate = truncate)
                return(logBFr - logBFo)
            }
            res <- try(stats::uniroot(f = BFdiff, lower = 0, upper = gMinBFo)$root,
                       silent = TRUE)
            if (inherits(res, "try-error")) {
                vss <- NaN
            } else {
                vss <- res
            }
        }
    }

    ## return BFs (and vss if specified)
    if (is.nan(vss)) {
        bfs <- NaN
    } else {
        bfs <- BFo(zo = zo, g = vss)
    }
    return(bfs)

}

#' @title Sceptical Bayes factor
#'
#' @description Computes the sceptical Bayes factor
#'
#' @details The sceptical Bayes factor is a summary measure of the following
#'     two-step reverse-Bayes procedure for assessing replication success:
#'
#' 1. Use the data from the original study to determine the standard deviation
#' \eqn{\tau_{\gamma}}{tau_gamma} of a sceptical normal prior \eqn{\theta \sim
#' \mathrm{N}(0, \tau_{\gamma}^2)}{HS: theta ~ N(0, tau_gamma^2)} such that the
#' Bayes factor contrasting the null hypothesis \eqn{H_0: \theta = 0}{H0: theta
#' = 0} to the sceptic's hypothesis \eqn{H_{\mathrm{S}}: \theta \sim
#' \mathrm{N}(0, \tau_{\gamma}^2)}{HS: theta ~ N(0, tau_gamma^2)} equals a
#' specified level \eqn{\gamma \in (0, 1]}{gamma in (0 1]}. This prior
#' represents a sceptic who remains unconvinced about the presence of an effect
#' at level \eqn{\gamma}{gamma}.
#'
#' 2. Use the data from the replication study to compare the sceptic's
#' hypothesis \eqn{H_{\mathrm{S}}: \theta \sim \mathrm{N}(0,
#' \tau_{\gamma}^2)}{HS: theta ~ N(0, tau_gamma^2)} to the advocate's hypothesis
#' \eqn{H_{\mathrm{A}}: \theta \sim f(\theta \, | \,
#' \mathrm{original~study})}{HA: theta ~ f(theta | original study)}. The prior of
#' the effect size under \eqn{H_{\mathrm{A}}}{HA} is its posterior based on the
#' original study and a uniform prior, thereby representing the position of an
#' advocate of the original study. Replication success at level
#' \eqn{\gamma}{gamma} is achieved if the Bayes factor contrasting
#' \eqn{H_{\mathrm{S}}}{HS} to \eqn{H_{\mathrm{A}}}{HA} is smaller than
#' \eqn{\gamma}{gamma}, which means that the replication data favour the
#' advocate over the sceptic at a higher level than the sceptic's initial
#' objection. The sceptical Bayes factor \eqn{\mathrm{BF}_{\mathrm{S}}}{BF_S} is
#' the smallest level \eqn{\gamma}{gamma} at which replication success can be
#' established.
#'
#' The function can be used with two input parametrizations, either on the
#' absolute effect scale (\code{to}, \code{so}, \code{tr}, \code{sr}) or
#' alternatively on the relative *z*-scale (\code{zo}, \code{zr}, \code{c}). If
#' an argument on the effect scale is missing, the *z*-scale is automatically
#' used and the other non-missing arguments on the effect scale ignored.
#'
#' @param to Original effect estimate
#' @param so Standard error of the original effect estimate
#' @param tr Replication effect estimate
#' @param sr Standard error of the replication effect estimate
#' @param truncate Logical indicating whether advocacy prior should be truncated
#'     to direction of the original effect estimate (i.e., a one-sided test).
#'     Defaults to \code{FALSE}
#' @param zo Original *z*-value \code{zo} = \code{to}/\code{so} (alternative
#'     parametrization for \code{to} and \code{so})
#' @param zr Replication *z*-value \code{zr} = \code{tr}/\code{sr} (alternative
#'     parametrization for \code{tr} and \code{sr})
#' @param c Relative variance \code{c = so^2/sr^2} (alternative parametrization
#'     for \code{so} and \code{sr})
#'
#' @return The sceptical Bayes factor \eqn{\mathrm{BF}_{\mathrm{S}}}{BF_S}.
#'     \eqn{\mathrm{BF}_{\mathrm{S}} < 1}{BF_S < 1} indicates replication
#'     success, the smaller the value of \eqn{\mathrm{BF}_{\mathrm{S}}}{BF_S}
#'     the higher the degree of replication success. It is possible that the
#'     result of the replication is so inconclusive that replication success
#'     cannot be established at any level. In this case, the sceptical Bayes
#'     factor does not exist and the function returns \code{NaN}.
#'
#'
#' @references Pawel, S. and Held, L. (2022). The sceptical Bayes factor for the
#'     assessment of replication success. Journal of the Royal Statistical
#'     Society Series B: Statistical Methodology, 84(3): 879-911.
#'     \doi{10.1111/rssb.12491}
#'
#' @author Samuel Pawel
#'
#' @seealso \code{\link{BFsSMD}}, \code{\link{BFslogOR}}
#'
#' @examples
#' to <- 2
#' tr <- 2.5
#' so <- 1
#' sr <- 1
#' BFs(to = to, so = so, tr = tr, sr = sr)
#' BFs(zo = to/so, zr = tr/sr, c = so^2/sr^2)
#'
#' @export
BFs <- Vectorize(.BFs)

## non-vectorized version of BFslogOR
.BFslogOR <- function(ao, bo, nTo = ao + bo, co, do, nCo = co + do, ar, br,
                      nTr = ar + br, cr, dr, nCr = cr + dr,
                      method = c("integration", "hypergeo")) {

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

        !is.null(method)
    )
    method <- match.arg(method)

    ## compute minimum BF for original data
    minBFoptim <- stats::optim(par = 0, fn = function(ss2) {
        ## for some reasons this optimization works better on the variance scale
        log(BFologOR(ao = ao, nTo = nTo, co = co, nCo = nCo, ss = sqrt(ss2)))
    }, method = "L-BFGS-B", lower = 0, upper = Inf)
    ssminBFo <- sqrt(minBFoptim$par)
    minBFo <- exp(minBFoptim$value)

    ## when ssmin = 0, check by hand whether success at level = 1
    if (ssminBFo == 0) {
        BFr0 <- BFrlogOR(ao = ao, nTo = nTo, co = co, nCo = nCo, ar = ar,
                         nTr = nTr, cr = cr, nCr = nCr, ss = 0, method = method)
        if (is.nan(BFr0)) {
            warnMessage <- paste("numerical problems when computing posterior density, try method =",
                                 ifelse(method == "integration", "'hypergeo'",
                                        "'integration'"))
            warning(warnMessage)
            ssSceptical <- NaN
            bfs <- NaN
        } else if (BFr0 <= 1) {
            ssSceptical <- 0
            bfs <- 1
        } else {
            ssSceptical <- NaN
            bfs <- NaN
        }
    } else {
        ## check by hand whether success at level minBFo
        BFrminBFo <- BFrlogOR(ao = ao, nTo = nTo, co = co, nCo = nCo, ar = ar,
                              nTr = nTr, cr = cr, nCr = nCr, ss = ssminBFo,
                              method = method)
        if (is.nan(BFrminBFo)) {
            warnMessage <- paste("numerical problems when computing posterior density, try method =",
                                 ifelse(method == "integration", "'hypergeo'",
                                        "'integration'"))
            warning(warnMessage)
            ssSceptical <- NaN
            bfs <- NaN
        } else if (BFrminBFo <= minBFo) {
            ssSceptical <- ssminBFo
            bfs <- minBFo
        } else {
            ## otherwise use uniroot
            rootFun <- function(ss2) {
                res <- BFrlogOR(ao = ao, nTo = nTo, co = co, nCo = nCo, ar = ar,
                                nTr = nTr, cr = cr, nCr = nCr, ss = sqrt(ss2),
                                method = method) -
                    BFologOR(ao = ao, nTo = nTo, co = co, nCo = nCo, ss = sqrt(ss2))
                return(res)
            }
            rootRes <- try(stats::uniroot(f = rootFun, interval = c(0, ssminBFo^2)),
                           silent = TRUE)
            if (inherits(rootRes, "try-error")) {
                ssSceptical <- NaN
                bfs <- NaN
            } else {
                ssSceptical <- sqrt(rootRes$root)
                bfs <- BFologOR(ao = ao, nTo = nTo, co = co, nCo = nCo, ss = ssSceptical)
            }
        }
    }
    return(bfs)
}


#' @title Sceptical Bayes factor for logOR effect sizes
#'
#' @description Computes the sceptical Bayes factor for logOR effect sizes
#'
#' @details This function computes the sceptical Bayes factor for log odds ratio
#'     (logOR) effect sizes using an exact binomial likelihood for the data
#'     instead of the normal approximation used in \code{\link{BFs}} (for
#'     details, see Section 4 in Pawel and Held, 2022).
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
#' @param method Method to compute posterior density. Either
#'     \code{"integration"} (default) or \code{"hypergeo"}
#'
#' @return The sceptical Bayes factor \eqn{\mathrm{BF}_{\mathrm{S}}}{BF_S}.
#'     \eqn{\mathrm{BF}_{\mathrm{S}} < 1}{BF_S < 1} indicates replication
#'     success, the smaller the value of \eqn{\mathrm{BF}_{\mathrm{S}}}{BF_S}
#'     the higher the degree of replication success. It is possible that the
#'     result of the replication is so inconclusive that replication success
#'     cannot be established at any level. In this case, the sceptical Bayes
#'     factor does not exist and the function returns \code{NaN}.
#'
#' @author Samuel Pawel
#'
#' @references Pawel, S. and Held, L. (2022). The sceptical Bayes factor for the
#'     assessment of replication success. Journal of the Royal Statistical
#'     Society Series B: Statistical Methodology, 84(3): 879-911.
#'     \doi{10.1111/rssb.12491}
#'
#' @seealso \code{\link{BFs}}, \code{\link{BFslogOR}}
#'
#' @examples
#' data("SSRPexact")
#' balafoutas2012 <- subset(SSRPexact, study == "Balafoutas and Sutter (2012), Science")
#' with(balafoutas2012,
#'      BFslogOR(ao = ao, bo = bo, co = co, do = do, ar = ar, br = br, cr = cr, dr = dr))
#'
#' @export
BFslogOR <- Vectorize(.BFslogOR)


## non-vectorized version of BFsSMD
.BFsSMD <- function(to, no, n1o = no, n2o = no, tr, nr, n1r = nr, n2r = nr,
                     type = c("two.sample", "one.sample", "paired")) {

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

        !is.null(type)
    )
    type <- match.arg(type)
    if (type != "two.sample") {
        startpar <- max(c(to^2 - 1, 0))/(1/n1o + 1/n2o)
        if (n1o != n2o) {
            warning(paste0('different n1o and n2o supplied but type set to "', type,
                           '", using no = n1o'))
        }
        if (n1r != n2r) {
            warning(paste0('different n1r and n2r supplied but type set to "', type,
                           '", using nr = n1r'))
        }
    } else {
        startpar <- max(c(to^2 - 1, 0))/n1o
    }

    suppressWarnings({
        ## compute minimum BF for original data
        minBFoptim <- stats::optim(par = startpar, fn = function(ss) {
            BFoSMD(to = to, n1o = n1o, n2o = n2o, ss = ss, type = type)
        }, method = "L-BFGS-B", lower = 0, upper = Inf)
        ssminBFo <- minBFoptim$par
        minBFo <- minBFoptim$value

        ## when ssmin = 0, check by hand whether success at level = 1
        if (ssminBFo == 0) {
            BFr0 <- BFrSMD(to = to, n1o = n1o, n2o = n2o, tr = tr, n1r = n1r, n2r = n2r,
                           ss = 0, type = type)
            if (BFr0 <= 1) {
                ssSceptical <- 0
                bfs <- 1
            } else {
                ssSceptical <- NaN
                bfs <- NaN
            }
        } else {
            ## check by hand whether success at level = minBFo
            BFrminBFo <- BFrSMD(to = to, n1o = n1o, n2o = n2o, tr = tr, n1r = n1r, n2r = n2r,
                                ss = ssminBFo, type = type)
            if (BFrminBFo <= minBFo) {
                ssSceptical <- ssminBFo
                bfs <- minBFo
            } else {
                ## otherwise use uniroot
                rootFun <- function(ss2) {
                    res <- BFrSMD(to = to, n1o = n1o, n2o = n2o, tr = tr, n1r = n1r, n2r = n2r,
                                  ss = sqrt(ss2), type = type) -
                        BFoSMD(to = to, n1o = n1o, n2o = n2o, ss = sqrt(ss2), type = type)
                    return(res)
                }
                rootRes <- try(stats::uniroot(f = rootFun, interval = c(0, ssminBFo^2)),
                               silent = TRUE)
                if (inherits(rootRes, "try-error")) {
                    ssSceptical <- NaN
                    bfs <- NaN
                } else {
                    ssSceptical <- sqrt(rootRes$root)
                    bfs <- BFoSMD(to = to, n1o = n1o, n2o = n2o, ss = ssSceptical, type = type)
                }
            }
        }

        return(bfs)
    })
}


#' @title Sceptical Bayes factor for SMD effect sizes
#'
#' @description Computes the sceptical Bayes factor for standardized mean
#'     difference (SMD) effect sizes
#'
#' @details This function computes the sceptical Bayes factor for standardized
#'     mean difference (SMD) effect sizes using an exact *t*-likelihood for the
#'     data instead of the normal approximation used in \code{\link{BFs}} (for
#'     details, see Section 4 in Pawel and Held, 2022). Data from both studies
#'     are summarized by \eqn{t}-statistics and sample sizes. The following
#'     types of \eqn{t}-tests are accepted:
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
#' @param type Type of \eqn{t}-test associated with \eqn{t}-statistic. Can be
#'     `"two.sample"`, `"one.sample"`, `"paired"`. Defaults to `"two.sample"`.
#'
#' @return The sceptical Bayes factor \eqn{\mathrm{BF}_{\mathrm{S}}}{BF_S}.
#'     \eqn{\mathrm{BF}_{\mathrm{S}} < 1}{BF_S < 1} indicates replication
#'     success, the smaller the value of \eqn{\mathrm{BF}_{\mathrm{S}}}{BF_S}
#'     the higher the degree of replication success. It is possible that the
#'     result of the replication is so inconclusive that replication success
#'     cannot be established at any level. In this case, the sceptical Bayes
#'     factor does not exist and the function returns \code{NaN}.
#'
#' @author Samuel Pawel
#'
#' @references Pawel, S. and Held, L. (2022). The sceptical Bayes factor for the
#'     assessment of replication success. Journal of the Royal Statistical
#'     Society Series B: Statistical Methodology, 84(3): 879-911.
#'     \doi{10.1111/rssb.12491}
#'
#' @seealso \code{\link{BFs}}, \code{\link{BFslogOR}}
#'
#' @examples
#' data("SSRPexact")
#' morewedge2010 <- subset(SSRPexact, study == "Morewedge et al. (2010), Science")
#' with(morewedge2010,
#'      BFsSMD(to = to, n1o = n1o, n2o = n2o, tr = tr, n1r = n1r, n2r = n2r))
#'
#' @export
#'
BFsSMD <- Vectorize(.BFsSMD)

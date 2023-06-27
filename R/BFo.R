## non-vectorized version of BFo
.BFo <- function(zo, g, log = FALSE, truncate = FALSE) {
    ## check inputs
    stopifnot(
        length(zo) == 1,
        is.numeric(zo),
        is.finite(zo),

        length(g) == 1,
        is.numeric(g),
        is.finite(g),
        0 <= g,

        length(log) == 1,
        is.logical(log),

        length(truncate) == 1,
        is.logical(truncate)
    )

    ## compute bf
    logbf <- stats::dnorm(x = zo, mean = 0, sd = 1, log = TRUE) -
        stats::dnorm(x = zo, mean = 0, sd = sqrt(1 + g), log = TRUE)

    ## add correction factor when truncated alternative
    if (truncate == TRUE) {
        logbf <- logbf - log(2) -
            stats::pnorm(q = zo*sqrt(g/(1 + g)), log.p = TRUE)
    }

    if (log == TRUE) return(logbf)
    else return(exp(logbf))
}

#' @title Bayes factor for effect estimate from original study
#' 
#' @description Computes Bayes factor contrasting \deqn{H_0: \theta = 0}{H0:
#'     theta = 0} to \deqn{H_S: \theta \sim \mathrm{N}(0, g \cdot
#'     \sigma_o^2)}{HS: theta ~ N(0, g*sigma_o^2)} with respect to the effect
#'     estimate from an original study \eqn{\hat{\theta}_o}{hat(theta)_o},
#'     assumed to be approximately normally distributed, i.e.
#'     \eqn{\hat{\theta}_o \, | \, \theta \sim \mathrm{N}(\theta, \sigma_o^2),}{
#'     hat(theta)_o | theta ~ N(theta, sigma_o^2),} with known standard error
#'     \eqn{\sigma_o}{sigma_o}.
#' 
#' @param zo \eqn{z}{z}-value from original study, \eqn{z_o =
#'     \hat{\theta}_o/\sigma_o}{zo = hat(theta)_o/sigma_o}, i.e. original effect
#'     estimate divided by standard error
#' @param g Relative variance of \eqn{\mathrm{N}(0, g \cdot \sigma_o^2)}{N(0,
#'     g*sigma_o^2)} prior for the effect size under \eqn{H_S}{HS}
#' @param log Logical indicating whether the natural logarithm of the Bayes
#'     factor should be returned)
#' @param truncate Logical indicating whether a truncated alternative should be
#'     used (truncated in the direction of the data).
#' 
#' @return Bayes factor (BF \eqn{< 1} indicates that data favour \eqn{H_S}{HS},
#'     while BF \eqn{> 1} indicates that data favour \eqn{H_0}{H0})
#' 
#' @author Samuel Pawel
#'
#' @references Pawel, S. and Held, L. (2022). The sceptical Bayes factor for the
#'     assessment of replication success. Journal of the Royal Statistical
#'     Society Series B: Statistical Methodology. 84(3): 879-911.
#'     \doi{10.1111/rssb.12491}
#'
#' @seealso \code{\link{BFoSMD}}, \code{\link{BFologOR}}
#' 
#' @examples BFo(zo = 3, g = 1)
#'
#' @noRd
BFo <- Vectorize(.BFo)

## ## truncation test
## margf <- function(t, s, g) {
##     integrate(f = function(theta) {dnorm(x = t, mean = theta, sd = s)*
##                                 dnorm(x = theta, mean = 0, sd = sqrt(g)*s)},
##     lower = 0, upper = Inf)$value
## }
## analytf <- function(t, s, g) {
##     2*dnorm(x = t, mean = 0, sd = s*sqrt(1 + g))*2*pnorm(t/s*sqrt(g/(1 + g)))
## }
## t <- 2
## s <- 0.5
## g <- 1.3
## margf(t, s, g)
## analytf(t, s, g)
## integrate(f = analytf, lower = -Inf, upper = Inf, s = s, g = g)

## non-vectorized version of BFologOR
.BFologOR <- function(ao, bo, nTo = ao + bo, co, do, nCo = co + do, ss) {
    ## check inputs
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

        length(ss) == 1,
        is.numeric(ss),
        is.finite(ss),
        0 <= ss
    )

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

    ## if prior variance is exactly zero, BF is 1
    if (ss == 0) {
        bf0S <- 1
    } else {
        ## compute marginal likelihood under H_S
        mS <- try(stats::integrate(f = function(logOR) {
            vapply(X = logOR, FUN = function(logOR) {
                likExact(a = ao, nT = nTo, c = co, nC = nCo, logOR = logOR) *
                    stats::dnorm(x = logOR, mean = 0, sd = ss)
            }, FUN.VALUE = 1)
        }, lower = -Inf, upper = Inf)$value, silent = TRUE)

        ## compute bf
        if (inherits(mS, "try-error")) return(NaN)
        bf0S <- likExact(a = ao, nT = nTo, c = co, nC = nCo, logOR = 0)/mS
    }
    ## return bf
    return(bf0S)
}

#' @title Bayes factor for logOR effect estimate from original study
#'
#' @description Computes Bayes factor contrasting \deqn{H_0: \log \mathrm{OR} =
#'     0}{H0: logOR = 0} to \deqn{H_1: \log \mathrm{OR} \sim \mathrm{N}(0,
#'     \code{ss}^2)}{H1: logOR ~ N(0, ss^2)} with respect to the data from the
#'     original study (summarized by the entries of a standard 2\eqn{\times}{x}2
#'     table).
#'
#' @param ao Number of cases in treatment group
#' @param bo Number of non-cases in treatment group
#' @param nTo Number of patients in treatment group (specify alternatively to b)
#' @param co Number of cases in control group
#' @param do Number of non-cases in control group
#' @param nCo Number of patients in control group (specify alternatively to d)
#' @param ss Standard devation of the sceptical prior under
#'     \eqn{H_\mathrm{S}}{HS}. Defaults to \code{0}
#'
#' @return Bayes factor (BF \eqn{< 1} indicates that data favour \eqn{H_S}{HS},
#'     while BF \eqn{> 1} indicates that data favour \eqn{H_0}{H0})
#'
#' @author Samuel Pawel
#'
#' @references
#' Pawel, S. and Held, L. (2022). The sceptical Bayes factor for the
#'     assessment of replication success. Journal of the Royal Statistical
#'     Society Series B: Statistical Methodology. 84(3): 879-911.
#'     \doi{10.1111/rssb.12491}
#'
#' @seealso \code{\link{BFo}}, \code{\link{BFoSMD}}
#'
#' @examples
#' ao <- 4
#' bo <- 25
#' co <- 12
#' do <- 24
#'
#' ## compare to normal approximation
#' est <- log(ao*do/bo/co)
#' se <- sqrt(1/ao + 1/bo + 1/co + 1/do)
#' zo <- est/se
#' ssseq <- sqrt(seq(0, 1.5, length.out = 25))
#' comp <- cbind("normal" = BFo(zo = zo, g = ssseq^2/se^2),
#'               "exact" =  BFologOR(ao = ao, bo = bo, co = co, do = do, ss = ssseq))
#' matplot(ssseq, comp, type = "l", lty = 1, lwd = 1.5, log = "y", col = c(1, 2),
#'         xlab = "Prior variance", ylab = bquote(BF["0S"]))
#' legend("topright", c("normal", "exact"), lty = 1, lwd = 1.5, col = c(1,2),
#'        bty = "n")
#'
#' @noRd
BFologOR <- Vectorize(.BFologOR)


## non-vectorized version of BFoSMD
.BFoSMD <- function(to, no, n1o = no, n2o = no, ss,
                    type = c("two.sample", "one.sample", "paired")) {
    ## input checks
    stopifnot(
        length(to) == 1,
        is.numeric(to),
        is.finite(to),

        length(n1o) == 1,
        is.numeric(n1o),
        is.finite(n1o),
        0 < n1o,

        length(n2o) == 1,
        is.numeric(n2o),
        is.finite(n2o),
        0 < n2o,

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
    }

    ## compute df and effective sample size depending on test type
    if (type == "two.sample") {
        df <- n1o + n2o - 2
        nstar <- n1o*n2o/(n1o + n2o)
    } else {
        df <- n1o - 1
        nstar <- n1o
    }

    ## if prior variance is exactly zero, BF is 1
    if (ss == 0) {
        bf0S <- 1
    } else {
        ## compute marginal likelihood under H_S
        mS <- try(stats::integrate(f = function(SMD) {
            suppressWarnings({
                stats::dt(x = to, df = df, ncp = SMD*sqrt(nstar)) *
                    stats::dnorm(x = SMD, mean = 0, sd = ss)
            })
        }, lower = -Inf, upper = Inf)$value, silent = TRUE)

        ## compute bf
        if (inherits(mS, "try-error")) return(NaN)
        bf0S <- stats::dt(x = to, df = df, ncp = 0)/mS
    }
    ## return bf
    return(bf0S)
}

#' @title Bayes factor for SMD effect estimate from original study
#'
#' @description Computes Bayes factor contrasting \deqn{H_0: \mathrm{SMD} =
#'     0}{H0: SMD = 0} to \deqn{H_S: \mathrm{SMD} \sim \mathrm{N}(0,
#'     \code{ss}^2)}{HS: SMD ~ N(0, ss^2)} with respect to the data from the
#'     original study (summarized by \eqn{t}-statistic from \eqn{t}-test and the
#'     corresponding sample size).
#'
#' \eqn{t}-statistics from the following types of \eqn{t}-tests are
#' accepted:
#'
#' - Two-sample \eqn{t}-test where the SMD represents the standardized
#' mean difference between two group means (assuming equal variances in
#' both groups)
#' - One-sample \eqn{t}-test where the SMD represents the standardized
#' mean difference to the null value.
#' - Paired \eqn{t}-test where the SMD represents the standardized mean
#' difference score.
#'
#'
#' @param to \eqn{t}-statistic
#' @param no Sample size (per group)
#' @param n1o Sample size in group 1 (only required for two-sample \eqn{t}-test
#'     with unequal group sizes)
#' @param n2o Sample size in group 2 (only required for two-sample \eqn{t}-test
#'     with unequal group sizes)
#' @param ss Standard devation of the sceptical prior under
#'     \eqn{H_\mathrm{S}}{HS}. Defaults to \code{0}
#' @param type Type of \eqn{t}-test associated with \eqn{t}-statistic. Can be
#'     `"two.sample"`, `"one.sample"`, `"paired"`. Defaults to `"two.sample"`.
#'
#' @return Bayes factor (BF \eqn{< 1} indicates that data favour \eqn{H_S}{HS},
#'     while BF \eqn{> 1} indicates that data favour \eqn{H_0}{H0})
#'
#' @author Samuel Pawel
#'
#' @references Pawel, S. and Held, L. (2022). The sceptical Bayes factor for the
#'     assessment of replication success. Journal of the Royal Statistical
#'     Society Series B: Statistical Methodology. 84(3): 879-911.
#'     \doi{10.1111/rssb.12491}
#'
#' @seealso \code{\link{BFo}}, \code{\link{BFologOR}}
#'
#' @examples
#' to <- 2.5
#' n1 <- 8
#' n2 <- 10
#'
#' ## compare to normal approximation
#' est <- to*sqrt(1/n1 + 1/n2)
#' se <- sqrt(1/n1 + 1/n2)
#' z <- est/se
#' ssseq <- sqrt(seq(0, 1.5, length.out = 25))
#' comp <- cbind("normal" = BFo(zo = z, g = ssseq^2/se^2),
#'               "exact" =  BFoSMD(to = to, n1o = n1, n2o = n2, ss = ssseq))
#' matplot(ssseq, comp, type = "l", lty = 1, lwd = 1.5, log = "y", col = c(1, 2),
#'          xlab = "Prior standard deviation", ylab = bquote(BF["0S"]))
#' legend("topright", c("normal", "exact"), lty = 1, lwd = 1.5, col = c(1,2),
#'        bty = "n")
#'
#' @noRd
BFoSMD <- Vectorize(.BFoSMD)

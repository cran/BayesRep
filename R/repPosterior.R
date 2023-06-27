#' @title Effect size posterior distribution
#' 
#' @description Computes the posterior distribution of the effect size based on
#'     the original and replication effect estimates and their standard errors,
#'     assuming a common underlying effect size and an initial flat prior.
#' 
#' @param to Original effect estimate
#' @param so Standard error of the original effect estimate
#' @param tr Replication effect estimate
#' @param sr Standard error of the replication effect estimate
#' @param lower Lower bound of range for which distribution should computed.
#'     Defaults to minimum of \code{to} and \code{tr} minus four times the
#'     pooled standard error
#' @param upper Upper bound of range for which distribution should computed.
#'     Defaults to maximum of \code{to} and \code{tr} plus four times the pooled
#'     standard error
#' @param nGrid Number of grid points. Defaults to \code{1000}
#' @param plot Logical indicating whether posterior distribution should be
#'     plotted. If \code{FALSE}, only data used for plotting are returned.
#'     Defaults to \code{TRUE}
#' @param CI Logical indicating whether 95% highest posterior credible interval
#'     should be plotted. Defaults to \code{TRUE}
#' @param ... Additional arguments passed to \code{matplot}
#' 
#' @return Plots posterior distribution of the effect size, invisibly returns a
#'     list with the data for the plot
#' 
#' @author Samuel Pawel
#'
#' @examples
#' ## Example from Reproducibility Project Cancer Biology
#' ## Aird: Data from https://elifesciences.org/articles/21253 Fig4B
#' hro <- 25.93
#' lhro <- log(hro)
#' hroCI <- c(5.48, 122.58)
#' se_lhro <- diff(log(hroCI))/(2*qnorm(0.975))
#' hrr <- 3.75
#' lhrr <- log(hrr)
#' hrrCI <- c(1.19, 11.81)
#' se_lhrr <- diff(log(hrrCI))/(2*qnorm(0.975))
#' repPosterior(to = lhro, so = se_lhro, tr = lhrr, sr = se_lhrr)
#'
#'  
#' @export 
repPosterior <- function(to, so, tr, sr,
                         lower = min(c(to, tr)) - 4/sqrt(1/so^2 + 1/sr^2),
                         upper = max(c(to, tr)) + 4/sqrt(1/so^2 + 1/sr^2),
                         nGrid = 1000, plot = TRUE, CI = TRUE, ...) {
    ## input checks
    stopifnot(
        length(to) == 1,
        is.numeric(to),
        is.finite(to),

        length(tr) == 1,
        is.numeric(tr),
        is.finite(tr),

        length(so) == 1,
        is.numeric(so),
        is.finite(so),
        0 < so,

        length(sr) == 1,
        is.numeric(sr),
        is.finite(sr),
        0 < sr,

        length(lower) == 1,
        is.numeric(lower),
        is.finite(lower),
        length(upper) == 1,
        is.numeric(upper),
        is.finite(upper),
        lower < upper,

        length(nGrid) == 1,
        is.numeric(nGrid),
        is.finite(nGrid),
        0 < nGrid,

        length(plot) == 1,
        is.logical(plot),
        !is.na(plot),

        length(CI) == 1,
        is.logical(CI),
        !is.na(CI)
    )

    ## Define appropriate range
    x <- seq(from = lower, to = upper, length.out = nGrid)

    ## Compute posterior
    s2Post <- 1/(1/so^2 + 1/sr^2)
    muPost <- s2Post*(to/so^2 + tr/sr^2)
    hpdLower <- muPost - stats::qnorm(p = 0.975)*sqrt(s2Post)
    hpdUpper <- muPost + stats::qnorm(p = 0.975)*sqrt(s2Post)

    ## Compoute prior, likelihood, and posterior
    prior <- function(x) {
        d <- stats::dnorm(x = x, mean = to, sd = so)
        return(d)
    }
    likelihood <- function(x) {
        stats::dnorm(x = tr, mean = x, sd = sr)
    }
    posterior <- function(x) {
        d <- stats::dnorm(x = x, mean = muPost, sd = sqrt(s2Post))
        return(d)
    }
    priorDF <- data.frame(x = x, density = prior(x))
    posteriorDF <- data.frame(x = x, density = posterior(x))
    likelihoodDF <- data.frame(x = x, density = likelihood(x))

    if (plot == TRUE) {
        graphics::matplot(x = x,
                          y = cbind(posteriorDF$density, priorDF$density,
                                    likelihoodDF$density),
                          type = "l", lty = c(1, 2, 3), col = 1,
                          ylim = c(0, posterior(x = muPost)*1.05),
                          las = 1,
                          xlab = "Effect size", ylab = "Density", ...)
        graphics::legend("topright",
                         legend = c("Posterior", "Original study (prior)",
                                    "Replication study (likelihood)"),
                         lty = c(1, 2, 3), bty = "n")
        if (CI == TRUE) {
            graphics::arrows(x0 = hpdLower, x1 = hpdUpper,
                             y0 = posterior(x = muPost)*1.025,
                             length = 0.1, angle = 90, code = 3)
        }
    }
    out <- list("priorDF" = priorDF, "posteriorDF" = posteriorDF,
                "likelihoodDF" = likelihoodDF, "CI" = c(hpdLower, hpdUpper),
                "prior" = prior, "likelihood" = likelihood, "posterior" = posterior)
    invisible(out)
}

.BFe <- function(to, so, tr, sr, tau, log = FALSE) {
    ## check inputs
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

        length(tau) == 1,
        is.numeric(tau),
        is.finite(tau),
        0 <= tau,

        length(log) == 1,
        is.logical(log)
    )

    ## compute bf
    logbf <- stats::dnorm(x = to - tr, mean = 0, sd = sqrt(so^2 + sr^2),
                          log = TRUE) -
        stats::dnorm(x = to - tr, mean = 0, sd = sqrt(so^2 + sr^2 + 2*tau^2),
                     log = TRUE)
    if (log == TRUE) return(logbf)
    else return(exp(logbf))
}

#' @title Equality of effect size Bayes factor
#'
#' @description Computes the equality of effect size Bayes factor
#' 
#' @details The equality of effect size Bayes factor is the Bayes factor
#'     contrasting the hypothesis of equal original and replication effect sizes
#'     \eqn{H_0: \theta_o = \theta_r}{H0: theta_o = theta_r} to the hypothesis
#'     of unequal effect sizes \eqn{H_1: \theta_o \neq \theta_r}{H1: theta_o !=
#'     theta_r}. Under the hypothesis of unequal effect sizes \eqn{H_1}{H1} the
#'     study specific effect sizes are assumed to be normally distributed around
#'     an overall effect size with heterogeneity standard deviation \code{tau}.
#' 
#' @param to Original effect estimate
#' @param so Standard error of the original effect estimate
#' @param tr Replication effect estimate
#' @param sr Standard error of the replication effect estimate
#' @param tau The heterogeneity standard deviation \eqn{\tau}{tau} under the
#'     hypothesis of unequal effect sizes \eqn{H_1}{H1}
#' @param log Logical indicating whether the natural logarithm of the Bayes
#'     factor should be returned. Defaults to \code{FALSE}
#' 
#' @return The equality of effect size Bayes factor
#'     \eqn{\mathrm{BF}_{01}}{BF01}. \eqn{\mathrm{BF}_{01} > 1}{BF01 > 1}
#'     indicates that the data favour the hypothesis of equal effect sizes
#'     \eqn{H_0}{H0} (replication success), whereas \eqn{\mathrm{BF}_{01} <
#'     1}{BF01 < 1} indicates that the data favour the hypothesis of unequal
#'     effect sizes \eqn{H_1}{H1} (replication failure).
#'
#' 
#' @author Samuel Pawel
#'
#' @references Bayarri, M. and Mayorall, A. (2002). Bayesian Design of
#'     "Successful" Replications. The American Statistician, 56(3): 207-214.
#'     \doi{10.1198/000313002155}
#'
#' Verhagen, J. and Wagenmakers, E. J. (2014). Bayesian tests to quantify the result
#' of a replication attempt. Journal of Experimental Psychology: General,
#' 145:1457-1475. \doi{10.1037/a0036731}
#' 
#' @examples
#' ## strong evidence for unequal effect sizes
#' BFe(to = 1, tr = 0.5, so = sqrt(1/100), sr = sqrt(1/100), tau = 0.3)
#'
#' ## some evidence for equal effect sizes
#' BFe(to = 1, tr = 1, so = sqrt(1/200), sr = sqrt(1/200), tau = 0.3)
#'
#'  
#' @export 
BFe <- Vectorize(.BFe)

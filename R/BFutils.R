## non-vectorized version of formatBF
.formatBF <- function(BF, digits = "default") {
    ## check inputs
    stopifnot(
        length(BF) == 1,
        is.numeric(BF),
        (is.finite(BF) && 0 < BF) || is.na(BF),

        length(digits) == 1,
        (is.character(digits) && digits == "default") ||
        (is.numeric(digits) && 0 <= digits)
    )

    ## return NA if input NA
    if (is.na(BF) || is.nan(BF)) {
        result <- NA
    } else {
        ## format BF
        if (digits == "default") {
            if (BF < 1/1000)
                result <- "< 1/1000"
            if ((BF >= 1/1000) && (BF <= 1/10))
                result <- paste0("1/", as.character(round(1/BF)))
            if ((BF > 1/10) && (BF < 1))
                result <- paste0("1/", as.character(round(1/BF, digits = 1)))
            if ((BF < 10) && (BF >= 1))
                result <- as.character(round(BF, digits = 1))
            if ((BF >= 10) && (BF <= 1000))
                result <- as.character(round(BF))
            if (BF > 1000)
                result <- "> 1000"
        } else {
            if (BF < 1) {
                result <- paste0("1/", as.character(round(1/BF, digits = digits)))
            } else {
                result <- as.character(round(BF, digits = digits))
            }
        }
        ## when 1/1 return 1
        if (result == "1/1") result <- "1"
    }
    return(result)
}

#' @title Formatting of Bayes factors
#'
#' @description Formats Bayes factors such that Bayes factors smaller than 1 are
#'     represented as ratios \eqn{1/x}, where \eqn{x} is rounded to the
#'     specified number of digits, while Bayes factors larger than 1 are only
#'     rounded to the specified number of digits.
#'
#' @param BF Bayes factor
#' @param digits either \code{"default"} (see Details) or a positive integer
#'     specifiying the number of decimal places to round the Bayes factor (for
#'     Bayes factors \eqn{\geq 1}{>= 1}) or its inverse (for Bayes factors
#'     \eqn{< 1}{< 1})
#'
#' @return A character vector of ratios (for inputs \eqn{< 1}{< 1}) or rounded
#'     numeric values (for inputs \eqn{\geq 1}{>= 1}) ).
#'
#' @details The default formatting, which is recommended in Held and Ott (2018),
#'     is as follows: For very small Bayes factors BF < 1/1000, "< 1/1000" is
#'     returned. Bayes factors BF with 1/1000 \eqn{\leq}{<=} BF \eqn{\leq}{<=}
#'     1/10 are formatted as \eqn{1/x} where \eqn{x} is an integer and Bayes
#'     factors BF with \eqn{1/10} \eqn{<} BF \eqn{<} 1 as \eqn{1/x}, where
#'     \eqn{x} is rounded to one decimal place. Accordingly, Bayes factors
#'     \eqn{\leq}{<=} BF \eqn{<} 10 are rounded to one decimal place, Bayes
#'     factors 10 \eqn{\leq}{<=} BF \eqn{\leq}{<=} 1000 are rounded to the next
#'     integer and for larger Bayes factors, "> 1000" is returned.
#'
#' If digits is specified, the Bayes factor (if it is \eqn{\geq}{>=} 1) or its
#' inverse (if the Bayes factor is \eqn{<} 1) is rounded to the number of
#' decimal places specified and returned as a ratio if the Bayes factor is
#' \eqn{<} 1.
#'
#' @references Held, L. and Ott, M. (2018). On \eqn{p}-values and Bayes factors.
#'     Annual Review of Statistics and Its Application, 5, 393-419.
#'     \doi{10.1146/annurev-statistics-031017-100307}
#'
#' @author Manuela Ott (creator of package \code{pCalibrate}), Leonhard Held
#'     (contributor of package \code{pCalibrate}), Samuel Pawel (made small
#'     changes to \code{pCalibrate::formatBF})
#'
#' @examples
#' (bf <- BFr(to = 2, so = 0.5, tr = 2.5, sr = 0.9))
#' formatBF(BF = bf)
#'
#' @export
formatBF <- Vectorize(.formatBF)


#' @title Density of truncated normal distribution
#'
#' @description Computes density of normal distribution truncated to interval
#' \eqn{[a, b]}{`[`a, b`]`}
#'
#' @param x Quantile
#' @param mean Mean
#' @param sd Standard deviation
#' @param a Lower truncation bound
#' @param b Upper truncation bound
#' @param log Logical indicating whether natural logarithm of density
#' should be returned
#'
#' @return Numeric vector of (log) densities
#'
#' @author Samuel Pawel
#'
#' @examples
#' ## verify that density integrates to one
#' stats::integrate(f = dnormtrunc, lower = -0.5, upper = 2, a = -0.5, b = 2)
#'
#' @noRd
dnormtrunc <- function(x, mean = 0, sd = 1, a = -Inf, b = Inf, log = FALSE) {
  ## compute normalizing constant k and indicator I
  k <- diff(stats::pnorm(q = c(a, b), mean = mean, sd = sd))
  I <- as.numeric(x >= a & x <= b)

  ## compute and return truncated density
  dens <- stats::dnorm(x = x, mean = mean, sd = sd)*I/k
  if (log == TRUE) dens <- log(dens)
  return(dens)
}

## non-vectorized version of vss
.vss <- function(x, gamma, jeffreys = FALSE) {
    ## input checks
    stopifnot(length(x) == 1,
              is.numeric(x),
              is.finite(x),

              length(gamma) == 1,
              is.numeric(gamma),
              is.finite(gamma),
              0 < gamma, gamma <= 1,

              length(jeffreys) == 1,
              is.logical(jeffreys))

    ## if |x| <= 1, vss exists only for gamma = 1
    if (abs(x) <= 1) {
      if (gamma == 1) res <- 0
      else res <- NaN
    } else {
      ## compute sufficiently sceptical relative prior variance with W-function
      y <- -x^2*exp(-x^2)/gamma^2
      if (jeffreys == FALSE) {
        res <- as.numeric(-x^2/lamW::lambertWm1(x = y) - 1)
      } else {
        res <- as.numeric(-x^2/lamW::lambertW0(x = y) - 1)
      }
      if (!is.nan(res) && res < 0) res <- NaN
    }
    return(res)
}
#' @title Sufficiently sceptical relative prior variance
#'
#' @description Computes sufficiently sceptical relative prior variance at level
#'     \eqn{\gamma}{gamma}. It is defined as the relative prior variance
#'     \eqn{g_\gamma}{g_gamma} such that the Bayes factor contrasting \deqn{H_0:
#'     \theta = 0}{H0: theta = 0} to \deqn{H_1: \theta \sim \mathrm{N}(0,
#'     g_\gamma \cdot \sigma^2)}{H1: theta ~ N(0, g_gamma*sigma^2)} for data
#'     \eqn{x \, | \, \theta \sim \mathrm{N}(\theta, \sigma^2)}{ x | theta ~
#'     N(theta, sigma^2)} with known variance \eqn{\sigma^2}{sigma^2}, is fixed
#'     at the specified threshold \eqn{\gamma}{gamma}, i.e.
#'     \deqn{\mathrm{BF}_{01}(x) = \frac{f(x \, | \, H_0)}{f(x \, | \, H_1)} =
#'     \gamma.}{ BF01(x) = f(x|H0)/f(x|H1) = gamma.}
#'
#'  If the sufficiently
#'     sceptical relative prior variance exists, there are always two solutions
#'     \eqn{g_\gamma}{g_gamma} and \eqn{g^\prime_\gamma}{g'_gamma}, however,
#'     only the smaller of the two is usually of interest, as the second
#'     solution merely exists due to the Jeffreys-Lindley paradox. If desired,
#'     also the larger solution \eqn{g^\prime_\gamma}{g'_gamma} can be computed
#'     with the argument `jeffreys` set to `TRUE`.
#'
#' @param x Observed data value
#' @param gamma Bayes factor threshold \eqn{\gamma}{gamma}
#' @param jeffreys Logical indicating whether Jeffreys-Lindley solution
#' should be returned. Default is `FALSE`
#' @return Sufficiently sceptical relative prior variance (relative to
#' variance of the data) if existant, otherwise NaN
#'
#'
#' @references Pawel, S. and Held, L. (2022). The sceptical Bayes factor for the
#'     assessment of replication success. Journal of the Royal Statistical
#'     Society Series B: Statistical Methodology, 84(3): 879-911.
#'     \doi{10.1111/rssb.12491}
#'
#' @author Samuel Pawel
#'
#' @examples
#' zo <- 3
#' gamma <- 1/10
#' gss <- vss(x = zo, gamma = gamma)
#' g <- seq(0, 4, 0.01)
#'
#' plot(g, BFo(zo = zo, g = g), type = "l", log = "y", yaxt = "n", ylab = "BF")
#' bf_breaks <- c(100, 30, 10, 3, 1, 1/3, 1/10, 1/30, 1/100)
#' axis(side = 2, at = bf_breaks, labels = formatBF(bf_breaks), las = 1)
#' abline(h = gamma, lty = 2)
#' segments(x0 = gss, y0 = 0.001, y1 = gamma, lty = 2, col = 2)
#'
#' @noRd
vss <- Vectorize(.vss)

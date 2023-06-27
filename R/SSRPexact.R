#' Data from the Social Sciences Replication Project
#'
#' @description Data from the Social Sciences Replication Project.
#' The variables are as follows:
#' \describe{
#'   \item{`study`}{Authors, year, and journal of the original study}
#'   \item{`type`}{Type of effect size. Either `"logOR"` for log oddds ratio effect
#' size, `"SMD1"` for standardized mean difference from one-sample or paired
#' \eqn{t}-test, or `"SMD2"` for standardized mean difference from two-sample
#' \eqn{t}-test}
#'   \item{`to`}{\eqn{t}-statistic from the original study
#' (only available for `"SMD1"` and `"SMD2"`)}
#'   \item{`n1o`}{Sample size in group 1 of the original study
#' (only available for `"SMD1"` and `"SMD2"`)}
#'   \item{`n2o`}{Sample size in group 2 of the original study
#' (only available for `"SMD2"`)}
#'   \item{`tr`}{\eqn{t}-statistic from the replication study
#' (only available for `"SMD1"` and `"SMD2"`)}
#'   \item{`n1r`}{Sample size in group 1 of the replication study
#' (only available for `"SMD1"` and `"SMD2"`)}
#'   \item{`n2r`}{Sample size in group 2 of the replication study
#'    (only available for `"SMD2"`)}
#'   \item{`ao`}{Number of cases in original study treatment group
#' (only available for `"logOR"`)}
#'   \item{`bo`}{Number of non-cases in original study treatment group
#' (only available for `"logOR"`)}
#'   \item{`co`}{Number of cases in original study control group
#' (only available for `"logOR"`)}
#'   \item{`do`}{Number of non-cases in original study control group
#' (only available for `"logOR"`)}
#'   \item{`ar`}{Number of cases in replication study treatment group
#' (only available for `"logOR"`)}
#'   \item{`br`}{Number of cases in replication study control group
#' (only available for `"logOR"`)}
#'   \item{`cr`}{Number of cases in replication study control group
#' (only available for `"logOR"`)}
#'   \item{`dr`}{Number of non-cases in replication study control group
#' (only available for `"logOR"`)}
#' }
#' 
#' @name SSRPexact
#'
#' @docType data
#'
#' @author Samuel Pawel
#'
#' @usage data(SSRPexact)
#'
#' @format A data frame with 21 rows and 16 variables
#'
#' @source The data were manually extracted from the Bayesian supplement of the
#'     SSRP (<https://osf.io/nsxgj/>). The data are licensed under CC0 1.0
#'     Universal.
#'
#' @references Camerer, C. F., Dreber, A., Holzmeister, F., Ho, T.-H., Huber,
#'     J., Johannesson, M., ... Wu, H. (2018). Evaluating the replicability of
#'     social science experiments in Nature and Science between 2010 and 2015.
#'     Nature Human Behaviour, 2, 637-644. \doi{10.1038/s41562-018-0399-z}
#' 
#' 
#' @keywords data
"SSRPexact"

## ## data extracted from Bayesian supplement of SSRP (https://osf.io/nsxgj/)
## study <- c("Ackerman et al. (2010), Science",
##            "Aviezer et al. (2012), Science",
##            "Balafoutas and Sutter (2012), Science",
##            "Derex et al. (2013), Nature",
##            "Duncan et al. (2012), Science",
##            "Gervais and Norenzayan (2012), Science",
##            "Gneezy et al. (2014), Science",
##            "Hauser et al. (2014), Nature",
##            "Janssen et al. (2010), Science",
##            "Karpicke and Blunt (2011), Science",
##            "Kidd and Castano (2013), Science",
##            "Kovacs et al. (2010), Science",
##            "Lee and Schwarz (2010), Science",
##            "Morewedge et al. (2010), Science",
##            "Nishi et al. (2015), Nature",
##            "Pyc and Rawson (2010), Science",
##            "Ramirez and Beilock (2011), Science",
##            "Rand et al. (2012), Nature",
##            "Shah et al. (2012), Science",
##            "Sparrow et al. (2011), Science",
##            "Wilson et al. (2014), Science")
## type <- c("SMD2",  "SMD1", "logOR",
##           NA, # Derex: H0/H1: equal/ordered probabilities
##           "SMD1", "SMD2", "logOR", "logOR",
##           NA, # Janssen: non-parametric Mann-Whitney test
##           "SMD2", "SMD2", "SMD1", "SMD2", "SMD2", "SMD2", "SMD2", "SMD2",
##           "SMD2", "SMD2", "SMD1", "SMD2")
## ## results from studies with t-tests (one-sample or two-sample)
## to <- c(2.02, 13.07, NA, NA, 3.41, 2.24, NA, NA, NA, 4.65, 2.53, 2.42, 2.6, 2.78,
##         2.68, 2.37, 5.53, 2.45, 2.04, 3.26, 4.83)
## n1o <- c(26, 15, NA, NA, 15, 26, NA, NA, NA, 20, 43, 24, 21, 16, 10, 18, 10, 175,
##          26, 69, 15)
## n2o <- c(28, NA, NA, NA, NA, 31, NA, NA, NA, 20, 43, NA, 19, 16, 10, 18, 10, 168,
##          30, NA, 15)
## tr <- c(1.5351, 5.342, NA, NA, 4.6276, -0.82, NA, NA, NA, 2.8825, -0.726, 7.0152,
##         -0.78, 3.54, 2.53, 2.64, -0.7352, 1.01, -0.373, 0.7579, 4.49)
## n1r <- c(296, 14, NA, NA, 92, 262, NA, NA, NA, 23, 349, 95, 147, 44, 24, 156, 45,
##          1058, 298, 234, 20)
## n2r <- c(303, NA, NA, NA, NA, 269, NA, NA, NA, 26, 365, NA, 139, 45, 24, 150, 34,
##          1078, 321, NA, 19)
## ## results from studies with binary outcomes
## ao <- c(NA, NA, 21, NA, NA, NA, 65, 20, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
##         NA, NA, NA)
## bo <- c(NA, NA, 15, NA, NA, NA, 26, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
##         NA, NA, NA)
## co <- c(NA, NA, 11, NA, NA, NA, 43, 4, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
##         NA, NA, NA)
## do <- c(NA, NA, 25, NA, NA, NA, 44, 16, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
##         NA, NA, NA)
## ar <- c(NA, NA, 63, NA, NA, NA, 147, 11, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
##         NA, NA, NA)
## br <- c(NA, NA, 60, NA, NA, NA, 55, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
##         NA, NA, NA)
## cr <- c(NA, NA, 44, NA, NA, NA, 113, 2, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
##         NA, NA, NA)
## dr <- c(NA, NA, 76, NA, NA, NA, 92, 9, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
##         NA, NA, NA)
## SSRPexact <- data.frame(study = study, type = type,
##                         to, n1o, n2o, tr, n1r, n2r,
##                         ao, bo, co, do, ar, br, cr, dr)
## ## save as rds and use version = 2, otherwise the pkg will depend on R > 3.5
## save(SSRPexact, file = "../data/SSRPexact.rda", version = 2)

#' robsurvey: Robust survey statistics.
#'
#' The package robsurvey is a collection of functions for robust survey
#' statistics.
#'
#' @section robsurvey functions:
#' robust Horvitz-Thompson M-estimator of mean and total in \code{svymean_huber()}
#' and \code{svytotal_huber()}, robust trimmed Horvitz-Thompson estimator of mean
#' and total in \code{svymean_trimmed()} and \code{svytotal_trimmed()}, robust winsorized
#' Horvitz-Thompson estimator of mean and total in \code{svymean_winsorized()} and
#' \code{svytotal_winsorized()}, weighted median estimator in \code{weighted_median()},
#' weighted quantile estimator in \code{weighted_quantile()}, weighted median
#' absolute deviation in \code{weighted_mad()}, weighted mean and total
#' estimators in \code{weighted_mean()} and \code{weighted_total()}.
#'
#' @references
#' Hulliger, B. (1995). Outlier Robust Horvitz-Thompson Estimators.
#' Survey Methodology, 21, 79 - 87.
#'
#' Hulliger, B. (2011). Main Results of the AMELI Simulation Study
#' on Advanced Methods for Laeken Indicators. In Proceedings of
#' NTTS2011, Brussels.
#'
#' @docType package
#' @name robsurvey
NULL



#' Weighted median with weighted interpolation
#'
#' \code{weighted.median} computes a weighted median where the
#' exact location corresponds exactly to a cumulative weight of 0.5.
#' This yields a symmetric median.
#'
#' TBD
#'
#' @param x a numeric vector whose weighted sample median is wanted
#' @param w a numeric vector of weights
#' @param na.rm a logical value indicating whether \code{NA} values should be
#'        stripped before the computation proceeds.
#' @return weighted sample median
#' @examples
#' x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
#' weighted_median(x, x)
#' @seealso \code{\link{weighted_quantile}}
#' @export weighted_median
weighted_median <- function(x, w, na.rm=FALSE) {
  # checking and dealing with missingness
  if(missing(w)) stop('Argument w (weights) is missing, with no default.')
  sel <- is.finite(x) & is.finite(w)
  if (sum(sel) < length(x) & na.rm==FALSE) {
    return("There are missing values in x or w. \n")
  }
  x <- x[sel]
  w <- w[sel]

  ord <- order(x)
  w <- w[ord]
  x <- x[ord]

  n <- length(x)
  cumw <- cumsum(w) / sum(w)

  # warning if one of the weights accounts for half of total weight
  if(max(w)/sum(w) > 0.5) {
    cat("\n Dominance of one observation!\n")
  }

  lower.k <- max(which(cumw <= 0.5))
  upper.k <- min(which(cumw > 0.5))

  if (cumw[lower.k] < 0.5) return(x[upper.k])
  else return(0.5 * x[lower.k] + 0.5 * x[upper.k])

}



#' Weighted lower sample quantiles
#'
#' \code{weighted_quantile} computes the weighted lower sample quantile
#'
#' Weighted lower quantiles are computed using an algorithm with \eqn{O(n*log(n))}
#' in worst-case time. There exist superior algorithms; see Cormen et al.
#' (2009, Problem 9.2).
#'
#' @param x a numeric vector whose weighted sample quantiles are wanted
#' @param w a numeric vector of weights
#' @param probs a numeric vector of probabilities with values in [0,1]
#' @param na.rm a logical value indicating whether \code{NA} values should be
#'        stripped before the computation proceeds.
#' @return Weighted sample quantiles
#' @examples
#' x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
#' weighted_quantile(x, x, probs = c(0.25, 0.5, 0.75))
#' @references Cormen,T.H., Leiserson, C.E., Rivest, R.L., and Stein, C. (2009):
#'    Introduction to Algorithms, 3rd ed., Cambridge: MIT Press.
#' @seealso \code{\link{weighted_median}}
#' @export weighted_quantile
#' @importFrom stats na.omit
#' @useDynLib robsurvey wquantile
weighted_quantile <- function(x, w, probs, na.rm = FALSE){
   if(missing(w)) stop('Argument w (weights) is missing, with no default.')
   if (is.factor(x) || is.factor(w) || is.data.frame(x)){
      stop("Arguments 'x' and 'w' must be numeric vectors\n")
   }
   n <- length(x); nw <- length(w)
   if (nw != n) stop("Vectors 'x' and 'w' are not of the same dimension\n")
   if (n == 0){
      return(NA)
   }
   dat <- cbind(x, w)
   if (na.rm){
      dat <- na.omit(dat)
   }else if(any(is.na(dat))) {
      return(NA)
   }
   if (any(probs < 0) | any(probs > 1)){
      stop("Argument 'probs' must satisfy (elementwise): 0 <= probs <= 1!\n")
   }
   res <- NULL
   for (i in 1:length(probs)){
      tmp <- .C("wquantile", x = as.double(dat[, 1]), w = as.double(dat[, 2]),
                probs = as.double(probs[i]), q = as.double(numeric(1)),
                n = as.integer(n))
      res <- c(res, tmp$q)
   }
   return(res)
}



#' Weighted median absolute deviation from the median (MAD)
#'
#' \code{weighted_mad} computes weighted median absolute deviation from the
#' weighted median
#'
#' The weighted MAD is computed as the (normalized) weighted median of the
#' absolute deviation from the weighted median; the median is computed as the
#' weighted lower sample median (see \code{\link{weighted_median}}); the MAD
#' is normalized to be an unbiased estimate of scale at the Gaussian core model.
#'
#' @param x a numeric vector
#' @param w a numeric vector of weights
#' @param na.rm a logical value indicating whether \code{NA} values should be
#'        stripped before the computation proceeds.
#' @param constant (scale factor, default: 1.4826)
#' @return Weighted median absolute deviation from the (weighted) median
#' @examples
#' x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
#' weighted_mad(x, x)
#' @seealso \code{\link{weighted_median}}
#' @export weighted_mad
#' @importFrom stats na.omit
#' @useDynLib robsurvey wmad
weighted_mad <- function(x, w, na.rm = FALSE, constant = 1.4826){
   if(missing(w)) stop('Argument w (weights) is missing, with no default.')
   if (is.factor(x) || is.factor(w) || is.data.frame(x)){
      stop("Arguments 'x' and 'w' must be numeric vectors\n")
   }
   n <- length(x); nw <- length(w)
   if (nw != n) stop("Vectors 'x' and 'w' are not of the same dimension\n")
   if (n <= 1){
      return(NA)
   }
   dat <- cbind(x, w)
   if (na.rm){
      dat <- na.omit(dat)
   }else if(any(is.na(dat))) {
      return(NA)
   }
   tmp <- .C("wmad", x = as.double(dat[, 1]), w = as.double(dat[, 2]),
             mad = as.double(numeric(1)), n = as.integer(n))
   return(tmp$mad * constant)
}



#' Control function for M-estimation (tuning parameters etc.)
#'
#' This function is called internally.
#'
#' Tuning parameters for \code{\link{weighted_mean_huber}},
#' \code{\link{weighted_total_huber}}, \code{\link{svymean_huber}},
#' \code{\link{svytotal_huber}}.
#'
#' @param acc numeric tolerance, stoping rule in the iterative
#' updating scheme (default: \code{1e-5})
#' @param maxit maximum number of updating iterations
#' @param psi psi-function (\code{Huber} or \code{asymHuber})
#' @param ... additional arguments
#' @return List
#' @export rht_control
rht_control <- function(acc = 1e-5, maxit = 100, psi = "Huber", ...){
   if(!(psi %in% c("Huber", "asymHuber"))) stop("Function 'psi' must be
      either 'Huber' or 'asymHuber'\n")
   psi0 <- switch(psi,
      "Huber" = 0L,
      "asymHuber" = 1L)
   list(acc = unname(acc), maxit = unname(maxit), psi = unname(psi0))
}



#' @name wgtmeantotal
#' @aliases weighted_total
#' @aliases weighted_mean
#'
#' @title Weighted total and mean (Horvitz-Thompson and Hajek estimators)
#'
#' @description TBD
#'
#' @details TBD
#'
#' @note \code{wgtmeantotal} is a generic name for the functions documented.
#'
#' @param x a numeric vector
#' @param w a numeric vector of weights
#' @param na.rm a logical value indicating whether \code{NA} values
#' should be stripped before the computation proceeds.
#' @return Estimate (scalar)
#'
#' @rdname wgtmeantotal
#' @examples
#' x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
#' weighted_total(x, x)
#' @export weighted_total
#' @importFrom stats na.omit
weighted_total <- function(x, w, na.rm = FALSE){
   if(missing(w)) stop('Argument w (weights) is missing, with no default.')
   if (is.factor(x) || is.factor(w) || is.data.frame(x)){
      stop("Arguments 'x' and 'w' must be numeric vectors\n")
   }
   n <- length(x); nw <- length(w)
   if (nw != n) stop("Vectors 'x' and 'w' are not of the same dimension\n")
   if (n == 0){
      return(NA)
   }
   dat <- cbind(x, w)
   if (na.rm){
      dat <- na.omit(dat)
   }else if(any(is.na(dat))) {
      return(NA)
   }
   return(sum(dat[, 1] * dat[, 2]))
}
#' @rdname wgtmeantotal
#' @examples
#' x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
#' weighted_mean(x, x)
#' @export weighted_mean
#' @importFrom stats na.omit
weighted_mean <- function(x, w, na.rm = FALSE){
   if(missing(w)) stop('Argument w (weights) is missing, with no default.')
   if (is.factor(x) || is.factor(w) || is.data.frame(x)){
      stop("Arguments 'x' and 'w' must be numeric vectors\n")
   }
   n <- length(x); nw <- length(w)
   if (nw != n) stop("Vectors 'x' and 'w' are not of the same dimension\n")
   if (n == 0){
      return(NA)
   }
   dat <- cbind(x, w)
   if (na.rm){
      dat <- na.omit(dat)
   }else if(any(is.na(dat))) {
      return(NA)
   }
   return(sum(dat[ ,1] * dat[, 2]) / sum(dat[, 2]))
}



#' @name huberwgt
#' @aliases weighted_mean_huber
#' @aliases weighted_total_huber
#' @aliases svymean_huber
#' @aliases svytotal_huber
#'
#' @title Huber M-estimators of the weighted mean and weighted total
#'
#' @description Weighted Huber M-estimators of the mean and total are available in two forms:
#' \itemize{
#'    \item \strong{bare-bone} functions: \code{weighted_mean_huber} and
#'	 \code{weighted_total_huber},
#'    \item estimation \strong{methods}: \code{svymean_huber} and
#'	 \code{svytotal_huber} (incl. variance estimation
#'	 based on the functionality of the \pkg{survey} package).
#' }
#'
#' @details \describe{
#'    \item{\emph{Overview}}{
#'    Robust M-estimator of the Horvitz--Thompson total or the Hajek mean
#'	 \itemize{
#'	 \item bare-bone functions: return the estimate (no variance estimation)
#'	 \item estimation methods on the basis of \pkg{survey} (incl. variance estimation)
#'	 }
#'    }
#'    \item{\emph{Type}}{
#'    Two \code{type}s of estimation methods are available:
#'	 \describe{
#'	    \item{\code{rht}}{(robust) Horvitz-Thompson M-estimator of the
#'	    total/mean
#'	    }
#'	    \item{\code{rwm}}{(robust) weighted mean estimator of
#'	    a Hajek-type estimator of the mean.
#'	    }
#'	 }
#'    If the study variable \code{x} is positively correlated with the inclusion
#'    probabilities, type \code{"rht"} tends to be superior.
#'    }
#'    \item{\emph{Scale}}{
#'    M-estimators of location are not scale invariant. The unkown scale is
#'    estimated simultaneously with the estimate of location (mean or toal) as
#'    the weighted median absolute deviation from the weighted median (MAD, see
#'    \code{\link{weighted_mad}}).
#'
#'    }
#'    \item{\emph{Variance}}{
#'    Variance estimates of the mean or total estimator are computed as first-order
#'    linearization using the design-based-estimation capabilities available
#'    in package \pkg{survey}.
#'    }
#'    \item{\emph{Tuning}}{
#'    Additional arguments can be passed (via \dots) to specify the control
#'    parameters (e.g. number of iterations, psi-function, etc.); see
#'    \code{\link{rht_control}} for details.
#'    }
#'    \item{\emph{Domain estimation}}{
#'    Estimates for domains can be obtained using the \link[survey]{svyby}
#'    wrapper in the \pkg{survey} package (see examples).
#'    }
#' }
#'
#' @section Utility functions:
#' For the methods \code{svymean_huber} and \code{svytotal_huber}, the following
#' utility functions can be used
#' \itemize{
#'    \item \code{summary} gives a summary of the estimation properties
#'    \item \code{\link{robweights}} retrieves the robustness weights
#'    \item \code{coef}, \code{vcov}, \code{residuals}, and \code{fitted}
#'	 retrieve the estimate, variance, residuals and fitted
#'	 values, respectively
#' }
#'
#' @note \code{huberwgt} is a generic name for the functions documented.
#'
#' @param x a numeric vector (\code{weighted.[total/mean].huber} or
#' \code{weighted.[total/mean].huber}); a formula object or variable
#' name (\code{svymean_huber} or \code{svytotal_huber})
#' @param w a numeric vector of weights
#' @param design a \code{survey.design} object (see \link[survey]{svydesign}
#' in \pkg{survey})
#' @param k a robustness tuning constant, \eqn{k} in \eqn{[0, \infty)}
#' @param type type of estimator: \code{"rht"} (default) or \code{"rwm"}
#' @param info logical (default: \code{FALSE}); if \code{TRUE} further
#' estimation details are returned
#' @param na.rm a logical value indicating whether \code{NA} values should
#' be stripped before the computation proceeds.
#' @param ... additional arguments passed to the control object
#' (see \code{\link{rht_control}})
#'
#' @return \itemize{
#'    \item An estimate (scalar) for \code{weighted.[total/mean].huber}
#'    (unless \code{info=TRUE})
#'    \item An object of class \code{svystat.rob} for functions of the type
#'    \code{msvy[total/mean]}, i.e. a list including the following components:
#'    \code{characteristic}, \code{estimator}, \code{estimate}, \code{variance},
#'    \code{robust}, \code{optim}, \code{residuals}, \code{model}, \code{design},
#'    and \code{call}.
#' }
#'
#' @references Hulliger, B. (1995). Outlier Robust Horvitz-Thompson Estimators,
#' \emph{Survey Methodology} 21(1): 79-87.
#'
#' @seealso \code{\link{svymean_trimmed}}, \code{\link{svytotal_trimmed}},
#' \code{\link{svymean_winsorized}}, \code{\link{svytotal_winsorized}},
#' \code{\link{weighted_mean_trimmed}}, \code{\link{weighted_total_trimmed}}
#' \code{\link{weighted_mean_winsorized}}, \code{\link{weighted_total_winsorized}}
#'
#' @rdname huberwgt
#' @export weighted_mean_huber
#' @importFrom stats na.omit
#' @useDynLib robsurvey rwlslm
weighted_mean_huber <- function(x, w, k = 1.5, type = "rht", info = FALSE,
   na.rm = FALSE, ...){
   if(missing(w)) stop('Argument w (weights) is missing, with no default.')
   ctrl <- rht_control(...)
   if(!(type %in% c("rht", "rwm"))) stop("Argument 'type' must be either 'rht' or 'rwm'\n")
   n <- length(x)
   if (length(w) != n) stop("Vectors 'x' and 'w' are not of the same dimension\n")
   if (is.factor(x) || is.factor(w) || is.data.frame(x)){
      stop("Arguments 'x' and 'w' must be numeric vectors\n")
   }
   n <- length(x); nw <- length(w)
   if (nw != n) stop("Vectors 'x' and 'w' are not of the same dimension\n")
   if (n == 0){
      return(NA)
   }
   dat <- cbind(x, w)
   if (na.rm){
      dat <- na.omit(dat)
      n <- nrow(dat)
   }else if(any(is.na(dat))) {
      return(NA)
   }
   x <- switch(type,
      "rht" = rep(1, n),
      "rwm" = mean(dat[,2]) / dat[,2])
   # estimate by irwls
   tmp <- .C("rwlslm", x = as.double(x), y = as.double(dat[, 1]),
      w = as.double(dat[, 2]), resid = as.double(numeric(n)),
      infl = as.double(numeric(n)), robwgt = as.double(numeric(n)),
      n = as.integer(n), p = as.integer(1), k = as.double(k),
      beta = as.double(numeric(1)), scale = as.double(numeric(1)),
      maxit = as.integer(ctrl$maxit), tol = as.double(ctrl$acc),
      psi = as.integer(ctrl$psi))
   if(tmp$maxit == 0) cat("IRWLS algorithm did not converge!\n")
   if(info){
      res <- list(characteristic = "mean",
	 estimator = paste0("M-estimator (", type, ", ", ifelse(ctrl$psi == 0,
	    "Huber", "asymHuber"), ", k = ", tmp$k,")"),
	 estimate = tmp$beta,
	 scale = tmp$scale,
	 optim = ifelse(tmp$maxit == 0, paste0("did NOT converge in ",
	 ctrl$maxit ," iterations"), paste0("converged in ", tmp$maxit,
	 " iterations")),
	 robweights = tmp$robwgt)
   }else{
      res <- tmp$beta
   }
   return(res)
}
#' @rdname huberwgt
#' @export weighted_total_huber
weighted_total_huber <- function(x, w, k = 1.5, type = "rht", info = FALSE,
   na.rm = FALSE, ...){
   if(missing(w)) stop('Argument w (weights) is missing, with no default.')
   res <- weighted_mean_huber(x, w, k, type, info, na.rm, ...)
   if(length(res) == 1){
      res <- res * sum(w)
   }else{
      res$characteristic <- "total"
      res$estimate <- res$estimate * sum(w)
   }
   return(res)
}
#' @rdname huberwgt
#' @examples
#' library(survey)
#' data(api)
#' dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
#' svymean_huber(~api00, dstrat, k = 2)
#' # Domain estimates
#' svyby(~api00, by = ~stype, design = dstrat, svymean_huber, k = 1.34)
#' @export svymean_huber
#' @importFrom stats model.frame na.fail
#' @importFrom stats weights
#' @useDynLib robsurvey rwlslm
svymean_huber <- function(x, design, k = 1.5, type = "rht", ...){
   ctrl <- rht_control(...)
   if (class(x) == "formula"){
      mf <- model.frame(x, design$variables, na.action = na.fail)
      n <- nrow(mf)
      if (ncol(mf) > 1) stop("Argument 'y' must be a formula of one single variable")
      yname <- names(mf)
      y <- mf[[1]]
   }else{
      if (is.character(x)){
	 yname <- x
	 y <- design$variables[, x]
	 n <- length(y)
	 if (any(is.na(y))) stop(paste0("Variable '", yname, "' must not contain NA's\n"))
      }else{
	 stop("svymean_huber is not defined for object of class: ", class(x), "\n")
      }
   }
   w <- as.numeric(weights(design))
   x <- switch(type,
      "rht" = rep(1, n),
      "rwm" = mean(w) / w)
   # estimate by irwls
   tmp <- .C("rwlslm", x = as.double(x), y = as.double(y), w = as.double(w),
      resid = as.double(numeric(n)), infl = as.double(numeric(n)),
      robwgt = as.double(numeric(n)), n = as.integer(n), p = as.integer(1),
      k = as.double(k), beta = as.double(numeric(1)),
      scale = as.double(numeric(1)), maxit = as.integer(ctrl$maxit),
      tol = as.double(ctrl$acc), psi = as.integer(ctrl$psi))
   if(tmp$maxit == 0) cat("IRWLS algorithm did not converge!\n")
   names(tmp$beta) <- yname
   # compute variance (using influence function values)
   design$variables$zz <- tmp$infl
   v <- as.numeric(attr(survey::svymean(~zz, design), "var"))
   design$variables$zz <- NULL
   robweights <- tmp$robwgt
   outliers <- 1 * (abs(robweights) < 1)
   res <- list(characteristic = "mean",
      estimator = paste0("M-estimator (", type, ")"),
      estimate = tmp$beta,
      variance = v,
      robust = list(psifunction = ifelse(ctrl$psi == 0, "Huber", "asymHuber"),
	 k = k, robweights = robweights, outliers = outliers, scale = tmp$scale),
      optim = list(converged = (tmp$maxit != 0), niter = ifelse(tmp$maxit == 0,
	 ctrl$maxit, tmp$maxit)),
      residuals = tmp$resid,
      model = list(y = tmp$y, x = tmp$x, w = tmp$w),
      design = design,
      call = match.call())
   class(res) <- "svystat.rob"
   res
}
#' @rdname huberwgt
#' @export svytotal_huber
#' @importFrom stats weights
svytotal_huber <- function(x, design, k = 1.5, ...){
   tmp <- svymean_huber(x, design, k, type = "rht", ...)
   tmp$characteristic <- "total"
   sumw <- sum(weights(design))
   tmp$estimate <- tmp$estimate * sumw
   tmp$variance <- tmp$variance * sumw^2
   tmp
}



#' @name trimwgt
#' @aliases weighted_mean_trimmed
#' @aliases weighted_total_trimmed
#' @aliases svymean_trimmed
#' @aliases svytotal_trimmed
#'
#' @title Weighted trimmed mean and trimmed total
#'
#' @description Weighted trimmed estimators of the mean and total are available in two forms:
#' \itemize{
#'    \item \strong{bare-bone} functions: \code{weighted_mean_trimmed} and
#'	 \code{weighted_total_trimmed},
#'    \item estimation \strong{methods}: \code{svymean_trimmed} and
#'	 \code{svytotal_trimmed} (incl. variance estimation
#'	 based on the functionality of the \pkg{survey} package).
#' }
#'
#' @details \describe{
#'    \item{\emph{Overview}}{
#'    Robust trimmed Horvitz--Thompson total or Hajek mean
#'	 \itemize{
#'	 \item bare-bone functions: return the estimate (no variance estimation)
#'	 \item estimation methods on the basis of \pkg{survey} (incl. variance estimation)
#'	 }
#'    }
#'    \item{\emph{Variance}}{
#'    Variance estimates of the mean or total estimator are computed as first-order
#'    linearization using the design-based-estimation capabilities available
#'    in package \pkg{survey}.
#'    }
#'    \item{\emph{Domain estimation}}{
#'    Estimates for domains can be obtained using the \link[survey]{svyby}
#'    wrapper in the \pkg{survey} package (see examples).
#'    }
#' }
#'
#' @section Utility functions:
#' For the methods \code{svymean_trimmed} and \code{svytotal_trimmed}, the following
#' utility functions can be used
#' \itemize{
#'    \item \code{summary} gives a summary of the estimation properties
#'    \item \code{\link{robweights}} retrieves the robustness weights
#'    \item \code{coef}, \code{vcov}, \code{residuals}, and \code{fitted}
#'	 retrieve, respectively, the estimate, variance, residuals and fitted
#'	 values
#' }
#'
#' @note \code{trimwgt} is a generic name for the functions documented.
#'
#' @param x numeric vector (\code{weighted_mean_trimmed} or \code{weighted_total_trimmed});
#' a formula object or variable name (\code{svymean_trimmed} or \code{svytotal_trimmed})
#' @param w numeric vector of weights
#' @param design a \code{survey.design} object (see \link[survey]{svydesign} in \pkg{survey})
#' @param LB lower bound of trimming, such that \eqn{0 \leq LB < UB \leq 1}
#' @param UB upper bound of trimming, such that \eqn{0 \leq LB < UB \leq 1}
#' @param ... additional arguments (not used)
#' @param na.rm a logical value indicating whether \code{NA} values should be
#' stripped before the computation proceeds.
#'
#' @return Estimate (scalar) or object of class \code{svystat.rob}
#'
#' @seealso \code{\link{svymean_huber}}, \code{\link{svytotal_huber}},
#' \code{\link{svymean_winsorized}}, \code{\link{svytotal_winsorized}},
#' \code{\link{weighted_mean_huber}}, \code{\link{weighted_total_huber}},
#' \code{\link{weighted_mean_winsorized}}, \code{\link{weighted_total_winsorized}}
#'
#' @rdname trimwgt
#' @export weighted_mean_trimmed
#' @importFrom stats na.omit
#' @useDynLib robsurvey wmeantrimmed
weighted_mean_trimmed <- function(x, w, LB = 0.05, UB = 1 - LB, na.rm = FALSE){
   if(missing(w)) stop('Argument w (weights) is missing, with no default.')
   n <- length(x)
   if (length(w) != n) stop("Vectors 'x' and 'w' are not of the same dimension\n")
   if (is.factor(x) || is.factor(w) || is.data.frame(x)){
      stop("Arguments 'x' and 'w' must be numeric vectors\n")
   }
   n <- length(x); nw <- length(w)
   if (nw != n) stop("Vectors 'x' and 'w' are not of the same dimension\n")
   if (n == 0){
      return(NA)
   }
   dat <- cbind(x, w)
   if (na.rm){
      dat <- na.omit(dat)
      n <- nrow(dat)
   }else if(any(is.na(dat))) {
      return(NA)
   }
   if (LB >= UB) stop("Argument 'LB' must be smaller than 'UB'!")
   if (LB < 0) stop("Argument 'LB' must not be < 0!")
   if (UB > 1) stop("Argument 'UB' must not be > 1!")
   tmp <- .C("wmeantrimmed", x = as.double(x), w = as.double(w),
             lb = as.double(LB), ub = as.double(UB), mean = as.double(numeric(1)),
             n = as.integer(n))
   return(tmp$mean)
}
#' @rdname trimwgt
#' @export weighted_total_trimmed
weighted_total_trimmed <- function(x, w, LB = 0.05, UB = 1 - LB, na.rm = FALSE){
   if(missing(w)) stop('Argument w (weights) is missing, with no default.')
   res <- weighted_mean_trimmed(x, w, LB, UB, na.rm)
   if(length(res) == 1){
      res <- res * sum(w)
   }else{
      res$characteristic <- "total"
      res$estimate <- res$estimate * sum(w)
   }
   return(res)
}
#' @rdname trimwgt
#' @examples
#' library(survey)
#' data(api)
#' dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
#' svymean_trimmed(~api00, dstrat, LB = 0.05)
#' # Domain estimates
#' svyby(~api00, by = ~stype, design = dstrat, svymean_trimmed, LB = 0.1)
#' @export svymean_trimmed
#' @importFrom stats model.frame na.fail
#' @importFrom stats weights
svymean_trimmed <- function(x, design, LB = 0.05, UB = 1 - LB, ...){
   if (class(x) == "formula"){
      mf <- model.frame(x, design$variables, na.action = na.fail)
      n <- nrow(mf)
      if (ncol(mf) > 1) stop("Argument 'y' must be a formula of one single variable")
      yname <- names(mf)
      y <- mf[[1]]
   }else{
      if (is.character(x)){
	 yname <- x
	 y <- design$variables[, x]
	 n <- length(y)
	 if (any(is.na(y))) stop(paste0("Variable '", yname, "' must not contain NA's\n"))
      }else{
	 stop("svymean_trimmed is not defined for object of class: ", class(x), "\n")
      }
   }
   w <- as.numeric(weights(design))
   est <- weighted_mean_trimmed(y, w, LB, UB)
   names(est) <- yname
   # compute influence function
   quant <- weighted_quantile(y, w, probs = c(LB, UB))
   below <- floor(LB * n)
   above <- ceiling(UB * n)
   mat <- c(rep((1 - LB) * quant[1] - (1 - UB) * quant[2], below),
      rep(-LB * quant[1] - (1 - UB) * quant[2], (above - below)),
      rep(UB * quant[2] - LB * quant[1], (n - above)))
   if(below != 0){
      y[1:below] <- 0
   }
   if(above != n){
      y[(above + 1):n] <- 0
   }
   infl <- (y + mat) * (1 / (UB - LB)) - est
   # compute variance (using influence function values)
   design$variables$zz <- infl
   v <- as.numeric(attr(survey::svymean(~zz, design), "var"))
   design$variables$zz <- NULL
   # return value
   res <- list(characteristic = "mean",
      estimator = paste0("Weighted trimmed estimator (LB = ", LB, ", UB = ",
	 UB, ")"),
      estimate = est,
      variance = v,
      robust = list(robweights = NULL),
      optim = NULL,
      residuals = y - est,
      model = list(y = y, w = w),
      design = design,
      call = match.call())
   class(res) <- "svystat.rob"
   res
}
#' @rdname trimwgt
#' @export svytotal_trimmed
#' @importFrom stats weights
svytotal_trimmed <- function(x, design, LB = 0.05, UB = 1 - LB, ...){
   tmp <- svymean_trimmed(x, design, LB, UB)
   tmp$characteristic <- "total"
   sumw <- sum(weights(design))
   tmp$estimate <- tmp$estimate * sumw
   tmp$variance <- tmp$variance * sumw^2
   tmp
}



#' @name winswgt
#' @aliases weighted_mean_winsorized
#' @aliases weighted_total_winsorized
#' @aliases svymean_winsorized
#' @aliases svytotal_winsorized
#'
#' @title Weighted winsorized mean and trimmed total
#'
#' @description Weighted winsorized estimators of the mean and total are available in two forms:
#' \itemize{
#'    \item \strong{bare-bone} functions: \code{weighted_mean_winsorized} and
#'	 \code{weighted_total_winsorized},
#'    \item estimation \strong{methods}: \code{svymean_winsorized} and
#'	 \code{svytotal_winsorized} (incl. variance estimation
#'	 based on the functionality of the \pkg{survey} package).
#' }
#'
#' @details \describe{
#'    \item{\emph{Overview}}{
#'    Robust winsorized Horvitz--Thompson total or Hajek mean
#'	 \itemize{
#'	 \item bare-bone functions: return the estimate (no variance estimation)
#'	 \item estimation methods on the basis of \pkg{survey} (incl. variance estimation)
#'	 }
#'    }
#'    \item{\emph{Variance}}{
#'    Variance estimates of the mean or total estimator are computed as first-order
#'    linearization using the design-based-estimation capabilities available
#'    in package \pkg{survey}.
#'    }
#'    \item{\emph{Domain estimation}}{
#'    Estimates for domains can be obtained using the \link[survey]{svyby}
#'    wrapper in the \pkg{survey} package (see examples).
#'    }
#' }
#'
#' @section Utility functions:
#' For the methods \code{svymean_winsorized} and \code{svytotal_winsorized}, the following
#' utility functions can be used
#' \itemize{
#'    \item \code{summary} gives a summary of the estimation properties
#'    \item \code{\link{robweights}} retrieves the robustness weights
#'    \item \code{coef}, \code{vcov}, \code{residuals}, and \code{fitted}
#'	 retrieve, respectively, the estimate, variance, residuals and fitted
#'	 values
#' }
#'
#' @note \code{winswgt} is a generic name for the functions documented.
#'
#' @param x numeric vector (\code{weighted_mean_winsorized} or
#' \code{weighted_total_winsorized}); a formula object or variable name
#' (\code{svymean_winsorized} or \code{svytotal_winsorized})
#' @param w numeric vector of weights
#' @param design a \code{survey.design} object (see \link[survey]{svydesign}
#' in \pkg{survey})
#' @param LB lower bound of winsorizing, such that \eqn{0 \leq LB < UB \leq 1}
#' @param UB upper bound of winsorizing, such that \eqn{0 \leq LB < UB \leq 1}
#' @param ... additional arguments (not used)
#' @param na.rm a logical value indicating whether \code{NA} values should be
#' stripped before the computation proceeds.
#'
#' @return Estimate (scalar) or object of class \code{svystat.rob}
#'
#' @seealso \code{\link{svymean_huber}}, \code{\link{svytotal_huber}},
#' \code{\link{svymean_trimmed}}, \code{\link{svytotal_trimmed}},
#' \code{\link{weighted_mean_huber}}, \code{\link{weighted_total_huber}},
#' \code{\link{weighted_mean_trimmed}}, \code{\link{weighted_total_trimmed}}
#'
#' @rdname winswgt
#' @export weighted_mean_winsorized
#' @importFrom stats na.omit
#' @useDynLib robsurvey wmeanwinsorized
weighted_mean_winsorized <- function(x, w, LB = 0.05, UB = 1 - LB,
   na.rm = FALSE){
   if(missing(w)) stop('Argument w (weights) is missing, with no default.')
   n <- length(x)
   if (length(w) != n) stop("Vectors 'x' and 'w' are not of the same dimension\n")
   if (is.factor(x) || is.factor(w) || is.data.frame(x)){
      stop("Arguments 'x' and 'w' must be numeric vectors\n")
   }
   n <- length(x); nw <- length(w)
   if (nw != n) stop("Vectors 'x' and 'w' are not of the same dimension\n")
   if (n == 0){
      return(NA)
   }
   dat <- cbind(x, w)
   if (na.rm){
      dat <- na.omit(dat)
      n <- nrow(dat)
   }else if(any(is.na(dat))) {
      return(NA)
   }
   if (LB >= UB) stop("Argument 'LB' must be smaller than 'UB'!")
   if (LB < 0) stop("Argument 'LB' must not be < 0!")
   if (UB > 1) stop("Argument 'UB' must not be > 1!")
   tmp <- .C("wmeanwinsorized", x = as.double(x), w = as.double(w),
      lb = as.double(LB), ub = as.double(UB), mean = as.double(numeric(1)),
      n = as.integer(n))
   return(tmp$mean)
}
#' @rdname winswgt
#' @export weighted_total_winsorized
weighted_total_winsorized <- function(x, w, LB = 0.05, UB = 1 - LB, na.rm = FALSE){
   if(missing(w)) stop('Argument w (weights) is missing, with no default.')
   res <- weighted_mean_winsorized(x, w, LB, UB, na.rm)
   if(length(res) == 1){
      res <- res * sum(w)
   }else{
      res$characteristic <- "total"
      res$estimate <- res$estimate * sum(w)
   }
   return(res)
}
#' @rdname winswgt
#' @examples
#' library(survey)
#' data(api)
#' dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
#' svymean_winsorized(~api00, dstrat, LB = 0.05)
#' # Domain estimates
#' svyby(~api00, by = ~stype, design = dstrat, svymean_winsorized, LB = 0.1)
#' @export svymean_winsorized
#' @importFrom stats model.frame na.fail
#' @importFrom stats weights
svymean_winsorized <- function(x, design, LB = 0.05, UB = 1 - LB, ...){
   if (class(x) == "formula"){
      mf <- model.frame(x, design$variables, na.action = na.fail)
      n <- nrow(mf)
      if (ncol(mf) > 1) stop("Argument 'y' must be a formula of one single variable")
      yname <- names(mf)
      y <- mf[[1]]
   }else{
      if (is.character(x)){
	 yname <- x
	 y <- design$variables[, x]
	 n <- length(y)
	 if (any(is.na(y))) stop(paste0("Variable '", yname, "' must not contain NA's\n"))
      }else{
	 stop("svymean_winsorized is not defined for object of class: ", class(x), "\n")
      }
   }
   w <- as.numeric(weights(design))
   est <- weighted_mean_winsorized(y, w, LB, UB)
   names(est) <- yname
   # compute influence function
   # FIXME: variance => density estimate (now, it is trimmed mean variance)
   quant <- weighted_quantile(y, w, probs = c(LB, UB))
   below <- floor(LB * n)
   above <- ceiling(UB * n)
   mat <- c(rep((1 - LB) * quant[1] - (1 - UB) * quant[2], below),
      rep(-LB * quant[1] - (1 - UB) * quant[2], (above - below)),
      rep(UB * quant[2] - LB * quant[1], (n - above)))
   if(below != 0){
      y[1:below] <- 0
   }
   if(above != n){
      y[(above + 1):n] <- 0
   }
   infl <- (y + mat) * (1 / (UB - LB)) - est
   # compute variance (using influence function values)
   design$variables$zz <- infl
   v <- as.numeric(attr(survey::svymean(~zz, design), "var"))
   design$variables$zz <- NULL
   # return value
   res <- list(characteristic = "mean",
      estimator = paste0("Weighted winsorized estimator (LB = ", LB, ", UB = ",
	 UB, ")"),
      estimate = est,
      variance = v,
      robust = list(robweights = NULL),
      optim = NULL,
      residuals = y - est,
      model = list(y = y, w = w),
      design = design,
      call = match.call())
   class(res) <- "svystat.rob"
   res
}
#' @rdname winswgt
#' @export svytotal_winsorized
#' @importFrom stats weights
svytotal_winsorized <- function(x, design, LB = 0.05, UB = 1 - LB, ...){
   tmp <- svymean_winsorized(x, design, LB, UB)
   tmp$characteristic <- "total"
   sumw <- sum(weights(design))
   tmp$estimate <- tmp$estimate * sumw
   tmp$variance <- tmp$variance * sumw^2
   tmp
}



#' Weighted robust line fitting
#'
#' \code{weighted_line} fits a robust line and allows weights.
#'
#' Uses different quantiles for splitting the sample than \code{line()}.
#' Is based on \code{weighted_median()}.
#'
#' @param x a numeric vector (explanatory variable)
#' @param y a numeric vector (response variable)
#' @param w a numeric vector of weights
#' @param iter number of iterations for enhancing the slope
#' @param na.rm a logical value indicating whether rows with \code{NA}
#' values should be stripped before the computation proceeds
#' @return intercept and slope of the fitted line
#' @examples
#' data(cars)
#' weighted_line(cars$speed, cars$dist, w=rep(1, length(cars$speed)))
#' weighted_line(cars$speed, cars$dist, w=rep(1:10, each=5))
#' @seealso \code{\link[stats]{line}}
#' @export weighted_line
#' @importFrom stats model.frame complete.cases
weighted_line <- function(x, y=NULL, w, na.rm=FALSE, iter = 1){

  if(missing(w)) stop('Argument w (weights) is missing, with no default.')

  # quantiles as implemented in line() but with weights
  # x and w sorted according to x
  line.quantile <- function(x, w, prob){
    n <- length(x)
    cumw <- n * cumsum(w)/sum(w)
    low <- min(which(cumw >= (n-1) * prob))
    high <- max(which(cumw <= (n-1) * prob))
    return((x[low+1] + x[high+1])/2)
  }

  if (inherits(x, "formula")) {
    mf <- model.frame(x)
    y <- mf[, 1]
    x <- mf[, -1]
  }
  if (NCOL(x) > 1) return("x contains more than 1 explanatory variables.")
  dat <- data.frame(x, y, w)
  ok <- complete.cases(dat$x, dat$y, dat$w)
  n <- sum(ok)
  if (n < length(x)) {
    if (na.rm) {x <- dat$x[ok]; y <- dat$y[ok]; w <- dat$w[ok]} else
      return("There are missing values in x, y or w.\n")
  }

  # Sort data according to x
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  # standardise weights to n
  w <- n * w[ord] / sum(w)

  # Groups
  lim1 <- line.quantile(x, w, 1/3)
  lim2 <- line.quantile(x, w, 2/3)
  groups <- (x <= lim1) * 1 + (x > lim1 & x < lim2) * 2 + (x >= lim2) * 3

  # Medians for x
  wmedx <- c(weighted_median(x[groups==1], w[groups==1]),
             weighted_median(x[groups==3], w[groups==3]))

  # polishing (affects only the slope)
  slope <- 0; r <- y; j <- 0
  while (j <= iter - 1) {
    wmedr <- c(weighted_median(r[groups==1], w[groups==1]),
               weighted_median(r[groups==3], w[groups==3]))
    slope <- slope + (wmedr[2] - wmedr[1]) / (wmedx[2] - wmedx[1])
    r <- y - slope * x
    j <- j + 1
  }

  # intercept and predicted values
  intercept <- weighted_median(r, w)
  yhat <- intercept + slope * x

  # return
  structure(list(call = sys.call(), coefficients = c(intercept, slope),
                 residuals = y - yhat, fitted.values = yhat), class = "tukeyline")
}



#' Robust simple linear regression based on medians
#'
#' For type medslopes the median individual ratios response/explanatory
#' is used as estimator of the slope. For version ratiomeds the ratio of
#' the median crossproduct to the median of squares of the explanatory
#' variable is used as the estimator of the slope. Survey weights may be used.
#' Missing values are neglected.
#'
#' Uses \code{weighted_median()}. The median of slopes (type="slopes")
#' uses \eqn{b1=M((y-M(y,w))/(x-M(x,w)), w)}. The median of crossproducts
#' by median of squares (type="products") uses
#' \eqn{b1=M((y-M(y,w))(x-M(x,w)), w )/ M((x-M(x,w )^2), w )}, where
#' \eqn{M(x, w)} is shorthand for the function \code{weighted_median(x, w)}.
#' The function allows weights and missing values.
#'
#' @param x a numeric vector (explanatory variable)
#' @param y a numeric vector (response variable)
#' @param w a numeric vector of weights
#' @param type either "slopes" (default) or "products"
#' @param na.rm a logical value indicating whether rows with \code{NA}
#' values should be stripped before the computation proceeds
#' @return a vector with two components: intercept and slope
#' @examples
#' x <- c(1, 2, 4, 5)
#' y <- c(3, 2, 7, 4)
#' weighted_line(y~x, w=rep(1, length(x)))
#' weighted_median_line(y~x, w=rep(1, length(x)))
#' weighted_median_line(y~x, w=rep(1, length(x)), type="prod")
#'
#' data(cars)
#' with(cars, weighted_median_line(dist ~ speed, w=rep(1, length(dist))))
#' with(cars, weighted_median_line(dist ~ speed, w=rep(1, length(dist)), type="prod"))
#'
#' # weighted
#' w <- c(rep(1,20), rep(2,20), rep(5, 10))
#' with(cars, weighted_median_line(dist ~ speed, w=w))
#' with(cars, weighted_median_line(dist ~ speed, w=w, type="prod"))
#'
#' # outlier in y
#' cars$dist[49] <- 360
#' with(cars, weighted_median_line(dist ~ speed, w=w))
#' with(cars, weighted_median_line(dist ~ speed, w=w, type="prod"))
#'
#' # outlier in x
#' data(cars)
#' cars$speed[49] <- 72
#' with(cars, weighted_median_line(dist ~ speed, w=w))
#' with(cars, weighted_median_line(dist ~ speed, w=w, type="prod"))
#' @seealso \code{\link[stats]{line}}, \code{\link{weighted_line}},
#' \code{\link{weighted_median_ratio}}
#' @export weighted_median_line
#' @importFrom stats complete.cases
#' @importFrom grDevices xy.coords
weighted_median_line <- function(x, y=NULL, w, type="slopes", na.rm=FALSE){

  if(missing(w)) stop('Argument w (weights) is missing, with no default.')

  if (inherits(x, "formula")) {
    dat <- xy.coords(x)
    x <- dat$x
    y <- dat$y
  }
  if (NCOL(x) > 1) return("x contains more than 1 explanatory variables.")

  # Robustification type
  stype <- pmatch(type, c("slopes"), nomatch=2)

  # Ensure case-wise completeness
  ok <- complete.cases(data.frame(x, y, w))
  n <- sum(ok)
  if (n < length(x)) {
    if (na.rm) {x <- x[ok]; y <- y[ok]; w <- w[ok]} else
      return("There are missing values in x, y or w.\n")
  }

  # univariate medians
  wmedx <- weighted_median(x, w)
  wmedy <- weighted_median(y, w)

  # slope (remove NA created due to division by 0)
  slope <- switch(stype,
                  weighted_median((y - wmedy) / (x - wmedx), w, na.rm=TRUE),
                  weighted_median((y - wmedy) * (x - wmedx), w, na.rm=TRUE) / weighted_median((x - wmedx)^2, w, na.rm=TRUE)
  )

  # residuals
  r <- y - slope * x

  # intercept
  intercept <- weighted_median(r, w)
  yhat <- intercept + slope * x

  # return
  structure(list(call = sys.call(), coefficients = c(intercept, slope),
                 residuals = y - yhat, fitted.values = yhat), class = "medline")

}



#' Weighted robust ratio based on median
#'
#' A weighted median of the ratios y/x determines the slope of a
#' regression through the origin.
#'
#' TBD
#'
#' @param x a numeric vector (explanatory variable)
#' @param y a numeric vector (response variable)
#' @param w a numeric vector of (optional) weights
#' @param na.rm a logical value indicating whether rows with \code{NA}
#' values should be stripped before the computation proceeds
#' @return a vector with two components: intercept and slope
#' @examples
#' x <- c(1,2,4,5)
#' y <- c(1,0,5,2)
#' weighted_median_ratio(y~x, w = rep(1, length(y)))
#' @seealso \code{\link[stats]{line}}, \code{\link{weighted_line}},
#' \code{\link{weighted_median_line}}
#' @export weighted_median_ratio
#' @importFrom stats complete.cases
#' @importFrom grDevices xy.coords
weighted_median_ratio <- function(x, y=NULL, w, na.rm=FALSE){

  if(missing(w)) stop('Argument w (weights) is missing, with no default.')

  if (inherits(x, "formula")) {
    dat <- xy.coords(x)
    y <- dat$y
    x <- dat$x
  }
  if (NCOL(x) > 1) return("x contains more than 1 explanatory variables.")

  # Ensure case-wise completeness
  ok <- complete.cases(data.frame(x, y, w))
  n <- sum(ok)
  if (n < length(x)) {
    if (na.rm) {x <- x[ok]; y <- y[ok]; w <- w[ok]} else
      return("There are missing values in x, y or w.\n")
  }

  # ratio
  ratio <- weighted_median(y / x, w, na.rm=TRUE)

  # fitted.varlues and residuals
  yhat <- ratio * x
  r <- y - yhat

  # return
  structure(list(call = sys.call(), coefficients = ratio,
                 residuals = y - yhat, fitted.values = yhat), class = "medline")

}

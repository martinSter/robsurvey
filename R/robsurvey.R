#' robsurvey: Robust survey statistics.
#'
#' The package robsurvey is a collection of functions for robust survey
#' statistics.
#'
#' @section robsurvey functions:
#' robust Horvitz-Thompson M-estimator of mean and total in \code{msvymean()}
#' and \code{msvytotal()}, robust trimmed Horvitz-Thompson estimator of mean
#' and total in \code{tsvymean()} and \code{tsvytotal()}, robust winsorized
#' Horvitz-Thompson estimator of mean and total in \code{wsvymean()} and
#' \code{wsvytotal()}, weighted median estimator in \code{weighted.median()},
#' weighted quantile estimator in \code{weighted.quantile()}, weighted median
#' absolute deviation in \code{weighted.mad()}, weighted mean and total
#' estimators in \code{weighted.mean()} and \code{weighted.total()}.
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



#' Bushfire scars.
#'
#' The bushfire data set was used by Campbell (1984, 1989) to locate bushfire scars.
#' The dataset contains satellite measurements on five frequency bands, corresponding
#' to each of 38 pixels.
#'
#' The data contains an outlying cluster of observations 33 to 38 a second outlier
#' cluster of observations 7 to 11 and a few more isolated outliers, namely observations
#' 12, 13, 31 and 32.
#'
#' For testing purposes weights are provided:
#' \code{bushfire.weights <- rep(c(1,2,5), length = nrow(bushfire))}
#'
#' @format A data frame with 38 rows and 5 variables.
#' @references Campbell, N. (1989) Bushfire Mapping using NOAA AVHRR Data. Technical Report.
#' Commonwealth Scientific and Industrial Research Organisation, North Ryde.
#' @examples
#' data(bushfire)
"bushfire"



#' Bushfire scars with missing data.
#'
#' The bushfire data set was used by Campbell (1984, 1989) to locate bushfire scars.
#' The dataset contains satellite measurements on five frequency bands, corresponding
#' to each of 38 pixels. However, this dataset contains missing values.
#'
#' The data contains an outlying cluster of observations 33 to 38 a second outlier
#' cluster of observations 7 to 11 and a few more isolated outliers, namely observations
#' 12, 13, 31 and 32.
#'
#' \code{bushfirem} is created from bushfire by setting a proportion of 0.2 of the values
#' to missing.
#'
#' For testing purposes weights are provided:
#' \code{bushfire.weights <- rep(c(1,2,5), length = nrow(bushfire))}
#'
#' @format A data frame with 38 rows and 5 variables.
#' @references Campbell, N. (1989) Bushfire Mapping using NOAA AVHRR Data. Technical Report.
#' Commonwealth Scientific and Industrial Research Organisation, North Ryde.
#' @examples
#' data(bushfirem)
"bushfirem"



#' Weights for Bushfire scars.
#'
#' The bushfire data set was used by Campbell (1984, 1989) to locate bushfire scars.
#' The dataset contains satellite measurements on five frequency bands, corresponding
#' to each of 38 pixels.
#'
#' For testing purposes, \code{bushfire.weights} provides artificial weights created
#' according to: \code{bushfire.weights <- rep(c(1,2,5), length = nrow(bushfire))}
#'
#' @format A vector of length 38.
#' @references Campbell, N. (1989) Bushfire Mapping using NOAA AVHRR Data. Technical Report.
#' Commonwealth Scientific and Industrial Research Organisation, North Ryde.
#' @examples
#' data(bushfire.weights)
"bushfire.weights"



#' Weighted lower sample median
#'
#' \code{weighted.median} computes the weighted lower sample median
#'
#' See \code{\link{weighted.quantile}} for further details.
#'
#' @param x a numeric vector whose weighted sample median is wanted
#' @param w a numeric vector of weights
#' @param na.rm a logical value indicating whether \code{NA} values should be
#'        stripped before the computation proceeds.
#' @return weighted sample median
#' @examples
#' x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
#' weighted.median(x, x)
#' @seealso \code{\link{weighted.quantile}}
#' @export weighted.median
#' @importFrom stats na.omit
#' @useDynLib robsurvey wquantile
weighted.median <- function(x, w, na.rm = FALSE){
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
   tmp <- .C("wquantile", x = as.double(dat[, 1]), w = as.double(dat[, 2]),
             probs = as.double(0.5), q = as.double(numeric(1)), n = as.integer(n))
   return(tmp$q)
}



#' Weighted lower sample quantiles
#'
#' \code{weighted.quantile} computes the weighted lower sample quantile
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
#' weighted.quantile(x, x, probs = c(0.25, 0.5, 0.75))
#' @references Cormen,T.H., Leiserson, C.E., Rivest, R.L., and Stein, C. (2009):
#'    Introduction to Algorithms, 3rd ed., Cambridge: MIT Press.
#' @seealso \code{\link{weighted.median}}
#' @export weighted.quantile
#' @importFrom stats na.omit
#' @useDynLib robsurvey wquantile
weighted.quantile <- function(x, w, probs, na.rm = FALSE){
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
#' \code{weighted.mad} computes weighted median absolute deviation from the
#' weighted median
#'
#' The weighted MAD is computed as the (normalized) weighted median of the
#' absolute deviation from the weighted median; the median is computed as the
#' weighted lower sample median (see \code{\link{weighted.median}}); the MAD
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
#' weighted.mad(x, x)
#' @seealso \code{\link{weighted.median}}
#' @export weighted.mad
#' @importFrom stats na.omit
#' @useDynLib robsurvey wmad
weighted.mad <- function(x, w, na.rm = FALSE, constant = 1.4826){
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
#' Tuning parameters for \code{\link{weighted.mean.huber}},
#' \code{\link{weighted.total.huber}}, \code{\link{msvymean}},
#' \code{\link{msvytotal}}.
#'
#' @param acc numeric tolerance, stoping rule in the iterative
#' updating scheme (default: \code{1e-5})
#' @param maxit maximum number of updating iterations
#' @param psi psi-function (\code{Huber} or \code{asymHuber})
#' @param ... additional arguments
#' @return List
#' @export rht.control
rht.control <- function(acc = 1e-5, maxit = 100, psi = "Huber", ...){
   if(!(psi %in% c("Huber", "asymHuber"))) stop("Function 'psi' must be
      either 'Huber' or 'asymHuber'\n")
   psi0 <- switch(psi,
      "Huber" = 0L,
      "asymHuber" = 1L)
   list(acc = unname(acc), maxit = unname(maxit), psi = unname(psi0))
}



#' @name wgtmeantotal
#' @aliases weighted.total
#' @aliases weighted.mean
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
#' weighted.total(x, x)
#' @export weighted.total
#' @importFrom stats na.omit
weighted.total <- function(x, w, na.rm = FALSE){
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
#' weighted.mean(x, x)
#' @export weighted.mean
#' @importFrom stats na.omit
weighted.mean <- function(x, w, na.rm = FALSE){
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
#' @aliases weighted.mean.huber
#' @aliases weighted.total.huber
#' @aliases msvymean
#' @aliases msvytotal
#'
#' @title Huber M-estimators of the weighted mean and weighted total
#'
#' @description Weighted Huber M-estimators of the mean and total are available in two forms:
#' \itemize{
#'    \item \strong{bare-bone} functions: \code{weighted.mean.huber} and
#'	 \code{weighted.total.huber},
#'    \item estimation \strong{methods}: \code{msvymean} and
#'	 \code{msvytotal} (incl. variance estimation
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
#'    \code{\link{weighted.mad}}).
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
#'    \code{\link{rht.control}} for details.
#'    }
#'    \item{\emph{Domain estimation}}{
#'    Estimates for domains can be obtained using the \link[survey]{svyby}
#'    wrapper in the \pkg{survey} package (see examples).
#'    }
#' }
#'
#' @section Utility functions:
#' For the methods \code{msvymean} and \code{msvytotal}, the following
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
#' name (\code{msvymean} or \code{msvytotal})
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
#' (see \code{\link{rht.control}})
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
#' @seealso \code{\link{tsvymean}}, \code{\link{tsvytotal}},
#' \code{\link{wsvymean}}, \code{\link{wsvytotal}},
#' \code{\link{weighted.mean.trimmed}}, \code{\link{weighted.total.trimmed}}
#' \code{\link{weighted.mean.winsorized}}, \code{\link{weighted.total.winsorized}}
#'
#' @rdname huberwgt
#' @export weighted.mean.huber
#' @importFrom stats na.omit
#' @useDynLib robsurvey rwlslm
weighted.mean.huber <- function(x, w, k, type = "rht", info = FALSE,
   na.rm = FALSE, ...){
   ctrl <- rht.control(...)
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
#' @export weighted.total.huber
weighted.total.huber <- function(x, w, k, type = "rht", info = FALSE,
   na.rm = FALSE, ...){
   res <- weighted.mean.huber(x, w, k, type, info, na.rm, ...)
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
#' msvymean(~api00, dstrat, k = 2)
#' # Domain estimates
#' svyby(~api00, by = ~stype, design = dstrat, msvymean, k = 1.34)
#' @export msvymean
#' @importFrom stats model.frame na.fail
#' @importFrom stats weights
#' @useDynLib robsurvey rwlslm
msvymean <- function(x, design, k, type = "rht", ...){
   ctrl <- rht.control(...)
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
	 stop("msvymean is not defined for object of class: ", class(x), "\n")
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
#' @export msvytotal
#' @importFrom stats weights
msvytotal <- function(x, design, k, ...){
   tmp <- msvymean(x, design, k, type = "rht", ...)
   tmp$characteristic <- "total"
   sumw <- sum(weights(design))
   tmp$estimate <- tmp$estimate * sumw
   tmp$variance <- tmp$variance * sumw^2
   tmp
}



#' @name trimwgt
#' @aliases weighted.mean.trimmed
#' @aliases weighted.total.trimmed
#' @aliases tsvymean
#' @aliases tsvytotal
#'
#' @title Weighted trimmed mean and trimmed total
#'
#' @description Weighted trimmed estimators of the mean and total are available in two forms:
#' \itemize{
#'    \item \strong{bare-bone} functions: \code{weighted.mean.trimmed} and
#'	 \code{weighted.total.trimmed},
#'    \item estimation \strong{methods}: \code{tsvymean} and
#'	 \code{tsvytotal} (incl. variance estimation
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
#' For the methods \code{tsvymean} and \code{tsvytotal}, the following
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
#' @param x numeric vector (\code{weighted.mean.trimmed} or \code{weighted.total.trimmed});
#' a formula object or variable name (\code{tsvymean} or \code{tsvytotal})
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
#' @seealso \code{\link{msvymean}}, \code{\link{msvytotal}},
#' \code{\link{wsvymean}}, \code{\link{wsvytotal}},
#' \code{\link{weighted.mean.huber}}, \code{\link{weighted.total.huber}},
#' \code{\link{weighted.mean.winsorized}}, \code{\link{weighted.total.winsorized}}
#'
#' @rdname trimwgt
#' @export weighted.mean.trimmed
#' @importFrom stats na.omit
#' @useDynLib robsurvey wmeantrimmed
weighted.mean.trimmed <- function(x, w, LB = 0.05, UB = 1 - LB, na.rm = FALSE){
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
#' @export weighted.total.trimmed
weighted.total.trimmed <- function(x, w, LB = 0.05, UB = 1 - LB, na.rm = FALSE){
   res <- weighted.mean.trimmed(x, w, LB, UB, na.rm)
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
#' tsvymean(~api00, dstrat, k = 2)
#' # Domain estimates
#' svyby(~api00, by = ~stype, design = dstrat, tsvymean, LB = 0.1)
#' @export tsvymean
#' @importFrom stats model.frame na.fail
#' @importFrom stats weights
tsvymean <- function(x, design, LB = 0.05, UB = 1 - LB, ...){
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
	 stop("tsvymean is not defined for object of class: ", class(x), "\n")
      }
   }
   w <- as.numeric(weights(design))
   est <- weighted.mean.trimmed(y, w, LB, UB)
   names(est) <- yname
   # compute influence function
   quant <- weighted.quantile(y, w, probs = c(LB, UB))
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
#' @export tsvytotal
#' @importFrom stats weights
tsvytotal <- function(x, design, LB = 0.05, UB = 1 - LB, ...){
   tmp <- tsvymean(x, design, LB, UB)
   tmp$characteristic <- "total"
   sumw <- sum(weights(design))
   tmp$estimate <- tmp$estimate * sumw
   tmp$variance <- tmp$variance * sumw^2
   tmp
}



#' @name winswgt
#' @aliases weighted.mean.winsorized
#' @aliases weighted.total.winsorized
#' @aliases wsvymean
#' @aliases wsvytotal
#'
#' @title Weighted winsorized mean and trimmed total
#'
#' @description Weighted winsorized estimators of the mean and total are available in two forms:
#' \itemize{
#'    \item \strong{bare-bone} functions: \code{weighted.mean.winsorized} and
#'	 \code{weighted.total.winsorized},
#'    \item estimation \strong{methods}: \code{wsvymean} and
#'	 \code{wsvytotal} (incl. variance estimation
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
#' For the methods \code{wsvymean} and \code{wsvytotal}, the following
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
#' @param x numeric vector (\code{weighted.mean.winsorized} or
#' \code{weighted.total.winsorized}); a formula object or variable name
#' (\code{wsvymean} or \code{wsvytotal})
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
#' @seealso \code{\link{msvymean}}, \code{\link{msvytotal}},
#' \code{\link{tsvymean}}, \code{\link{tsvytotal}},
#' \code{\link{weighted.mean.huber}}, \code{\link{weighted.total.huber}},
#' \code{\link{weighted.mean.trimmed}}, \code{\link{weighted.total.trimmed}}
#'
#' @rdname winswgt
#' @export weighted.mean.winsorized
#' @importFrom stats na.omit
#' @useDynLib robsurvey wmeanwinsorized
weighted.mean.winsorized <- function(x, w, LB = 0.05, UB = 1 - LB,
   na.rm = FALSE){
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
#' @export weighted.total.winsorized
weighted.total.winsorized <- function(x, w, LB = 0.05, UB = 1 - LB, na.rm = FALSE){
   res <- weighted.mean.winsorized(x, w, LB, UB, na.rm)
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
#' wsvymean(~api00, dstrat, k = 2)
#' # Domain estimates
#' svyby(~api00, by = ~stype, design = dstrat, wsvymean, LB = 0.1)
#' @export wsvymean
#' @importFrom stats model.frame na.fail
#' @importFrom stats weights
wsvymean <- function(x, design, LB = 0.05, UB = 1 - LB, ...){
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
	 stop("wsvymean is not defined for object of class: ", class(x), "\n")
      }
   }
   w <- as.numeric(weights(design))
   est <- weighted.mean.winsorized(y, w, LB, UB)
   names(est) <- yname
   # compute influence function
   # FIXME: variance => density estimate (now, it is trimmed mean variance)
   quant <- weighted.quantile(y, w, probs = c(LB, UB))
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
#' @export wsvytotal
#' @importFrom stats weights
wsvytotal <- function(x, design, LB = 0.05, UB = 1 - LB, ...){
   tmp <- wsvymean(x, design, LB, UB)
   tmp$characteristic <- "total"
   sumw <- sum(weights(design))
   tmp$estimate <- tmp$estimate * sumw
   tmp$variance <- tmp$variance * sumw^2
   tmp
}



#' Weighted robust line fitting
#'
#' \code{weighted.line} fits a robust line and allows weights.
#'
#' Uses different quantiles for splitting the sample than \code{line()}.
#' Is based on \code{weighted.median()}.
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
#' weighted.line(cars$speed, cars$dist)
#' weighted.line(cars$speed, cars$dist, w=rep(1:10, each=5))
#' @seealso \code{\link[stats]{line}}
#' @export weighted.line
#' @importFrom stats model.frame complete.cases
weighted.line <- function(x, y=NULL, w, na.rm=FALSE, iter = 1){

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
  if (missing(w)) w <- rep(1, length(x))
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
  wmedx <- c(weighted.median(x[groups==1], w[groups==1]),
             weighted.median(x[groups==3], w[groups==3]))

  # polishing (affects only the slope)
  slope <- 0; r <- y; j <- 0
  while (j <= iter - 1) {
    wmedr <- c(weighted.median(r[groups==1], w[groups==1]),
               weighted.median(r[groups==3], w[groups==3]))
    slope <- slope + (wmedr[2] - wmedr[1]) / (wmedx[2] - wmedx[1])
    r <- y - slope * x
    j <- j + 1
  }

  # intercept and predicted values
  intercept <- weighted.median(r, w)
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
#' Uses \code{weighted.median()}. The median of slopes (type="slopes")
#' uses \eqn{b1=M((y-M(y,w))/(x-M(x,w)), w)}. The median of crossproducts
#' by median of squares (type="products") uses
#' \eqn{b1=M((y-M(y,w))(x-M(x,w)), w )/ M((x-M(x,w )^2), w )}, where
#' \eqn{M(x, w)} is shorthand for the function \code{weighted.median(x, w)}.
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
#' weighted.line(y~x)
#' weighted.median.line(y~x)
#' weighted.median.line(y~x, type="prod")
#'
#' data(cars)
#' with(cars, weighted.median.line(dist ~ speed))
#' with(cars, weighted.median.line(dist ~ speed, type="prod"))
#'
#' # weighted
#' w <- c(rep(1,20), rep(2,20), rep(5, 10))
#' with(cars, weighted.median.line(dist ~ speed, w=w))
#' with(cars, weighted.median.line(dist ~ speed, w=w, type="prod"))
#'
#' # outlier in y
#' cars$dist[49] <- 360
#' with(cars, weighted.median.line(dist ~ speed))
#' with(cars, weighted.median.line(dist ~ speed, type="prod"))
#'
#' # outlier in x
#' data(cars)
#' cars$speed[49] <- 72
#' with(cars, weighted.median.line(dist ~ speed))
#' with(cars, weighted.median.line(dist ~ speed, type="prod"))
#' @seealso \code{\link[stats]{line}}, \code{\link{weighted.line}},
#' \code{\link{weighted.median.ratio}}
#' @export weighted.median.line
#' @importFrom stats complete.cases
#' @importFrom grDevices xy.coords
weighted.median.line <- function(x, y=NULL, w, type="slopes", na.rm=FALSE){

  if (inherits(x, "formula")) {
    dat <- xy.coords(x)
    x <- dat$x
    y <- dat$y
  }
  if (NCOL(x) > 1) return("x contains more than 1 explanatory variables.")

  if (missing(w)) w <- rep(1, length(x))

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
  wmedx <- weighted.median(x, w)
  wmedy <- weighted.median(y, w)

  # slope (remove NA created due to division by 0)
  slope <- switch(stype,
                  weighted.median((y - wmedy) / (x - wmedx), w, na.rm=TRUE),
                  weighted.median((y - wmedy) * (x - wmedx), w, na.rm=TRUE) / weighted.median((x - wmedx)^2, w, na.rm=TRUE)
  )

  # residuals
  r <- y - slope * x

  # intercept
  intercept <- weighted.median(r, w)
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
#' weighted.median.ratio(y~x)
#' @seealso \code{\link[stats]{line}}, \code{\link{weighted.line}},
#' \code{\link{weighted.median.line}}
#' @export weighted.median.ratio
#' @importFrom stats complete.cases
#' @importFrom grDevices xy.coords
weighted.median.ratio <- function(x, y=NULL, w, na.rm=FALSE){

  if (inherits(x, "formula")) {
    dat <- xy.coords(x)
    y <- dat$y
    x <- dat$x
  }
  if (NCOL(x) > 1) return("x contains more than 1 explanatory variables.")

  if (missing(w)) w <- rep(1, length(x))

  # Ensure case-wise completeness
  ok <- complete.cases(data.frame(x, y, w))
  n <- sum(ok)
  if (n < length(x)) {
    if (na.rm) {x <- x[ok]; y <- y[ok]; w <- w[ok]} else
      return("There are missing values in x, y or w.\n")
  }

  # ratio
  ratio <- weighted.median(y / x, w, na.rm=TRUE)

  # fitted.varlues and residuals
  yhat <- ratio * x
  r <- y - yhat

  # return
  structure(list(call = sys.call(), coefficients = ratio,
                 residuals = y - yhat, fitted.values = yhat), class = "medline")

}



#' Extraction of robustness weights (M-estimators)
#'
#' \code{robweights} retrieves the robustness weights from
#' an M-estimator of class \code{svystat.rob}.
#'
#' Extracts the robustness weights.
#'
#' @param object class of type \code{svystat.rob}
#' @return Vector of robustness weights
#' @export
robweights <- function(object){
   if(!inherits(object, "svystat.rob")){
      stop("robweights is not a valid method for this object!\n")
   }
   object$robust$robweights
}



#' @export
coef.svystat.rob <- function(object, ...){
   object$estimate
}

#' @export
print.svystat.rob <- function(x, digits = 3, ...){
   conv <- TRUE
   if(!is.null(x$optim)){
      conv <- x$optim$converged
   }
   if(conv){
      m <- cbind(x$estimate, sqrt(x$variance))
      colnames(m) <- c(x$characteristic, "SE")
      print(round(m, digits))
   }else{
      cat(paste0(x$call[[1]], ": failure of convergence in ", x$optim$niter,
	 " steps\n"))
      cat("(you may use the 'summary' method to see more details)\n")
   }
}

#' @export
residuals.svystat.rob <- function(object, ...){
   object$residuals
}

#' @export
summary.svystat.rob <- function(object, digits = 3, ...){
   cat(paste0("SUMMARY: ", object$estimator, " of the sample ",
      object$characteristic, "\n"))
   cat("\n")
   est <- cbind(round(object$estimate, digits), round(sqrt(object$variance),
      digits), length(object$residuals))
   colnames(est) <- c(object$characteristic, "SE", "n")
   print(est)
   cat("\n")
   if(!is.null(object$optim)){
      cat("ROBUSTNESS PROPERTIES\n")
      cat(paste0("  Psi-function: ", object$robust$psifunction, " with k = ",
	 object$robust$k, "\n"))
      cat(paste0("  mean of robustness weights: ", round(mean(object$robust$robweights),
	 digits), "\n"))
      cat("\n")
      cat("ALGORITHM PERFORMANCE \n")
      if (object$optim$converged){
	 cat(paste0("  IRLS converged in ", object$optim$niter, " iterations \n"))
	 cat(paste0("  with residual scale (MAD): ", round(object$robust$scale,
	    digits), "\n"))
      }else{
	  cat(paste0("  FAILURE of convergence in ", object$optim$niter, " iterations \n"))
	 cat(paste0("  with residual scale (MAD): ", round(object$robust$scale,
	    digits), "\n"))
      }
      cat("\n")
   }
   cat("SAMPLING DESIGN\n")
   print(object$design)
}

#' @export
vcov.svystat.rob <- function(object, ...){
   v <- as.matrix(object$variance)
   rownames(v) <- names(object$estimate)
   colnames(v) <- "Variance"
   v
}

#' @export
fitted.svystat.rob <- function(object, ...){
   object$model$y - object$residuals
}

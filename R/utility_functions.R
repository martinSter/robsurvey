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

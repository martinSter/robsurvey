summary.svystat.rob <-
function(object, digits = 3, ...){
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

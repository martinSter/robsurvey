print.svystat.rob <-
function(x, digits = 3, ...){
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

svymean.huber <-
function(x, design, k, type = "rht", ...){
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
	 stop("svymean is not defined for object of class: ", class(x), "\n")
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
   v <- as.numeric(attr(svymean(~zz, design), "var"))
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

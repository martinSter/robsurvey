svymean.trimmed <-
function(x, design, LB = 0.05, UB = 1 - LB, ...){
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
   v <- as.numeric(attr(svymean(~zz, design), "var"))
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

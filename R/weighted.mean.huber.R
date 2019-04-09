weighted.mean.huber <-
function(x, w, k, type = "rht", info = FALSE, 
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

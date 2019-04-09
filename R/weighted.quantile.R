weighted.quantile <-
function(x, w, probs, na.rm = FALSE){
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

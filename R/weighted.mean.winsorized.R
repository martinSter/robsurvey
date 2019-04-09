weighted.mean.winsorized <-
function(x, w, LB = 0.05, UB = 1 - LB, 
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

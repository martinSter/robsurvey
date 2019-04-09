weighted.mad <-
function(x, w, na.rm = FALSE, constant = 1.4826){
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

weighted.mean <-
function(x, w, na.rm = FALSE){
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

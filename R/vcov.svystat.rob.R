vcov.svystat.rob <-
function(object, ...){
   v <- as.matrix(object$variance)
   rownames(v) <- names(object$estimate)
   colnames(v) <- "Variance"
   v 
}

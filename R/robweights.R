robweights <-
function(object){
   if(!inherits(object, "svystat.rob")){
      stop("robweights is not a valid method for this object!\n")
   }
   object$robust$robweights
}

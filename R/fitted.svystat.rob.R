fitted.svystat.rob <-
function(object, ...){
   object$model$y - object$residuals
}

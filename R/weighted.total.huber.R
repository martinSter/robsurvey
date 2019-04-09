weighted.total.huber <-
function(x, w, k, type = "rht", info = FALSE, 
   na.rm = FALSE, ...){
   res <- weighted.mean.huber(x, w, k, type, info, na.rm, ...)
   if(length(res) == 1){
      res <- res * sum(w) 
   }else{
      res$characteristic <- "total"
      res$estimate <- res$estimate * sum(w)
   }
   return(res)
}

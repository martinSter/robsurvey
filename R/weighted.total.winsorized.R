weighted.total.winsorized <-
function(x, w, LB = 0.05, UB = 1 - LB, na.rm = FALSE){
   res <- weighted.mean.winsorized(x, w, LB, UB, na.rm)
   if(length(res) == 1){
      res <- res * sum(w) 
   }else{
      res$characteristic <- "total"
      res$estimate <- res$estimate * sum(w)
   }
   return(res)
}

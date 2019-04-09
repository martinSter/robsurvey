svytotal.huber <-
function(x, design, k, ...){
   tmp <- svymean.huber(x, design, k, type = "rht", ...)
   tmp$characteristic <- "total"
   sumw <- sum(weights(design))
   tmp$estimate <- tmp$estimate * sumw 
   tmp$variance <- tmp$variance * sumw^2 
   tmp
}

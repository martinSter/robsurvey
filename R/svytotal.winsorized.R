svytotal.winsorized <-
function(x, design, LB = 0.05, UB = 1 - LB, ...){
   tmp <- svymean.winsorized(x, design, LB, UB)
   tmp$characteristic <- "total"
   sumw <- sum(weights(design))
   tmp$estimate <- tmp$estimate * sumw 
   tmp$variance <- tmp$variance * sumw^2 
   tmp
}

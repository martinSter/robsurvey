rht.control <-
function(acc = 1e-5, maxit = 100, psi = "Huber", ...){
   if(!(psi %in% c("Huber", "asymHuber"))) stop("Function 'psi' must be 
      either 'Huber' or 'asymHuber'\n")
   psi0 <- switch(psi, 
      "Huber" = 0L,
      "asymHuber" = 1L)
   list(acc = unname(acc), maxit = unname(maxit), psi = unname(psi0))
}

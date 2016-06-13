PoissonVariantCall <- function(mod.comb, name, variant, error.rate) {
  
  n <- seq(0, mod.comb)
  lambda <-  error.rate * mod.comb
  aux <- exp(-lambda) * lambda^n / factorial(n)
  pval <- 1 - cumsum(aux)
  p.val <- p.adjust(pval, method = "bonferroni")
  cutoff <- which(p.val < 0.001)[1]
  aux.ret <- table(unlist(variant[names(variant) == name]))
  if (length(aux.ret) > 0){
    ret <- aux.ret[aux.ret >= cutoff]
  } else {
    ret <- integer(0)
  }
  return(ret)
  
}
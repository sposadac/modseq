MemTrace <- function() {
  bit <- 8L * .Machine$sizeof.pointer
  #if (!(bit == 32L || bit == 64L)) {
  #  stop("Unknown architecture", call. = FALSE)
  #}
  node_size <- if (bit == 32L) 28L else 56L
  usage <- gc() #Vcells: memory used by vectors. Ncells: memory used by everything else. 
  sum(usage[, 1] * c(node_size, 8)) / (1024 ^ 2)
}
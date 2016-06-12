ParseInfoResultsList <- function(res.list,  mod.comb = character(0), 
                                 mod.comb.len = NULL, num.cores = numeric(0)) {
  
  ## Function arguments
  #  mod.comb     Variable of type DNAStringSet containing sequences of the module
  #               combinations. Alternatively path to fasta file.
  
  if (length(num.cores) == 0) {
    num.cores <- detectCores()
  }
  
  if (is.character(mod.comb)) {
    
    if (file.exits(mod.comb)) {
      cat("Reading reference sequences ... \n")
      mod.comb <- readDNAStringSet(mod.comb)
    } else {
      stop("File \"", mod.comb, "\" not found.\n")
    }
    
  }
  if (is.null(mod.comb.len)) {
    mod.comb.len <- length(mod.comb)
  }
  
  
  res.list.lengths <- 
    unlist(mclapply(res.list, length, mc.cores = num.cores), 
           use.names = FALSE)
  
  ## Getting the IDs of the module combinations. 
  # An IDs is formed as the concatenation of IDs of modular variants
  if (sum(names(res.list) == mod.comb) != mod.comb.len) {
    
    ## Only needed if the order is not preserved
    mod.comb.seq2names <- names(mod.comb)
    names(mod.comb.seq2names) <- mod.comb
    res.list.ids <- 
      unlist(mcmapply(function (x,n) rep(x,n), 
                      x = mod.comb.seq2names[names(res.list)], 
                      n = res.list.lengths, USE.NAMES = FALSE,
                      mc.cores = num.cores), use.names = FALSE)
    
  } else {
    res.list.ids <- 
      unlist(mcmapply(function (x,n) rep(x,n), x = names(mod.comb), 
                      n = res.list.lengths, USE.NAMES = FALSE,
                      mc.cores = num.cores), use.names = FALSE)
  }
  
  return(list(res.list.lengths, res.list.ids))
}
ParseInfoResultsList <- function(res.list,  mod.comb = character(0), seq2id = FALSE,
                                 num.cores = numeric(0)) {
  
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
  
  
  res.list.lengths <- 
    unlist(mclapply(res.list, length, mc.cores = num.cores), 
           use.names = FALSE)
  
  ## Getting the IDs of the module combinations. 
  # An IDs is formed as the concatenation of IDs of modular variants
  if (any(names(res.list) != mod.comb)) {
    
    ## Only needed if the order is not preserved
    mod.comb.seq2names <- names(mod.comb)
    names(mod.comb.seq2names) <- mod.comb
    res.list.ids <- as.character(mod.comb.seq2names[names(res.list)])
    
  } else {
    
    res.list.ids <- names(mod.comb)
    
  }
  
  return(list(res.list.lengths, res.list.ids))
}
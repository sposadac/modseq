GetVariantsNames <- function(patterns, modules.tot = NULL, modseq.dir = NULL, 
                             num.cores = numeric(0)) {
  
  ## Function arguments
  # patterns      Data frame containing the sequences of the various variants
  #               per module. Columns correspond to different modules and rows 
  #               correspond to variants identifiers. Alternatively, it can 
  #               represent the path to the csv file.
  # modules.tot   (optional) total number of modules.
  # num.cores     (optional) number of cores available for performing parallel 
  #               tasks.
  
  ## Whenever modseq path is not specified, assumed to be the current working
  #  directory
  if (is.null(modseq.dir)) {
    modseq.dir <- getwd()
    warning("Object \'modseq.dir\' not found, set to: \"", modseq.dir, "\".")
  }
  
  if (!is.data.frame(patterns)) {
    source(file.path(modseq.dir, "R/functions/LoadModuleTable.R"))
    if (is.character(patterns)) {
      patterns <- LoadModuleTable(in.modulesDir = patterns) 
    } else {
      stop("Invalid data type for \'patterns\', \"", class(patterns), "\".\n")
    }
  }
  
  if (is.null(modules.tot)) {
    modules.tot <- ncol(patterns) 
  }
  
  if (length(num.cores) == 0) {
    num.cores <- detectCores()
  }
  
  
  ret <- unlist(mclapply(seq_len(modules.tot), 
                  function(x) row.names(na.omit(patterns[x])), 
                  mc.cores = num.cores))
  
  return(ret)
}
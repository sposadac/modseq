LoadModuleTable <- function(in.modulesDir, modules.filename=NULL, list=FALSE) {
  
  ## Function arguments
  # in.modulesDir     Path to table of modules. It can be either the path to  
  #                   table of modules or the absolute path, including the name 
  #                   and file extension of the input file.
  
  if (is.null(modules.filename)) {
    in.file <- file.path(in.modulesDir)
  } else{
    in.file <- 
      file.path(in.modulesDir, paste(modules.filename, ".csv", sep = ""))
  }
  
  if (file.exists(in.file)) {
    
    patterns <- 
      data.frame(
        read.csv(file = in.file, header = TRUE, sep = ",", dec = "."), 
        stringsAsFactors = FALSE
      )
    cat("Table of modules used: \"", in.file, "\". \n", sep = "")
    
    if (list) {
      mod.number <- sapply(patterns, function(x) length(na.omit(x)))
      
      patterns.list <- as.character(na.omit(as.vector(as.matrix(patterns))))
      if (length(unique(row.names(patterns))) == sum(mod.number)) {
        ## Apply when each modular variant has a unique name
        names(patterns.list) <- row.names(patterns)
        
      } else {
        ## Create unique names for each modular variant (ignoring common  
        #  feature among variants or different modules).
        names(patterns.list) <- 
          tolower(names(unlist(lapply(mod.number, function(x) seq(1:x)))))
        
      }
      
      return(list(patterns, patterns.list, mod.number))
      
    } else {
      return(patterns)
      
    }
    
  } else {
    
    stop("File \"", in.file, "\" not found. \n")
    
  }
  
}


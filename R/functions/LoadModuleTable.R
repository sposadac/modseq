LoadModuleTable <- function(in.modulesDir, modules.filename = NULL) {
  
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
    
    return(patterns)
    
  } else {
    
    stop("File \"", in.file, "\" not found. \n")
    
  }
  
}
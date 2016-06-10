LoadModuleTable <- function(in.modulesDir, modules.filename = NULL) {
  
  if (is.null(modules.filename)) {
    in.file <- file.path(in.modulesDir)
  } else{
    in.file <- 
      file.path(in.modulesDir, paste(modules.filename, ".csv", sep = ""))
  }
  
  cat("Module-table used: \"", in.file, "\". \n", sep = "")
  patterns <- 
    data.frame(
      read.csv(file = in.file, header = TRUE, sep = ",", dec = "."), 
      stringsAsFactors = FALSE
    )
  
  return(patterns)
  
}
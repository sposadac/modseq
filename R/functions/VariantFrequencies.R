VariantFrequencies <- function(data, patterns, num.hits, out.filename, 
                               modseq.dir=NULL, out.dir=NULL, 
                               num.cores=numeric(0)) {
  
  ## Function arguments
  # patterns      Complete path to the csv file (table of modules).
  # out.dir       Path to output files.
  # modseq.dir    (optional) modseq's head directory, by default current working
  #               directory.
  # num.cores     (optional) number of cores available for performing parallel 
  #               tasks.
  
  ## Whenever modseq path or output directory are not specified, assumed to be 
  #  the current working directory
  if (is.null(modseq.dir)) {
    modseq.dir <- getwd()
    warning("Object \'modseq.dir\' not found, set to: \"", modseq.dir, "\".")
  } 
  if (is.null(out.dir)) {
    out.dir <- getwd()
  }
  
  if (!existsFunction("GetVariantsNames")) {
    source(file.path(modseq.dir, "R/functions/GetVariantsNames.R"))
  }
  if (!existsFunction("PlotVariantFrequencies")) {
    source(file.path(modseq.dir, "R/functions/PlotVariantFrequencies.R"))
  }
  
  if (length(num.cores) == 0) {
    num.cores <- detectCores()
  }
  
  ## Loading module table
  if (!existsFunction("LoadModuleTable")) {
    source(file.path(modseq.dir, "R/functions/LoadModuleTable.R"))
  }
  
  if (is.character(patterns)) {
    retList <- LoadModuleTable(in.modulesDir = patterns, list = TRUE) 
    patterns <- retList[[1]]
    patterns.list <- retList[[2]]
    mod.number <- retList[[3]]
    
  } else {
    stop("Invalid data type for \'patterns\', \"", class(patterns), "\".\n")
  }
  mod.tot <- ncol(patterns)  
  
  ## Frequencies of modular variants
  aux.modDist <- table(
    factor(unlist(strsplit(data, ":")), levels = names(patterns.list))
  )
  
  mod.id2names <- 
    GetVariantsNames(patterns = patterns, modules.tot = mod.tot, 
                     modseq.dir = modseq.dir, num.cores = num.cores)
  names(mod.id2names) <- names(patterns.list)
  
  df.modDist <- data.frame(
    "module" =
      factor(rep(names(patterns), mod.number), levels = names(patterns)),
    "variant" = mod.id2names[names(aux.modDist)], 
    "counts"  = as.numeric(aux.modDist), row.names = NULL
  )
  
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, "_FreqVariantsPerModule_table.csv", sep = ""))
  cat("Exporting distribution of modular variants: \"", out.file, "\". \n",
      sep = "")
  write.csv(df.modDist, file = out.file)
  
  df.modDist$variantNumber <- unlist(
    mclapply(mod.number, seq_len, mc.cores = num.cores), use.names = FALSE)
  df.modDist$counts <- df.modDist$counts * 100 / num.hits 
  
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, "_FreqVariantsPerModule_graph.pdf", sep = ""))
  cat("Plotting distribution of modular variants in the last search round: \"",
      out.file, "\". \n", sep = "")
  PlotVariantFrequencies(data = df.modDist[df.modDist$counts != 100, ], 
                         x.var = "variantNumber", y.var = "counts", 
                         z.var = "variant", facet = "module", 
                         out.file = out.file)
  
  return(aux.modDist)
  
}
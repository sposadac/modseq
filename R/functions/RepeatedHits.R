# **NOTE**: Read can be mapped to multiple module combinations (even when 
#           gls.mma = 0) if paired-end reads were almost stitched, rather than
#           merged. This can happen when a very short overlap is found. Then,   
#           if you have both copies which differ in few positions, such read  
#           can be mapped to different module combinations.
RepeatedHits <- function(res.list, data = NULL, patterns = NULL, num.reads, 
                         in.modDir, mod.filename, res.listName, out.dir, 
                         modseq.dir = NULL, num.cores = numeric(0)) {
  
  ## Function arguments
  # res.list      List containing IDs of reads assigned to each module 
  #               combination.
  # data          List of the read IDs. Each ID is encountered as many times as
  #               the corresponding read has been assigned to a module 
  #               combination.
  # patterns      Data frame containing the sequences of the various variants 
  #               per module. Columns correspond to different modules and rows 
  #               correspond to variant identifiers. 
  # in.modDir     Path to file containing sequences of modular variants.
  # mod.filename  Name of file containing sequences of modular variants 
  #               (without file extension).
  # res.listName
  # out.dir       Path to output directory
  # modseq.dir    (optional) modseq's head directory, by default current working
  #               directory.
  # num.cores     (optional) number of cores available for performing parallel 
  #               tasks.
  
  ## Whenever modseq path is not specified, assumed to be the current working
  #  directory
  if (is.null(modseq.dir)) {
    modseq.dir <- getwd()
    warning("Object \'modseq.dir\' not found, set to: \"", modseq.dir, "\".")
  } 
  
  source(file.path(modseq.dir, "R/functions/PlotRepeatedHits.R"))
  source(file.path(modseq.dir, "R/functions/ParseInfoResultsList.R"))
  
  if (length(num.cores) == 0) {
    num.cores <- detectCores()
  }
  
  if (!is.data.frame(patterns)) {
    
    if (!existsFunction("LoadModuleTable")) {
      source(file.path(modseq.dir, "R/functions/LoadModuleTable.R"))
    }
  
    patterns <- 
      LoadModuleTable(in.modulesDir = in.modDir, 
                      modules.filename = mod.filename, list = FALSE) 
  }
  
  if (is.null(data)) {
    data <- unlist(res.list, use.names = FALSE)
  } 
  rep.tab <- table(data)
  
  out.file <- file.path(
    out.dir, paste(res.listName, "_repeated_hits_graph.pdf", sep = ""))
  cat("Printing plot of repeated hits in: \"", out.file, "\".\n", sep = "")
  PlotRepeatedHits(data = rep.tab, out.file = out.file)
  
  ## Generating module combinations and loading fasta file 
  mod.file <- 
    file.path(in.modDir, paste(mod.filename, "_modComb.fasta", sep = ""))
  
  if (!file.exists(mod.file)) {
    
    if (!existsFunction("ModuleCombinationsGen")) {
      source(file.path(modseq.dir, "R/functions/ModuleCombinationsGen.R"))
    }
    
    ModuleCombinationsGen(modules.filename = mod.filename, pattern = patterns, 
                          in.modDir = in.modDir, modseq.dir = modseq.dir, 
                          out.file = mod.file, num.cores = num.cores)
  }
  mod.comb <- readDNAStringSet(mod.file)
  mod.comb.len <- length(mod.comb)
  
  ## Number of hits per module combination
  retList <- 
    ParseInfoResultsList(res.list = res.list, mod.comb = mod.comb, 
                         mod.comb.len = mod.comb.len, num.cores = num.cores)
  res.list.lengths <- retList[[1]]
  res.list.ids <- retList[[2]]

  rep.hits <- names(rep.tab[rep.tab > 1])
  aux.rep.hits <- 
    unlist(mclapply(rep.hits, function(x) res.list.ids[which(data == x)],
                    mc.cores= num.cores))
  df.rep.hits <- data.frame("readID" = rep(rep.hits, rep.tab[rep.tab > 1]),
                            "refID"  = aux.rep.hits,
                            "refSEQ" = as.character(mod.comb[aux.rep.hits]))
  
  # Exporting reads which have been assigned to multiple module combinations
  # as csv file
  out.file <- 
    file.path(out.dir, 
              paste(res.listName, "_repeated_hits_table.csv", sep = ""))
  cat("Exporting repeated hits as table: \"", out.file, ".\" \n", sep = "")
  write.csv(df.rep.hits, out.file, row.names = TRUE)
  print(
    paste(rep.hits, "- read was counted", rep.tab[rep.tab > 1], "times"))
  
  # Exporting list cotaining reads which have been assigned to multiple module
  # combinations as a R object 
  out.file <- 
    file.path(out.dir, paste(res.listName, "_repeated_hits.rda", sep = ""))
  cat("Exporting list containing repeated hits: \"", out.file, "\".\n",
      sep = "")
  save(df.rep.hits, file = out.file)
  
  cat(sum(rep.tab == 1), "reads mapped exactly 1 time. \n")
  cat(sum(rep.tab > 1), " reads (", 
      round(sum(rep.tab > 1) * 100 / num.reads, digit = 2), 
      "%) mapped > 1 times. \n", sep = "")
  cat("Corrected for multiple hits: Found ", length(rep.tab), 
      " module combinations in ", num.reads, " reads (", 
      round(length(rep.tab) / num.reads * 100, digits = 2), "%).\n",
      sep = "")
  
  return(list(mod.comb, res.list.lengths, res.list.ids))
  
}
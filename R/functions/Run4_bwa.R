Run4_bwa <- function(patterns, mod.comb, res.sam.realn, bwa.dupl=TRUE,
                     coverage.left=0, coverage.right=0, mapQ.thold=8, 
                     editDist.thold=8, num.reads=NULL, out.filename, 
                     modseq.dir=NULL, out.dir=NULL, num.cores=numeric(0)) {
  
  ## Function arguments
  # patterns      Complete path to the csv file (table of modules).
  # modseq.dir    (optional) modseq's head directory, by default current working
  #               directory.
  
  ## Whenever modseq path or output directory are not specified, assumed to be 
  #  the current working directory
  if (is.null(modseq.dir)) {
    modseq.dir <- getwd()
    warning("Object \'modseq.dir\' not found, set to: \"", modseq.dir, "\".")
  } 
  if (is.null(out.dir)) {
    out.dir <- getwd()
  }
  
  if (length(num.cores) == 0) {
    num.cores <- detectCores()
  }
  
  if (!existsFunction("AlignmentsFiltering")) {
    source(file.path(modseq.dir, "R/functions/AlignmentsFiltering.R"))
  }
  if (!existsFunction("VariantFrequencies")) {
    source(file.path(modseq.dir, "R/functions/VariantFrequencies.R"))
  }
  if (!existsFunction("PlotModuleCombinationsDistribution")) {
    source(
      file.path(modseq.dir, "R/functions/PlotModuleCombinationsDistribution.R"))
  }
  
  mod.comb.len <- length(mod.comb)
  
  
  retList <- AlignmentsFiltering(
    mod.comb = mod.comb, res.sam.realn = res.sam.realn, bwa.dupl = bwa.dupl,
    coverage.left = coverage.left, coverage.right = coverage.right, 
    mapQ.thold = mapQ.thold, editDist.thold = editDist.thold, 
    num.reads =  num.reads, modseq.dir = modseq.dir, num.cores = num.cores
    )
  res.sam.filt <- retList[[1]]
  editDistance <- retList[[2]]
  ind.edit <- retList[[3]]
  
  ## Exporting object containing results of the alignment filtering 
  out.file <- 
    file.path(out.dir, paste("alignmentFilter_", out.filename, ".rda", sep = ""))
  cat("Exporting object containing results of the alignment filtering : \"", 
      out.file, "\".\n", sep = "")
  save(res.sam.filt, file = out.file)
  
  ###### 1: Distribution of modular variants
  aux.modDist <- VariantFrequencies(
    data = res.sam.filt[["refID"]][ind.edit], patterns = patterns, 
    num.hits = length(ind.edit), out.filename = out.filename,
    modseq.dir = modseq.dir, out.dir = out.dir, num.cores = num.cores
  )

  ##### 2: Distribution of module combinations
  mod.dstr <- table(res.sam.filt[["refID"]][ind.edit])
  mod.dstr <- sort(mod.dstr, decreasing = TRUE)
  length(mod.dstr) <- mod.comb.len
  df.mod.dstr <- data.frame("x" = seq_len(mod.comb.len) * 100 / mod.comb.len, 
                            "y" = mod.dstr)
  
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, "_ModComb_distribution_graph.pdf", sep = ""))
  cat("Plotting distribution of module combinations: \"", out.file, "\". \n", 
      sep = "")
  
  PlotModuleCombinationsDistribution(data = df.mod.dstr, x.var = "x", 
                                     y.var = "y", out.file = out.file)
  
  return(list(res.sam.filt, editDistance, aux.modDist))
  
}
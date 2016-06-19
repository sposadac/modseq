library(Rcpp)

Run5_bwa <- function(mod.comb, patterns, res.sam.filt=NULL, aux.modDist=NULL, 
                     bwa.dupl=TRUE, coverage.left=0, coverage.right=0, 
                     mapQ.thold=8, editDist.thold=8, mismatch.filter=TRUE, 
                     out.varFiles=FALSE, num.reads=NULL, out.filename, 
                     modseq.dir=NULL, out.dir=NULL, num.cores=numeric(0)) {
  
  ## Function arguments
  # mod.comb      Variable of type DNAStringSet containing sequences of the 
  #               module.
  # aux.modDist   Table containing the number of hits per modular variant.
  # res.sam.filt  List containing the results of the alignment filtering.
  # 
  
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
  
  if (is.character(mod.comb)) {
    
    if (file.exits(mod.comb)) {
      cat("Reading reference sequences ... \n")
      mod.comb <- readDNAStringSet(mod.comb)
    } else {
      stop("File \"", mod.comb, "\" not found.\n")
    }
    
  }
  
  if (!existsFunction("LoadModuleTable")) {
    source(file.path(modseq.dir, "R/functions/LoadModuleTable.R"))
  }
  if (!existsFunction("SeqTrimming")) {
    source(file.path(modseq.dir, "R/functions/SeqTrimming.R"))
  }
  if (!existsFunction("PlotMismatchesPerVariant")) {
    source(file.path(modseq.dir, "R/functions/PlotMismatchesPerVariant.R"))
  }
  if (!existsFunction("PlotMismatches")) {
    source(file.path(modseq.dir, "R/functions/PlotMismatches.R"))
  }
  if (!existsFunction("PlotIndels")) {
    source(file.path(modseq.dir, "R/functions/PlotIndels.R"))
  }
  if (!existsFunction("PoissonVariantCall")) {
    source(file.path(modseq.dir, "R/functions/PoissonVariantCall.R"))
  }
  
  mod.comb.len <- length(mod.comb)
  
  ## Loading module table
  if (is.character(patterns)) {
    retList <- LoadModuleTable(in.modulesDir = patterns, list = TRUE) 
    patterns <- retList[[1]]
    patterns.list <- retList[[2]]
    
  } else {
    stop("Invalid data type for \'patterns\', \"", class(patterns), "\".\n")
  }
  
  if (is.null(res.sam.filt) || is.null(aux.modDist)) {
    
    if (!existsFunction("AlignmentsFiltering")) {
      source(file.path(modseq.dir, "R/functions/AlignmentsFiltering.R"))
    }
    
    res.sam.realn <- 
      file.path(out.dir, paste(out.filename, "_realigned.sam", sep = ""))
    
    retList <- AlignmentsFiltering(
      mod.comb = mod.comb, res.sam.realn = res.sam.realn, bwa.dupl = bwa.dupl,
      coverage.left = coverage.left, coverage.right = coverage.right, 
      mapQ.thold = mapQ.thold, editDist.thold = editDist.thold, 
      num.reads =  num.reads, modseq.dir = modseq.dir, num.cores = num.cores
    )
    res.sam.filt <- retList[[1]]
    editDistance <- retList[[2]]
    ind.edit <- retList[[3]]
    
    aux.modDist <- table(
      factor(unlist(strsplit(res.sam.filt[["refID"]][ind.edit], split = ":")), 
             levels = names(patterns.list))
      )
    
  }
  
  mod.len <- nchar(patterns.list)
  ind.edit <- which(editDistance <= editDist.thold & editDistance > 0)
  
  ### 1. "Soft Clipping" reads and quality strings 
  cigar.ops <- 
    explodeCigarOps(res.sam.filt[["CIGAR"]][ind.edit]) 
  cigar.ops.len <- 
    explodeCigarOpLengths(res.sam.filt[["CIGAR"]][ind.edit])
  
  reads.trimmed <- 
    DNAStringSet(mcmapply(
      SeqTrimming, cigarOps = cigar.ops, cigarOpLengths = cigar.ops.len, 
      seq = res.sam.filt[["readSeq"]][ind.edit], mc.cores = num.cores))
  
  if (out.varFiles) {
    ## Exporting object containing trimmed reads as an R object
    out.file <- 
      file.path(out.dir, 
                paste(out.filename, "_trimmedReads.rda", sep = ""))
    cat("Exporting object containing trimmed reads : \"", 
        out.file, "\".\n", sep = "")
    save(reads.trimmed, file = out.file)
  }
  
  qual.trimmed <- 
    BStringSet(mcmapply(
      SeqTrimming, cigarOps = cigar.ops, cigarOpLengths = cigar.ops.len, 
      seq = res.sam.filt[["readQuality"]][ind.edit], mc.cores = num.cores))
  
  if (out.varFiles) {
    ## Exporting object containing trimmed quality-strings as R object
    out.file <- 
      file.path(out.dir, 
                paste(out.filename, "_trimmedQualityScores.rda", sep = ""))
    cat("Exporting object containing trimmed quality-strings : \"", 
        out.file, "\".\n", sep = "")
    save(qual.trimmed, file = out.file)
  }
  
  ## **TODO**: MD field is not recomputed for realigned reads, write module
  #  NULL as output for OC (stands for original cigar) - tag
  md.ops <- mclapply(
    strsplit(res.sam.filt[["mismatchPOS"]][ind.edit], split = ":"),
    function(x) strsplit(x[3], split = "[0-9]+")[[1]], mc.cores = num.cores)
  
  md.ops.len <-  mclapply(
    strsplit(res.sam.filt[["mismatchPOS"]][ind.edit], split = ":"),
    function(x) as.numeric(strsplit(x[3], split = "[ACTGN^]+")[[1]]), 
    mc.cores = num.cores)
  
  mod.len.cumsum <- mclapply(
    res.sam.filt[["refID"]][ind.edit], 
    function(x) cumsum(
      c(1, as.numeric(head(mod.len[strsplit(x, split = ":")[[1]]], n = -1)))), 
    mc.cores = num.cores)
  
  names.mod.len.cumsum <- mclapply(
    res.sam.filt[["refID"]][ind.edit],
    function(x) names(mod.len[strsplit(x, split = ":")[[1]]]), 
    mc.cores = num.cores)
  
  leftclipped <- as.integer(res.sam.filt[["refPOS"]][ind.edit] - 1)
  
  
  ### 2. Identification of mismatches (pileup approach)
  if (mismatch.filter) { 
    sourceCpp(
      file.path(modseq.dir, 
                "src/mdSubstitutionAlongReferenceSpaceQualFilter.cpp"))
    
    substitution.modPos <- mcmapply(
      FUN = mdSubstitutionAlongReferenceSpaceQualFilter, 
      md_ops_lengths = md.ops.len, md_ops = md.ops, cigar_ops = cigar.ops, 
      cigar_ops_lengths =  cigar.ops.len, read = as.character(reads.trimmed), 
      qual = as.character(qual.trimmed), mod_len_cumsum = mod.len.cumsum, 
      names_mod_len_cumsum = names.mod.len.cumsum, offset = leftclipped, 
      qthold = mismatch.qthold, mc.cores = num.cores
    )
    
  } else {
    
    sourceCpp(
      file.path(modseq.dir, "src/mdSubstitutionAlongReferenceSpace.cpp"))
    
    substitution.modPos <- mcmapply(
      FUN = mdSubstitutionAlongReferenceSpace, md_ops_lengths = md.ops.len, 
      md_ops = md.ops, cigar_ops = cigar.ops, 
      cigar_ops_lengths =  cigar.ops.len, read = as.character(reads.trimmed), 
      mod_len_cumsum = mod.len.cumsum, 
      names_mod_len_cumsum = names.mod.len.cumsum, offset = leftclipped, 
      mc.cores = num.cores
    )
    
    qualSubstitution.modPos <- mcmapply(
      FUN = mdSubstitutionAlongReferenceSpace, md_ops_lengths = md.ops.len,
      md_ops = md.ops, cigar_ops = cigar.ops,
      cigar_ops_lengths =  cigar.ops.len, read = as.character(qual.trimmed), 
      mod_len_cumsum = mod.len.cumsum, 
      names_mod_len_cumsum = names.mod.len.cumsum, offset = leftclipped,
      mc.cores = num.cores
    )
    
  }
  names(substitution.modPos) <- res.sam.filt[["refID"]][ind.edit]
  if (out.varFiles) {
    ## Exporting object containing mismatches in the context of modular variants 
    out.file <- 
      file.path(out.dir, 
                paste(out.filename, "_substitutionModPos.rda", sep =""))
    cat("Exporting object containing mismatches in the context of modular ", 
        "variants: \"", out.file, "\".\n", sep = "")
    save(substitution.modPos, file = out.file)
  }
  
  aux.mismatch.modDist <- unlist(mclapply(
    strsplit(unlist(substitution.modPos, use.names = FALSE), split = ":"), 
    function(x) x[1],  mc.cores = num.cores)
  )
  
  mismatch.modDist <- table(
    factor(aux.mismatch.modDist, levels = names(patterns.list))
  )
  
  cat("Plotting distribution of mismatches per position for each modular", 
      "variant. \n") 
  for (i in seq_along(patterns.list)) {
    variant <- names(patterns.list)[i] 
    PlotMismatchesPerVariant(
      substitution.modPos = substitution.modPos, 
      ind = c(aux.mismatch.modDist == variant), mod.len = mod.len, 
      var = variant, tot = aux.modDist[variant], patterns = patterns.list,
      filename = out.filename, out.dir = out.dir, num.cores = num.cores
    )
    
  }
  
  ## Distibution of mismatches per variant normalized by their lengths
  df.mismatch <- data.frame(
    "freq" = mismatch.modDist * 100 / (aux.modDist * mod.len)
  )
  
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, "_mismatchPerVariant_normalizedFreq_graph.pdf",
                    sep = ""))
  cat("Plotting distribution of mismatches per modular variant: \"", out.file,
      "\".\n", sep = "")
  PlotMismatches(data = df.mismatch, x.var = "freq.Var1", y.var = "freq.Freq",
                 out.file = out.file)
  
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, "_mismatchPerVariant_normalizedFreq_table.csv",
                    sep = ""))
  cat("Exporting: \"", out.file, "\" ... \n", sep = "")
  write.csv(as.data.frame(mismatch.modDist), file = out.file)
  
  
  ### 3. Identification of significant variants (Poisson model)
  ### Mismatches per module combination
  mod.dstr <- sort(table(res.sam.filt[["refID"]][ind.edit]), decreasing = TRUE)
  mod.dstr <- mod.dstr[mod.dstr >= min.coverage]
  
  substitution.modComb.poisson <-  mcmapply(
    FUN = PoissonVariantCall, mod.comb = mod.dstr, name = names(mod.dstr), 
    MoreArgs = list(variant = substitution.modPos, error.rate = seq.error), 
    mc.cores = num.cores
  )
  if (out.varFiles) {
    ## Exporting object containing nucleotide variants (Poisson model)
    out.file <- 
      file.path(out.dir, 
                paste(out.filename, "_mismatchPerModComb.rda", sep = ""))
    cat("Exporting object containing nucleotide variants (Poisson model): \"",
        out.file, "\".\n", sep = "")
    save(substitution.modComb.poisson, file = out.file)
  }
  
  aux.var <- unlist(mclapply(
    substitution.modComb.poisson, names, mc.cores = num.cores), 
    use.names = FALSE
  )
  aux.var.len <- unlist(mclapply(
    substitution.modComb.poisson, length, mc.cores = num.cores)
  )
  
  df.substitution.modComb <- data.frame(
    "moduleCombination" = rep(names(aux.var.len), aux.var.len),
    "mismatch"          = aux.var, 
    "count"             = unlist(substitution.modComb.poisson, use.names = FALSE)
  )
  
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, "_mismatchPerModComb_poisson_table.csv", 
                    sep = ""))
  cat("Exporting: \"", out.file, "\" ... \n", sep = "")
  write.csv(df.substitution.modComb, file = out.file)
  
  substitution.modDist.poisson <- factor(unlist(
    mclapply(strsplit(aux.var, split = ":"), function(x) x[1]), use.names = FALSE), 
    levels = names(patterns.list))
  
  substitution.modDist.poisson <- table(substitution.modDist.poisson)
  
  df.mismatch <- data.frame(
    "freq" = as.numeric(substitution.modDist.poisson) * 100 / 
      (aux.modDist * mod.len)
  )
  
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, 
                    "_mismatchPerVariant_normalizedFreq_poisson_graph.pdf",
                    sep = ""))
  cat("Plotting distribution of mismatches per modular variant (Poisson model):",
      " \"", out.file, "\".\n", sep = "")
  PlotMismatches(data = df.mismatch, x.var = "freq.Var1", y.var = "freq.Freq",
                 out.file = out.file)
  
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, 
                    "_mismatchPerVariant_normalizedFreq_poisson_table.csv",
                    sep = ""))
  cat("Exporting: \"", out.file, "\" ... \n", sep = "")
  write.csv(as.data.frame(substitution.modDist.poisson), file = out.file)
  
  
  ### 4. Identification of short indels (pileup approach)
  sourceCpp(file.path(modseq.dir, "src/indelsAlongReferenceSpace.cpp"))
  indel.modPos <- mcmapply(
    FUN = indelsAlongReferenceSpace, cigar_ops = cigar.ops,
    cigar_ops_lengths = cigar.ops.len, mod_len_cumsum = mod.len.cumsum, 
    names_mod_len_cumsum = names.mod.len.cumsum, offset = leftclipped, 
    mc.cores = num.cores
  )
  names(indel.modPos) <- res.sam.filt[["refID"]][ind.edit]
  if (out.varFiles) {
    ## Exporting object containing indels in the context of modular variants
    out.file <- 
      file.path(out.dir, paste(out.filename, "_indelsModPos.rda", sep = ""))
    cat("Exporting object containing indels in the context of modular variants:",
        " \"", out.file, "\".\n", sep = "")
    save(indel.modPos, file = out.file)
  }
  
  aux.indel.modDist <- unlist(mclapply(
    strsplit(unlist(indel.modPos, use.names = FALSE), split = ":"), 
    function(x) x[2], mc.cores = num.cores)
  )
  
  aux.indel.type <- unlist(mclapply(
    strsplit(unlist(indel.modPos, use.names = FALSE), split = ":"),
    function(x) x[1], mc.cores = num.cores)
  )
  
  indel.modDist <- factor(aux.indel.modDist, levels = names(patterns.list))
  
  df.indel <- 
    data.frame("module" = indel.modDist, "type" = aux.indel.type)
  
  df.indel <- as.data.frame(table(df.indel))  # To normalize it
  df.indel[df.indel$Freq == 0, ]$Freq <- NA
  
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, 
                    "_indelPerVariant_normalizedFreq_table.csv", sep = ""))
  cat("Exporting: \"", out.file, "\" ... \n", sep = "")
  write.csv(df.indel, file = out.file)
  
  df.indel$Freq <- 
    as.numeric(df.indel$Freq) * 100 / (rep(aux.modDist, 2) * rep(mod.len, 2))
  
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, 
                    "_IndelPerVariant_normalizedFreq_graph.pdf", sep = ""))
  cat("Plotting distribution of indels per modular variant: \"", out.file,
      "\".\n", sep = "")
  PlotIndels(data = df.indel,  x.var = "module", y.var = "Freq",
             z.var = "type", out.file = out.file)
  
  ### 5. Identification of significant indels as variants (Poisson model)
  ### Indels per module combination
  indel.modComb.poisson <- mcmapply(
    FUN = PoissonVariantCall, mod.comb = mod.dstr, name = names(mod.dstr), 
    MoreArgs = list(variant = indel.modPos, error.rate = seq.error), 
    mc.cores = num.cores
  )
  if (out.varFiles) {
    ## Exporting object containing indels as variants (Poisson model)
    out.file <- 
      file.path(out.dir, 
                paste(out.filename, "_indelsPerModComb.rda", sep = ""))
    cat("Exporting object containing indels as variants (Poisson model): \"",
        out.file, "\".\n", sep = "")
    save(indel.modComb.poisson, file = out.file)
  }
  
  aux.var <- unlist(mclapply(
    indel.modComb.poisson, names, mc.cores = num.cores),
    use.names = FALSE
  )
  aux.var.len <- unlist(mclapply(
    indel.modComb.poisson, length, mc.cores = num.cores)
  )
  
  df.indel.modComb <- data.frame(
    "moduleCombination" = rep(names(aux.var.len), aux.var.len),
    "indel"             = aux.var, 
    "count"             = unlist(indel.modComb.poisson, use.names = FALSE)
  )
  
  out.file <- 
    file.path(out.dir, paste(out.filename, 
                             "_indelPerModComb_poisson_table.csv", sep = ""))
  cat("Exporting: \"", out.file, "\" ... \n", sep = "")
  write.csv(df.indel.modComb, file = out.file)
  
  indel.modDist.poisson <- factor(unlist(
    mclapply(strsplit(aux.var, split = ":"), function(x) x[2])),
    levels = names(patterns.list)
  )
  
  aux.indel.type <- unlist(
    mclapply(strsplit(aux.var, split = ":"), function(x) x[1], 
             mc.cores = num.cores)
  )
  
  df.indel <- 
    data.frame("module" = indel.modDist.poisson, "type" = aux.indel.type)
  
  df.indel <- as.data.frame(table(df.indel))  
  df.indel[df.indel$Freq == 0, ]$Freq <- NA
  
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, 
                    "_indelPerVariant_normalizedFreq_poisson_table.csv",
                    sep = ""))
  cat("Exporting: \"", out.file, "\" ... \n", sep = "")
  write.csv(df.indel, file = out.file)
  
  df.indel$Freq <- 
    as.numeric(df.indel$Freq) * 100 / (rep(aux.modDist, 2) * rep(mod.len, 2))
  
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, 
                    "_indelPerVariant_normalizedFreq_poisson_graph.pdf", sep = ""))
  cat("Plotting distribution of indels per modular variant: \"", out.file,
      "\".\n", sep = "")
  PlotIndels(data = df.indel,  x.var = "module", y.var = "Freq",
             z.var = "type", out.file = out.file)
  
}


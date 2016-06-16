# **TODO**: Exclude duplicates for the distribution (res.sam.real: 1024, 1040)?
AlignmentsFiltering <- function(mod.comb, res.sam.realn, bwa.dupl=TRUE,
                                coverage.left=0, coverage.right=0, mapQ.thold=8, 
                                editDist.thold=8, num.reads=NULL, modseq.dir=NULL, 
                                num.cores=numeric(0)) {
  
  ## Function arguments
  #  mod.comb     Variable of type DNAStringSet containing sequences of the module
  #               combinations. Alternatively path to fasta file.
  
  ## Whenever modseq path or output directory are not specified, assumed to be 
  #  the current working directory
  if (is.null(modseq.dir)) {
    modseq.dir <- getwd()
    warning("Object \'modseq.dir\' not found, set to: \"", modseq.dir, "\".")
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
  
  ## Loading results from read mapping
  if (length(res.sam.realn) == 1) {
    source(file.path(modseq.dir, "R/functions/LoadAlignmentFile.R"))
    res.sam.realn <- 
      LoadAlignmentFile(sam.file = res.sam.realn, bwa.dupl = bwa.dupl)
  }
  
  ## Mapped reads (excluding secondary alignments)
  #  0 and 16      -> mapped read (forward or reverse strand, respectively)
  #  1024 and 1040 -> reads marked as duplicates (forward or reverse strand, 
  #                   respectively)
  ind.sam.realn <- 
    which(res.sam.realn[["flag"]] == 0 | res.sam.realn[["flag"]] == 16 |
            res.sam.realn[["flag"]] == 1024 | res.sam.realn[["flag"]] == 1040)
  
  ### FILTERING 
  # 1. By coverage (first nucleotide consensus 1 and last of consensus 2)
  leftclipped <- as.integer(res.sam.realn[["refPOS"]][ind.sam.realn] - 1) 
  rightclipped <- 
    width(mod.comb[res.sam.realn[["refID"]][ind.sam.realn]]) - 
    cigarWidthAlongReferenceSpace(res.sam.realn[["CIGAR"]][ind.sam.realn]) -
    leftclipped
  ind.coverage <- 
    which(leftclipped <= coverage.left & rightclipped <= coverage.right) 
  
  checkEmpty(ind.coverage)
  printResults("coverage", ind.coverage, num.reads)
  
  # 2. By mapping quality
  ind.mapQ <- 
    which(as.numeric(res.sam.realn[["mapQ"]][ind.sam.realn][ind.coverage]) >=
            mapQ.thold)
  
  checkEmpty(ind.mapQ)
  printResults("mapping quality", ind.mapQ, num.reads)
  
  # 3. By edit distance
  ## Correcting entry containign mismatching positions 
  md.aux <- unlist(mclapply(
    strsplit(res.sam.realn[["mismatchPOS"]][ind.sam.realn][ind.coverage][ind.mapQ],
             split = ":"), function(x) x[[1]], mc.cores = num.cores)
    )
  md.aux2 <- unlist(mclapply(
    strsplit(res.sam.realn[[13]][ind.sam.realn][ind.coverage][ind.mapQ], split = ":"),
    function(x) x[[1]], mc.cores = num.cores))
  
  res.sam.realn[["mismatchPOS"]][ind.sam.realn][ind.coverage][ind.mapQ][
    md.aux == "XA"] <- 
    res.sam.realn[[13]][ind.sam.realn][ind.coverage][ind.mapQ][
      md.aux == "XA"]

  res.sam.realn[["mismatchPOS"]][ind.sam.realn][ind.coverage][ind.mapQ][
    md.aux == "SA" & md.aux2 == "XA"] <- 
    res.sam.realn[[14]][ind.sam.realn][ind.coverage][ind.mapQ][
      md.aux == "SA" & md.aux2 == "XA"]
  
  res.sam.realn[["mismatchPOS"]][ind.sam.realn][ind.coverage][ind.mapQ][
    md.aux == "SA" & md.aux2 == "MD"] <- 
    res.sam.realn[[13]][ind.sam.realn][ind.coverage][ind.mapQ][
      md.aux == "SA" & md.aux2 == "MD"]
  
  ## Correcting edit-distance entry
  if (bwa.dupl) {
    
    res.sam.realn[["editDistance"]][ind.sam.realn][ind.coverage][ind.mapQ][
      md.aux == "XA"] <- 
      res.sam.realn[[16]][ind.sam.realn][ind.coverage][ind.mapQ][
        md.aux == "XA"] 
    
    res.sam.realn[["editDistance"]][ind.sam.realn][ind.coverage][ind.mapQ][
      md.aux == "SA" & md.aux2 == "XA"] <-
      res.sam.realn[[17]][ind.sam.realn][ind.coverage][ind.mapQ][
        md.aux == "SA" & md.aux2 == "XA"]   
    
    res.sam.realn[["editDistance"]][ind.sam.realn][ind.coverage][ind.mapQ][
      md.aux == "SA" & md.aux2 == "MD"] <- 
      res.sam.realn[[16]][ind.sam.realn][ind.coverage][ind.mapQ][
        md.aux == "SA" & md.aux2 == "MD"]

  } else {
    
    res.sam.realn[["editDistance"]][ind.sam.realn][ind.coverage][ind.mapQ][
      md.aux == "XA"] <- 
      res.sam.realn[[15]][ind.sam.realn][ind.coverage][ind.mapQ][
        md.aux == "XA"]   
    
    res.sam.realn[["editDistance"]][ind.sam.realn][ind.coverage][ind.mapQ][
      md.aux == "SA" & md.aux2 == "XA"] <-
      res.sam.realn[[16]][ind.sam.realn][ind.coverage][ind.mapQ][
        md.aux == "SA" & md.aux2 == "XA"]
    
    res.sam.realn[["editDistance"]][ind.sam.realn][ind.coverage][ind.mapQ][
      md.aux == "SA" & md.aux2 == "MD"] <- 
      res.sam.realn[[15]][ind.sam.realn][ind.coverage][ind.mapQ][
        md.aux == "SA" & md.aux2 == "MD"] 
    
  }
  
  editDistance <- 
    as.integer(unlist(mclapply(strsplit(
      res.sam.realn[["editDistance"]][ind.sam.realn][ind.coverage][ind.mapQ],
      split = "NM:i:"), function(x) x[2], mc.cores = num.cores)))
  ind.exactMatch <- which(editDistance == 0)
  ind.edit <- which(editDistance <= editDist.thold)
  
  checkEmpty(ind.edit)
  printResults("mapping quality", ind.edit, num.reads)
  
  res.sam.filt <- 
    mclapply(res.sam.realn, function(x) x[ind.sam.realn][ind.coverage][
      ind.mapQ], mc.cores = num.cores)
  
  return(list(res.sam.filt, editDistance, ind.edit))
}

checkEmpty <- function(ind) {
  
  if (length(ind) == 0) {
    stop("None of the alignments meets the filtering criteria")
  }
}

printResults <- function(filter, ind, num.reads) {
  
  if (!is.null(num.reads)) {
    cat("After filtering alignments by ", filter, ": found ", length(ind), 
        " module combinations in ", num.reads, " reads (", 
        round(length(ind.edit) * 100 / num.reads, digits = 2), "%).\n", 
        sep = "")
  } else {
    cat("After filtering alignments by ", filter, ": found ", length(ind.edit), 
        " module combinations.\n", sep = "")
  }
  
}


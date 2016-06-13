# **TODO**: Exclude duplicates for the distribution (res.sam.real: 1024, 1040)?
AlignmentsFiltering <- function(mod.comb, res.sam.realn, bwa.dupl=TRUE,
                                coverage.left=0, coverage.right=0, mapQ.thold=8, 
                                editDist.thold=8, num.reads=NULL, modseq.dir=NULL, 
                                num.cores=numeric(0)) {
  
  ## Whenever modseq path or output directory are not specified, assumed to be 
  #  the current working directory
  if (is.null(modseq.dir)) {
    modseq.dir <- getwd()
    warning("Object \'modseq.dir\' not found, set to: \"", modseq.dir, "\".")
  } 
  
  if (length(num.cores) == 0) {
    num.cores <- detectCores()
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
  
  # 2. By mapping quality
  ind.mapQ <- 
    which(as.numeric(res.sam.realn[["mapQ"]][ind.sam.realn][ind.coverage]) >=
            mapQ.thold)
  
  # 3. By edit distance
  ## Correcting entry containign mismatching positions 
  md.aux <- 
    unlist(mclapply(strsplit(res.sam.realn[["mismatchPOS"]][ind.sam.realn][
      ind.coverage][ind.mapQ], split = ":"), function(x) x[[1]], 
      mc.cores = num.cores))
  md.aux2 <- 
    unlist(mclapply(strsplit(res.sam.realn[[13]][ind.sam.realn][
      ind.coverage][ind.mapQ], split = ":"), function(x) x[[1]], 
      mc.cores = num.cores))
  
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
  
  if (is.null(num.reads)) {
    cat("After filtering: found ", length(ind.edit), " module combinations in ",
        num.reads, " reads (", 
        round(length(ind.edit) * 100 / num.reads, digits = 2), "%).\n", sep = "")
  } else {
    cat("After filtering: found ", length(ind.edit), " module combinations.\n", 
        sep = "")
  }
  
  
  res.sam.filt <- 
    mclapply(res.sam.realn, function(x) x[ind.sam.realn][ind.coverage][
      ind.mapQ], mc.cores = num.cores)
  
  return(list(res.sam.filt, editDistance, ind.edit))
}
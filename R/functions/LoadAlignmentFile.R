LoadAlignmentFile <- function(sam.file, bwa.dupl=TRUE) {
  
  paste(out.dir, out.filename.run3, 
        "_realigned.sam", sep = "")
  if (bwa.dupl) {
    
    res.sam.realn <- 
      scan(file = sam.file, 
           what = list(character(), integer(), character(), integer(), 
                       integer(), character(), character(), integer(), 
                       integer(), character(), character(), character(),
                       character(), character(), character(), character(),
                       character()), 
           skip = 0, flush = TRUE, fill = TRUE, quote = "")
    names(res.sam.realn) <- c("readID", "flag", "refID", "refPOS", "mapQ", 
                              "CIGAR", "", "", "", "readSeq", "readQuality",
                              "mismatchPOS", "", "", "editDistance", "", "")
    
  } else {
    
    res.sam.realn <- 
      scan(file = sam.file, 
           what = list(character(), integer(), character(), integer(),
                       integer(), character(), character(), integer(),
                       integer(), character(), character(), character(),
                       character(), character(), character(), character()),
           skip = 0, flush = TRUE, fill = TRUE, quote = "") 
    names(res.sam.realn) <- c("readID", "flag", "refID", "refPOS", "mapQ", 
                              "CIGAR", "", "", "", "readSeq", "readQuality",
                              "mismatchPOS", "", "editDistance", "", "")

  }  
  
  return(res.sam.realn)
  
}
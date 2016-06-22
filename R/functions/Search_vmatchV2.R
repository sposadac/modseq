SearchPatListID <- function(hitsList.names, hitsList, pat, input.data, mm=0, 
                            mode="f", ambiguous.match=TRUE,
                            num.cores=numeric(0)) {
  
  if (length(num.cores) == 0) {
    num.cores <- detectCores()
  }
  
  patList <- as.matrix(na.omit(pat))
  if (mode == "f") {
    patList <- expand.grid(hitsList.names, patList)
  } else if (mode == "r" ) {
    patList <- expand.grid(patList, hitsList.names)
  } else {
    stop("Invalid entry for \'mode\', \"", mode, "\".\n Object \'mode\' should",
         "be either \"f\" or \"r\".")
  }
  
  patList <- as.matrix(by(data = patList, INDICES = c(1:nrow(patList)), 
                          FUN = function(x) paste(unlist(x), collapse = "")))  
  
  if (hitsList.names == "") {
    hitsList <- unlist(hitsList, use.names = FALSE)
  } else {
    hitsList <- unlist(hitsList[hitsList.names], use.names = FALSE)
  }
  
  # Subject has to be set as an auxiliary varaible to keep the previous list
  # intact
  subject <- input.data[hitsList]  
  retList <- mclapply(
    patList, function(pat, subj, mm, ambiguous.match) 
      SearchPatID(pat, subj, mm, ambiguous.match), subj = subject, mm = mm, 
    ambiguous.match = ambiguous.match, mc.cores= num.cores)
  retList <- mclapply(
    retList, function(x, hits) hits[which(x > 0)], hits = hitsList, 
    mc.cores= num.cores)
  names(retList) <- patList
  
  return(retList)
}

SearchPatID <- function(pat, subj, mm, ambiguous.match){
  
  ## Function arguments
  # pat               Character.
  # subj              Object of class BioString.
  # mm                Integer -- maximum mismatches allowed.
  # ambiguous.match   Bool
  
  # Evaluate empty option
  if (pat != "") {  
    ## **TODO**: make sure it is returning integers
    temp <- 
      vcountPattern(pat, subj, max.mismatch = mm, fixed = !ambiguous.match)
    patRC <- toupper(SF_patternrc(pat))
    temp <- temp + vcountPattern(patRC, subj, max.mismatch = mm, 
                                 fixed = !ambiguous.match)
    return(temp)
    
  } else {
    # If pattern is empty return whole subject as hit -> next search will then 
    # be on whole list again
    return(rep(1, length(sub))) 
  } 
}

PatCounts <- function(list, filename, mod, mod.tot, out.dir=NULL, 
                      num.cores = numeric(0)) {
  
  ## Whenever output directory is not specified, assumed to be the current
  #  working directory
  if (is.null(out.dir)) {
    out.dir <- getwd()
  }
  
  if (length(num.cores) == 0) {
    num.cores <- detectCores()
  }
  
  res.aux <- mclapply(list, length, mc.cores = num.cores)
  if (length(list) != length(unlist(res.aux))) {
    stop("Invalid entry found in \'res.list\'. \n")
  }
  tot <- sum(unlist(res.aux))
  res <- data.frame("seq"    = c(names(res.aux), "total"), 
                    "counts" = c(unlist(res.aux, use.names = FALSE), tot))
  
  out.file <- 
    file.path(out.dir, 
              paste(filename, "_round_", mod, "_of_", mod.tot, ".csv", sep = ""))
  write.csv(res, out.file, row.names = FALSE)
  
  return(tot)
}

#SF_functions provided by S.Schmitt
SF_patternrc <- function(pattern) {
  
  # Make pattern to RC, keep as vector of character
  patternrc <- rev(comp(s2c(pattern))) 
  
  # Replace NA in case there was a undefinded character (e.g. an N) with an N
  patternrc[is.na(patternrc)] <- "N" 
  # Make String again
  patternrc <- c2s(patternrc) 
  
  return(patternrc)
}
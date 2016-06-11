 Run3_gls <- function(patterns=NULL, reads, num.reads=NULL, in.modDir, 
                      mod.filename, res.listName, res.counts.filename, 
                      gls.ambiguity=TRUE, gls.direction="f", gls.mma=0, 
                      mem.trace=FALSE, memTrace=NULL, run.info="ModSeq", 
                      modseq.dir=NULL, out.dir=NULL, num.cores=numeric(0)) {
  
  # **TODO**: verbose level (display info per search round?)
  # **TODO**: map.mode pairwise alignment
  # **TODO**: mismatch allowance per module
  
  ## Function arguments
  # patterns    data frame containing the sequences of the various variants per
  #             module. Columns correspond to different modules and rows 
  #             correspond to variants identifiers. 
  # num.cores   (optional) number of cores available for performing parallel 
  #             tasks.
  
  ## Whenever modseq path or output directory are not specified, assumed to be 
  #  the current working directory
  if (is.null(modseq.dir)) {
    modseq.dir <- getwd()
    warning("Object \'modseq.dir\' not found, set to: \"", modseq.dir, "\".")
  } 
  if (is.null(out.dir)) {
    out.dir <- getwd()
  }
   
  source(file.path(modseq.dir, "R/functions/Search_vmatchV2.R"))
  source(file.path(modseq.dir, "R/functions/PlotModuleCounts.R"))
  source(file.path(modseq.dir, "R/functions/VariantCounts.R"))
  
  if (length(num.cores) == 0) {
    num.cores <- detectCores()
  }
  
  if (!is.data.frame(patterns)) {
    source(file.path(modseq.dir, "R/functions/LoadModuleTable.R"))
    patterns <- 
      LoadModuleTable(in.modulesDir = in.modDir, modules.filename = mod.filename) 
  }
  mod.tot <- ncol(patterns)
  
  if (is.null(num.reads)) {
    num.reads <- length(reads)
  }
  
  ## Pattern search
  ti.search <- Sys.time()
  res.list <- list()
  res.list[[1]] <-  names(reads)
  names(res.list)[1] <- "" 
  res.counts <- numeric(mod.tot + 1)
  res.counts[1] <- num.reads
  if (mem.trace) {
    if (!exists("MemTrace")) {
      source(file.path(modseq.dir, "R/functions/MemTrace.R"))
    }
    memTrace <- c(memTrace, MemTrace())
    offsetMem <- length(memTrace)
    memList <- format(object.size(res.list), units = "MB")
  }
  cat("GLS - Starting pattern search ...\n")
  
  if (gls.direction == "f") {
    
    cat("Search mode: Forward - gls. \n")
    for (i in seq_len(mod.tot)) {
      
      res.list <- mclapply(
        names(res.list), FUN = SearchPatListID, hitsList = res.list, 
        pat = patterns[i], input.data = reads, mm = gls.mma,  
        mode = "f", ambiguous.match = gls.ambiguity, num.cores = num.cores, 
        mc.cores = num.cores)
      res.list <- unlist(res.list, recursive = FALSE)
      res.counts[i+1] <- PatCounts(res.list, res.counts.filename, i, mod.tot,
                                   num.cores) 
      cat("Mod. ", i, "/", mod.tot, "\t", res.counts[i+1], "\n", sep = "")
      
      if (mem.trace) {
        memList <- c(memList,format(object.size(res.list), units = "MB"))
        memTrace <- c(memTrace, MemTrace())
      }
      
    }
    
  } else if (gls.direction == "r") {
    
    cat("Search mode: Reverse - gls. \n")
    for (i in mod.tot:1) {
      
      res.list <- mclapply(
        names(res.list), FUN = SearchPatListID, hitsList = res.list, 
        pat = patterns[i], input.data = reads, mm = gls.mma, mode = "r", 
        ambiguous.match = gls.ambiguity, num.cores = num.cores, 
        mc.cores = num.cores)
      res.list <- unlist(res.list, recursive = FALSE)
      res.counts[i+1] <- PatCounts(res.list, res.counts.filename, i, 
                                   mod.tot, num.cores)
      cat("Mod. ", i, "/", mod.tot, "\t", res.counts[i+1], "\n", sep = "")
      
      if (mem.trace) {
        memList <- c(memList,format(object.size(res.list), units = "MB"))
        memTrace <- c(memTrace, MemTrace())
      } 
      
    }
  }
  
  cat("Runtime:", round(as.numeric(
    difftime(Sys.time(), ti.search, units = "mins")), digits = 4), "min.\n")
  cat("Pattern search completed!\n")
  
  ## Exporting list containing results of the read mapping as R object 
  out.file <- 
    file.path(out.dir, paste("resList_", res.listName, ".rda", sep = ""))
  cat("Exporting list containing results of the read mapping: \"", out.file,
      "\".\n", sep = "")
  ## **TESTING**
  #save(res.list, file = out.file)
  
  if (mem.trace) {
    memTrace <- c(memTrace, MemTrace())
  }
  
  cat("Saving and plotting results of the read mapping...\n")
  ## Reformatting counts as data.frame
  res.counts <- data.frame(
    "round" = 
      c("Input-data", paste("Mod. ", seq_len(mod.tot), "/", mod.tot, sep ="")),
    "counts" = res.counts)
  
  if (gls.direction == "r") { 
    
    res.counts$round <- 
      factor(res.counts$round, 
             levels = c(levels(res.counts$round)[1], 
                        rev(levels(res.counts$round)[2:(mod.tot+1)])))
    counts <- c(res.counts$counts[1], rev(tail(res.counts$counts, n = mod.tot)))
    cat("Found ", res.counts[2, 2], " module combinations in ", num.reads,
        " reads (", round(res.counts[2, 2] * 100 / num.reads, digits = 2), 
        "%).\n", sep = "")
    
  } else {
    
    counts <- res.counts$counts
    mod.comb.tot <- tail(res.counts, n = 1)[[2]] #res.counts[nrow(res.counts), 2]
    cat("Found ", mod.comb.tot, " module combinations in ", num.reads, " reads (", 
        round(mod.comb.tot * 100 / num.reads, digits = 2), "%).\n", sep = "")
    
  }
  
  ### Module-counts plot
  out.file <- 
    file.path(out.dir, paste(res.counts.filename, "_graph.pdf", sep = ""))
  cat("Printing module-counts plot as: \"", out.file, "\".\n", sep = "")
  PlotModuleCounts(res.counts, x.var = "round", y.var = "counts", 
                   ymax = max(res.counts$counts), out.file, 
                   plot.label = run.info, gls.ambiguity = gls.ambiguity)
  
  ### Distribution of modular variants per search round
  VariantCounts(patterns = patterns, num.reads = num.reads,
                counts = counts, labels = levels(res.counts$round), 
                file.prefix = res.counts.filename, out.dir = out.dir,
                gls.direction = gls.direction, modseq.dir = modseq.dir, 
                num.cores = num.cores)

#   if (sum(gls.mma > 0) > 0) {
#     cat("Maximum mismatch allowance for module", which(gls.mma > 0), "was", 
#         gls.mma[which(gls.mma > 0)], ".\n")
#   }
  
  ### Plot of memory consumption
  if (mem.trace) {
    source(file.path(modseq.dir, "R/functions/PlotMemProfile.R"))
    PlotMemProfile(memList = memList, memTrace = memTrace, offsetMem = offsetMem,
                   num.mod = mod.tot, out.file = res.listName, out.dir = out.dir)
  }
  
  ### Plot of repeated hits (reads that were mapped to several 
  #   module-combination)
  aux <- unlist(res.list, use.names = FALSE)
  
  if (sum(duplicated(aux)) > 0) {
    
    source(file.path(modseq.dir, "R/functions/RepeatedHits.R"))
    retList <- RepeatedHits(res.list, data = aux, patterns = patterns, 
                            num.reads = num.reads, in.modDir, mod.filename,
                            res.listName, out.dir, modseq.dir = modseq.dir,
                            num.cores = num.cores)
    mod.comb <- retList[[1]]
    res.list.lengths <- retList[[2]]
    res.list.name <- retList[[3]]
    
    return(list(res.list, mod.comb, res.list.lengths, res.list.name))
    
  } else {
    
    cat("All mapped reads are mapped exactly 1 time (i.e., mapped to 1 ",
        "module-combination). \n")
    return(list(res.list))
    
  }

}
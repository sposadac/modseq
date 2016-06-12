### MAIN SCRIPT ###
# Authors: 
# Susana Posada Cespedes <susana.posada@bsse.ethz.ch>
# Steven Schmitt
####

### DELETE WORKSPACE 
rm(list=ls())

## Get current work directory
wdir <- getwd()

## LOAD USER-DEFINED OPTIONS AND PACKAGES
# **TODO** add R and R/functions directory to search path
# Path to modseq head directory. At the moment, set to current directory
modseq.dir <- wdir
source(file.path(modseq.dir, "R/modseq_input.R"))
source(file.path(modseq.dir, "R/functions/InputCheck.R")) 
source(file.path(modseq.dir, "R/functions/NameGen.R"))
names(run) <- c(1, 2, 3, 4, 5)

### TRACKING TOTAL MEMORY USAGE
# MemTrace() Memory in Megabytes
memTrace <- NULL
if (mem.trace) {
  source(file.path(modseq.dir, "R/functions/MemTrace.R"))
  memTrace <- MemTrace()
}

### START RUN
Ti <- Sys.time()

# Log-file
log.file <- file.path(out.dir, 
                      paste("out_", out.filename.run2, "_run",
                            paste(names(run)[which(run == 1)], collapse = "_"),
                            ".txt", sep = ""))
if (run[3] == 1 | run[4] == 1) {
  log.file <- file.path(out.dir, 
                        paste("out_", out.filename.run3, "_run",
                              paste(names(run)[which(run == 1)], collapse = "_"),
                              ".txt", sep = ""))
}
if (file.exists(log.file)) {
  output <- file(log.file, open = "at")
} else {
  output <- file(log.file, open = "wt")
}
sink(output, append = TRUE)

cat("Run", paste(names(run)[which(run == 1)], collapse = ","), "started...\n")
keep <- ls()

# Read fastq-files
if (run[1] == 1 || (run[2] == 1 && qtrim.flag == 0)) { 
  
  source(file.path(modseq.dir, "R/functions/PlotQualityDistribution.R"))
  source(file.path(modseq.dir, "R/functions/PlotReadLengthDistribution.R"))
  
  cat("Sequencing mode: \"", seq.mode, "\".\n", sep = "")
  if (seq.mode == "PE") {
    
    readF1 <- readFastq(dirPath = in.seqDir, pattern = forward.filename)
    readF2 <- readFastq(dirPath = in.seqDir, pattern = reverse.filename)
    num.reads <- length(readF1)
    # checks
    if (length(readF2) != num.reads) {
      stop("Number of reads in input files does not match: \n Forward set: ",
           num.reads, " reads, \n Reverse set: ", length(readF2), " reads.\n")
    }
    cat("Number of reads: ", num.reads, " x 2 (mode: PE).\n", sep = "")
    PlotQualityDistribution(readF1, forward.filename, out.dir)
    PlotQualityDistribution(readF2, reverse.filename, out.dir)
    PlotReadLengthDistribution(readF1, forward.filename, out.dir)
    PlotReadLengthDistribution(readF2, reverse.filename, out.dir)
    
  } else if (seq.mode == "SE") {
    
    readF1 <- readFastq(dirPath = in.seqDir, pattern = in.filename)
    num.reads <- length(readF1)
    cat("Number of reads: ", num.reads, ".\n", sep = "")
    PlotQualityDistribution(readF1, in.filename, out.dir)
    PlotReadLengthDistribution(readF1, in.filename, out.dir)
    
  }
}

###############################################################################
## 1: QUALITY TRIMMING
###############################################################################
if (run[1] == 1) {  
  cat("====================================================================\n")  
  cat("Quality trimming. \n")
  cat("====================================================================\n")  
  
  
  ### QUALITY TRIMMING: 
  # Uncalled bases (i.e. N) are removed from either both ends (qtrim.3end = 0)
  # or only the right end (qtrim.3end = 1).
  # Bases are removed if quality score is less than or equal to qtrim.thold.
  source(file.path(modseq.dir, "R/functions/Run1_qualTrimming.R"))
  
  if (seq.mode == "SE") {
    
    readF1.nQtrim <- Run1_qualTrimming(
      readF1 = readF1, in.filename = in.filename, seq.mode = seq.mode, 
      qtrim.3end = qtrim.3end, qtrim.thold = qtrim.thold, 
      out.filename = out.filename.run1, out.ssplot = out.ssplot, 
      wdir = wdir, modseq.dir = modseq.dir, out.dir = out.dir)

  } else if (seq.mode == "PE") {
    
    reads.nQtrim <- Run1_qualTrimming(
      readF1 = readF1, in.filename = forward.filename, seq.mode = seq.mode,
      qtrim.3end = qtrim.3end, qtrim.thold = qtrim.thold, 
      out.filename = out.filename.run1, out.ssplot = out.ssplot,
      readF2 = readF2, reverse.filename = reverse.filename, wdir = wdir, 
      modseq.dir = modseq.dir, out.dir = out.dir)
    readF1.nQtrim <- reads.nQtrim[[1]]
    readF2.nQtrim <- reads.nQtrim[[2]]
    
  }

  qtrim.flag <- as.integer(1) 
  cat("Object \'qtrim.flag\' is set to ", qtrim.flag,".\n", sep = "")
  
  if (run[2] == 1) {
    keep <- append(keep, c("num.reads", "readF1.nQtrim", "IlluminaStat"))
    if (seq.mode == "PE") {
      keep <- append(keep, c("readF2.nQtrim"))
    }
  }
  if (mem.trace) {
    memTrace <- c(memTrace, MemTrace())
  }
  rm(list=ls()[!(ls() %in% keep)])
  
}

###############################################################################
## 2: PAIRED-END READ ASSEMBLY
###############################################################################
if (run[2] == 1) {
  cat("========================================================================\n")
  cat("Paired-end read assembly.\n")
  cat("========================================================================\n")
  
  if (!exists("keep")) {
    keep <- ls()
  }

  if (mem.trace) {
    memTrace <- c(memTrace, MemTrace())
  }
  
  if (!exists("num.reads")) {
    num.reads <- NULL
  }
  
  ## PE read assembly
  source(file.path(modseq.dir, "R/functions/Run2_peAssembly.R"))
  if (qtrim.flag == 0) {
    
    num.reads <- Run2_peAssembly(
      readF1 = readF1, readF2 = readF2, qtrim.flag = qtrim.flag, 
      forward.file = forward.filename, reverse.file = reverse.filename, 
      out.ssplot = out.ssplot, out.filename.run2 = out.filename.run2,
      len.rawData = num.reads, in.seqDir = in.seqDir,  wdir = wdir,
      modseq.dir = modseq.dir, pandaseq.path = pandaseq.path, out.dir = out.dir
      )
    
    rm(list = c("readF1", "readF2"))
    
  } else if (qtrim.flag == 1) {

    in.tr1 <- file.path(wdir, "data/processed", 
                        paste(out.filename.run1, "_1_nQtrim.fastq", sep = ""))
    in.tr2 <- file.path(wdir, "data/processed",
                        paste(out.filename.run1, "_2_nQtrim.fastq", sep = ""))
    
    ## Check if files are located in the defatult directory, 
    #  i.e. 'data/processed/'
    if (!file.exists(in.tr1) || !file.exists(in.tr2)) {
      warning("Files \"", in.tr1, "\" and \"", in.tr2, "\" where not found.\n", 
              "Instead, loading trimmed-reads using user-defined options: \'",
              file.path(in.seqDir, paste(forward.filename, ".fastq", sep = "")), 
              "\' and \'", 
              file.path(in.seqDir, paste(reverse.filename, ".fastq", sep = "")),
              "\'.\n")
      in.tr1 <- file.path(in.seqDir, paste(forward.filename, ".fastq", sep = ""))
      in.tr2 <- file.path(in.seqDir, paste(reverse.filename, ".fastq", sep = ""))
    }
    if (!exists("readF1.nQtrim")) {
      readF1.nQtrim <- readFastq(in.tr1) 
    } 
    if (!exists("readF2.nQtrim")) {
      readF2.nQtrim <- readFastq(in.tr2) 
    }
    
    num.reads <- Run2_peAssembly(
      readF1 = readF1.nQtrim, readF2 = readF2.nQtrim, qtrim.flag = qtrim.flag,
      forward.file = in.tr1, reverse.file = in.tr2, out.ssplot = out.ssplot, 
      out.filename.run2 = out.filename.run2, in.filename.run1 = out.filename.run1,
      len.rawData = num.reads, in.seqDir = in.seqDir, wdir = wdir, 
      modseq.dir = modseq.dir, pandaseq.path = pandaseq.path, out.dir = out.dir
      )
    
    rm(list = c("readF1.nQtrim", "readF2.nQtrim"))
    
  }
  paired.flag <- as.integer(1)
  cat("Object \'paired.flag\' is set to ", paired.flag, ". \n", sep = "")
  
  if (mem.trace) {
    memTrace <- c(memTrace, MemTrace())
  }
  
  if (run[3] == 1) {
    keep <- append(keep, "num.reads")
  }
  
  rm(list=ls()[!(ls() %in% keep)])
  if (exists("IlluminaStat")) {
    rm(IlluminaStat)
  }
  
}

###############################################################################
## 3: PATTERN SEARCH: Linear search or read mapping (bwa mem)
###############################################################################
if (run[3] == 1) {
  cat("========================================================================\n")
  cat("Read mapping. \n")
  cat("========================================================================\n")
  
  if (!exists("keep")) {
    keep <- ls()
  } 
  
  num.cores <- detectCores()
  ## Loading module table
  source(file.path(modseq.dir, "R/functions/LoadModuleTable.R"))
  patterns <- 
    LoadModuleTable(in.modulesDir = in.modDir, modules.filename = mod.filename,
                    list = FALSE)
  
  source((file.path(modseq.dir, "R/functions/GetReadsFile.R")))
  reads.file <- GetReadsFile(wdir)
  
  if (map.mode == "gls") {
    
    if (file.exists(reads.file[[1]])) {
      reads.subj <- readDNAStringSet(
        file = reads.file[[1]], format = reads.file[[2]], use.names = TRUE) 
      cat("Input reads file: \"", reads.file[[1]], "\".\n", sep = "")
      
    } else {
      stop("File \"", reads.file, "\" not found.\n")
      
    }
    num.reads <- length(reads.subj)
    
    source(file.path(modseq.dir, "R/functions/Run3_gls.R"))
    retList <- Run3_gls(
      patterns=patterns, reads=reads.subj, num.reads=num.reads, in.modDir, 
      mod.filename, res.listName, res.counts.filename,  
      gls.ambiguity = gls.ambiguity, gls.direction = gls.direction, 
      gls.mma = gls.mma, mem.trace = mem.trace, memTrace = memTrace, 
      run.info = run.info, modseq.dir = modseq.dir, out.dir = out.dir, num.cores = num.cores)
    
    res.list <- retList[[1]]
    if (length(retList) > 1) {
      mod.comb <- retList[[2]]
      res.list.lengths <- retList[[3]]
      res.list.ids <- retList[[4]]
      if (run[4] == 1) {
        keep <- append(keep, c("mod.comb", "res.list.lengths", "res.list.ids", 
                               "GetVariantsNames")) 
      }
    } 
  
    if (run[4] == 1) {
      keep <- append(keep, "res.list") 
    }
    
  } else if (map.mode == "bwa") {
    
    reads.file <- reads.file[[1]]
    if (!file.exists(reads.file)) {
      stop ("File \"", reads.file, "\" not found.\n")
    }
    
    source(file.path(modseq.dir, "R/functions/Run3_bwa.R"))
    sam.file <- 
      Run3_bwa(reads.file = reads.file, in.modDir = in.modDir, 
               mod.filename = mod.filename, bwa.path = bwa.path, 
               bwa.cVal = bwa.cVal, bwa.dupl = bwa.dupl, gatk.path = gatk.path, 
               picard.path = picard.path, samtools.path = samtools.path, 
               out.filename.run2 = out.filename.run2,  
               out.filename.run3 = out.filename.run3, 
               patterns = patterns, modseq.dir = modseq.dir, out.dir = out.dir,
               num.cores = num.cores)
    
    if (run[4] == 1) {
      
      source(file.path(modseq.dir, "R/functions/LoadAlignmentFile.R"))
      res.sam.realn <- LoadAlignmentFile(sam.file = sam.file, bwa.dupl = bwa.dupl)
      keep <- append(keep, "res.sam.realn")
      
    } else {
      
      res.sam.realn <- 
        scan(file = sam.file, what = list(character(), integer()), skip = 0,
             flush = TRUE, fill = TRUE, quote = "") 
      names(res.sam.realn) <- c("readID", "flag")
      
    }
    
    if (!exists("num.reads")) {
      reads.subj <- readFastq(reads.file)
      num.reads <- length(reads.subj)
    }
    
    ## Alignment stats.
    cat("Number of unmapped reads: ", sum(res.sam.realn[["flag"]] == 4), ".\n",
        sep ="")
    cat("Number of secondary alignments: ", 
        sum(res.sam.realn[["flag"]] == 256 | res.sam.realn[["flag"]] == 272), 
        ".\n", sep ="")
    
    ## Indentifying PCR or optical duplicates
    #  (1024, 1040 - forward strand and reverse strand, respec.)
    reads.dupl <- 
      sum(res.sam.realn[["flag"]] == 1024 | res.sam.realn[["flag"]] == 1040)
    ## Mapped reads (ignoring secondary alignments, flags 256 and 272)
    #  (0 - mapped in forward strand, 16 - mapped in reverse strand)
    reads.map <- 
      sum(res.sam.realn[["flag"]] == 0 | res.sam.realn[["flag"]] == 16) + 
      reads.dupl
    cat("Found ", reads.map, " module combinations in ", num.reads, 
        " reads (", round(reads.map * 100 / num.reads, digits = 2), 
        "%). \n", sep = "")
    cat("Found ", reads.dupl, " marked reads as duplicates (", 
        round(reads.dupl * 100 / num.reads, digits = 2), "%). \n", 
        sep = "") 
  }
  
  if (run[4] == 1 || run[5] == 1) {
    keep <- append(keep, c("num.reads", "LoadModuleTable"))
    
    if (existsFunction("ModuleCombinationsGen")){
      keep <- append(keep, "ModuleCombinationsGen")
    }
    
  }
  
  if (mem.trace) {
    memTrace <- c(memTrace, MemTrace())
  }
  rm(list=ls()[!(ls() %in% keep)]) 
}

###############################################################################
## 4: LIBRARY COMPOSITION
###############################################################################
if (run[4] == 1) {
  cat("========================================================================\n")
  cat("Library composition. \n")
  cat("========================================================================\n")
  
  if (!exists("keep")) {
    keep <- ls()
  }
  
  num.cores <- detectCores()
  
  patterns <- file.path(in.modDir, paste(mod.filename, ".csv", sep = ""))
  
  ## Generating module combinations and loading fasta file 
  if (!exists("mod.comb")) {    
    mod.file <- 
      file.path(in.modDir, paste(mod.filename, "_modComb.fasta", sep = ""))
    
    if (!file.exists(mod.file)) {
      
      if (!existsFunction("ModuleCombinationsGen")) {
        source(file.path(modseq.dir, "R/functions/ModuleCombinationsGen.R"))
      }
      
      ModuleCombinationsGen(modules.filename = mod.filename, pattern = patterns, 
                            in.modDir = in.modDir, modseq.dir = modseq.dir, 
                            num.cores = num.cores)
    }
    
    cat("Reading reference sequences ... \n")
    mod.comb <- readDNAStringSet(mod.file)
    
  }
  
  if (map.mode == "gls") {
    
    if (!exists("res.list")) {
      res.list <- NULL
    }   
    if (!exists("res.list.lengths")) {
      res.list.lengths <- NULL
    }
    if (!exists("res.list.ids")) {
      res.list.ids <- NULL
    }
    
    source(file.path(modseq.dir, "R/functions/Run4_gls.R"))
    Run4_gls(
      patterns = patterns, mod.comb = mod.comb, res.listName = res.listName, 
      res.list = res.list, res.list.ids = res.list.ids, 
      res.list.lengths = res.list.lengths, out.filename = out.filename.run3,
      modseq.dir = modseq.dir, out.dir = out.dir, num.cores = num.cores
      )

  } else if (map.mode == "bwa") { 
    
    
    if (!exists("num.reads")) {
      num.reads <- NULL
    }
    if (!exists("res.sam.realn")) {
      res.sam.realn <- 
        file.path(out.dir, paste(out.filename.run3, "_realigned.sam", sep = ""))
    }
    
    source(file.path(modseq.dir, "R/functions/Run4_bwa.R"))
    retList <- Run4_bwa(
      patterns = patterns, mod.comb = mod.comb, res.sam.realn = res.sam.realn,
      bwa.dupl = bwa.dupl, coverage.left = coverage.left, 
      coverage.right = coverage.right, mapQ.thold = mapQ.thold, 
      editDist.thold = editDist.thold, num.reads = num.reads,
      out.filename = out.filename.run3,  modseq.dir = modseq.dir,
      out.dir = out.dir, num.cores = num.cores
      )
    
    res.sam.filt <- retList[[1]]
    editDistance <- retList[[2]] 
    aux.modDist <- retList[[3]]
    
    if (run[5] == 1) {
      keep <- append(keep, c("mod.comb", "res.sam.filt", "editDistance", 
                             "aux.modDist"))
    }
  }
  
  if (mem.trace) {
    memTrace <- c(memTrace, MemTrace())
  }
  
  rm(list=ls()[!(ls() %in% keep)]) 
  
  if (exists("res.list")) {
    rm(res.list)
  }
}

###############################################################################
## 5: VARIANT CALLING
###############################################################################
if (run[5] == 1) {
  cat("========================================================================\n")
  cat("Detection of nucleotide variants and short indels. \n")
  cat("========================================================================\n")
  
  if (!exists("keep")) {
    keep <- ls()
  }
  
  num.cores <- detectCores()
  if (!exists("num.reads")) {
    num.reads <- NULL
  }
  
  if (map.mode == "gls") {
    ## **TODO**: variant calling (mismatches and indels), when gls.mma > 0
    
  } else if (map.mode == "bwa") {
    ## **TODO**: verbose level
    
    ## Loading module table
    if (!existsFunction("LoadModuleTable")) {
      source(file.path(modseq.dir, "R/functions/LoadModuleTable.R"))
    }
    
    retList <- LoadModuleTable(in.modulesDir = in.modDir,
                               modules.filename = mod.filename, list = TRUE) 
    patterns <- retList[[1]]
    patterns.list <- retList[[2]]
    
    ## Generating module combinations and loading fasta file 
    if (!exists("mod.comb")) {    
      mod.file <- 
        file.path(in.modDir, paste(mod.filename, "_modComb.fasta", sep = ""))
      
      if (!file.exists(mod.file)) {
        
        if (!existsFunction("ModuleCombinationsGen")) {
          source(file.path(modseq.dir, "R/functions/ModuleCombinationsGen.R"))
        }
        
        ModuleCombinationsGen(modules.filename = mod.filename, pattern = patterns, 
                              in.modDir = in.modDir, modseq.dir = modseq.dir, 
                              num.cores = num.cores)
      }
      
      cat("Reading reference sequences ... \n")
      mod.comb <- readDNAStringSet(mod.file)
      
    }
    mod.comb.len <- length(mod.comb)
    
    if (!exists("res.sam.filt")) {
      
      if (!existsFunction("AlignmentsFiltering")) {
        source(file.path(modseq.dir, "R/functions/AlignmentsFiltering.R"))
      }
      
      if (!exists("res.sam.realn")) {
        res.sam.realn <- 
          file.path(out.dir, paste(out.filename.run3, "_realigned.sam", sep = ""))
      }
      
      retList <- AlignmentsFiltering(
        mod.comb = mod.comb, res.sam.realn = res.sam.realn, bwa.dupl = bwa.dupl,
        coverage.left = coverage.left, coverage.right = coverage.right, 
        mapQ.thold = mapQ.thold, editDist.thold = editDist.thold, 
        num.reads =  num.reads, modseq.dir = modseq.dir, num.cores = num.cores
      )
      res.sam.filt <- retList[[1]]
      editDistance <- retList[[2]]
      ind.edit <- retList[[3]]
      
      if (!exists("aux.modDist")) {
        aux.modDist <- table(
          factor(unlist(strsplit(res.sam.filt[["refID"]][ind.edit], split = ":")), 
                 levels = names(patterns.list))
        )
      }
      
    }
    
    mod.len <- nchar(patterns.list)
    
    ### 1. "Soft Clipping" reads and quality strings 
    cigar.ops <- 
      explodeCigarOps(res.sam.filt[["CIGAR"]][
        editDistance <= editDist.thold & editDistance > 0]) 
    cigar.ops.len <- 
      explodeCigarOpLengths(res.sam.filt[["CIGAR"]][
        editDistance <= editDist.thold & editDistance > 0])
    
    source(file.path(modseq.dir, "R/functions/SeqTrimming.R"))
    reads.trimmed <- 
      DNAStringSet(mcmapply(
        SeqTrimming, cigarOps = cigar.ops, cigarOpLengths = cigar.ops.len, 
        seq = res.sam.filt[["readSeq"]][
          editDistance <= editDist.thold & editDistance > 0], 
        mc.cores = num.cores))
  
    ## Exporting object containing results of the alignment filtering 
    #out.file <- 
    #  file.path(out.dir, 
    #            paste(out.filename.run3, "_trimmedReads.rda", sep = ""))
    #cat("Exporting object containing trimmed reads : \"", 
    #    out.file, "\".\n", sep = "")
    #save(reads.trimmed, file = out.file)
    
    qual.trimmed <- 
      BStringSet(mcmapply(
        SeqTrimming, cigarOps = cigar.ops, cigarOpLengths = cigar.ops.len, 
        seq = res.sam.filt[["readQuality"]][
          editDistance <= editDist.thold & editDistance > 0], 
        mc.cores = num.cores))
    
    ## Exporting object containing results of the alignment filtering 
    #out.file <- 
    #  file.path(out.dir, 
    #            paste(out.filename.run3, "_trimmedQualityScores.rda", sep = ""))
    #cat("Exporting object containing trimmed quality-strings : \"", 
    #    out.file, "\".\n", sep = "")
    #save(qual.trimmed, file = out.file)
    
    ## **TODO**: MD field is not recomputed for realigned reads, write module
    # NULL as output for OC (stands for original cigar) - tag
    md.ops <- 
      mclapply(strsplit(
        res.sam.filt[["mismatchPOS"]][
          editDistance <= editDist.thold & editDistance > 0], split = ":"),
        function(x) strsplit(x[3], split = "[0-9]+")[[1]], 
        mc.cores = num.cores)
    md.ops.len <- 
      mclapply(strsplit(
        res.sam.filt[["mismatchPOS"]][
          editDistance <= editDist.thold & editDistance > 0], split = ":"),
        function(x) as.numeric(strsplit(x[3], split = "[ACTGN^]+")[[1]]), 
        mc.cores = num.cores)
    mod.len.cumsum <- 
      mclapply(res.sam.filt[["refID"]][
        editDistance <= editDist.thold & editDistance > 0], 
        function(x) cumsum(c(1, as.numeric(head(
          mod.len[strsplit(x, split = ":")[[1]]], n = -1)))), 
        mc.cores = num.cores)
    names.mod.len.cumsum <- 
      mclapply(res.sam.filt[["refID"]][
        editDistance <= editDist.thold & editDistance > 0],
        function(x) names(mod.len[strsplit(x, split = ":")[[1]]]), 
        mc.cores = num.cores)
    
    leftclipped <- as.integer(
      res.sam.filt[["refPOS"]][editDistance <= editDist.thold & editDistance > 0] - 1)
    
    library(Rcpp)
    
    ### Mismatches
    if (mismatch.filter) {
      
      sourceCpp(
        file.path(modseq.dir, 
                  "src/mdSubstitutionAlongReferenceSpaceQualFilter.cpp"))
      
      substitution.modPos <- 
        mcmapply(FUN = mdSubstitutionAlongReferenceSpaceQualFilter, 
                 md_ops_lengths = md.ops.len, md_ops = md.ops, 
                 cigar_ops = cigar.ops, cigar_ops_lengths =  cigar.ops.len,
                 read = as.character(reads.trimmed), 
                 qual = as.character(qual.trimmed), 
                 mod_len_cumsum = mod.len.cumsum, 
                 names_mod_len_cumsum = names.mod.len.cumsum, 
                 offset = leftclipped, 
                 qthold = mismatch.qthold, mc.cores = num.cores)
      
    } else {
      
      sourceCpp(
        file.path(modseq.dir, "src/mdSubstitutionAlongReferenceSpace.cpp"))
      
      substitution.modPos <- 
        mcmapply(FUN = mdSubstitutionAlongReferenceSpace, 
                 md_ops_lengths = md.ops.len, md_ops = md.ops, 
                 cigar_ops = cigar.ops, cigar_ops_lengths =  cigar.ops.len, 
                 read = as.character(reads.trimmed), 
                 mod_len_cumsum = mod.len.cumsum, 
                 names_mod_len_cumsum = names.mod.len.cumsum, 
                 offset = leftclipped[ind.coverage][ind.mapQ][
                   editDistance <= editDist.thold & editDistance > 0], 
                 mc.cores = num.cores)
      
      qualSubstitution.modPos <- 
        mcmapply(FUN = mdSubstitutionAlongReferenceSpace, 
                 md_ops_lengths = md.ops.len, md_ops = md.ops, 
                 cigar_ops = cigar.ops, cigar_ops_lengths =  cigar.ops.len, 
                 read = as.character(qual.trimmed), 
                 mod_len_cumsum = mod.len.cumsum, 
                 names_mod_len_cumsum = names.mod.len.cumsum, 
                 offset = leftclipped,
                 mc.cores = num.cores)
      
    }
    names(substitution.modPos) <- res.sam.filt[["refID"]][
      editDistance <= editDist.thold & editDistance > 0]
    # TODO: verbose level
    # save(substitution.modPos, file = paste(out.dir, out.filename.run3, '_substitutionModPos.rda', sep =""))
    
    aux.mismatch.modDist <- unlist(
      mclapply(strsplit(unlist(substitution.modPos), split = ":"), 
               function(x) x[1],  mc.cores = num.cores)
      )
    mismatch.modDist <- table(
      factor(aux.mismatch.modDist, levels = names(patterns.list))
      )
    
    source(file.path(modseq.dir, "R/functions/PlotMismatchesPerVariant.R"))
    for (i in seq_along(patterns.list)) {
      PlotMismatchesPerVariant(substitution.modPos = substitution.modPos, 
                               ind = c(aux.mismatch.modDist == names(patterns.list)[i]), 
                               mod.len = mod.len, var = names(patterns.list)[i], 
                               tot = aux.modDist[names(patterns.list)[i]], 
                               patterns = patterns.list, 
                               filename = out.filename.run3, 
                               subdirectory = '')
    }
    
    df.mismatch <- data.frame("freq" = mismatch.modDist * 100 / 
                                (aux.modDist * mod.len))
    
    out.file <- 
      file.path(out.dir, 
                paste(out.filename.run3, "_mismatchPerVariant_normalizedFreq_graph.pdf",
                      sep = "")
    pl.mismatch <- ggplot(df.mismatch, aes(x = freq.Var1, y = freq.Freq))
    pl.mismatch <- pl.mismatch + 
      geom_bar(stat = "identity", position = "identity", width = 0.7) +
      labs(x = "Modular variants", y = "Normalized mismatch frequencies [%]") +
      ggtitle("ModSeq | Distribution of mismatches per modular variants") + 
      theme_bw() + 
      theme(plot.title = element_text(face = "bold", size = 16),
            text = element_text(size = 14))
    ggsave(out.file, pl.mismatch, paper = "a4r", width = 12) 
    
    out.file <- 
      file.path(out.dir, 
                paste(out.filename.run3, "_mismatchPerVariant_normalizedFreq_table.csv",
                      sep = "")
    write.csv(as.data.frame(mismatch.modDist), file = out.file)
    
    ### Mismatches per module combination
    mod.dstr <- sort(table(
      res.sam.filt[["refID"]][editDistance <= editDist.thold & editDistance > 0]),
      decreasing = TRUE)
    
    ### Poisson distribution
    source(file.path(modseq.dir, "R/functions/PoissonSNVCall.R"))
    substitution.modComb.poisson <- 
      mcmapply(PoissonSNVCall, mod.dstr[mod.dstr >= min.coverage], 
               names(mod.dstr)[mod.dstr >= min.coverage], 
               MoreArgs = list(variant = substitution.modPos,
                               error.rate = seq.error), mc.cores = num.cores)
    # save(substitution.modComb.poisson, file = paste(out.dir, out.filename.run3, "_mismatchPerModComb.rda", sep = ""))
    
    aux <- unlist(mclapply(substitution.modComb.poisson, names, 
                           mc.cores = num.cores), use.names = FALSE)
    aux.len <- unlist(mclapply(substitution.modComb.poisson, length, 
                               mc.cores = num.cores))
    df.substitution.modComb <- 
      data.frame("moduleCombination" = rep(names(aux.len), aux.len),
                 "mismatch" = aux, 
                 "count" = unlist(substitution.modComb.poisson, 
                                  use.names = FALSE))
    write.csv(df.substitution.modComb, 
              file = paste(out.dir, out.filename.run3, 
                           '_mismatchPerModComb_poisson_table.csv', sep = ""))
    
    substitution.modDist.poisson <- 
      factor(unlist(mclapply(strsplit(unlist(
        mclapply(substitution.modComb.poisson, names), use.names = FALSE), 
        split = ":"), function(x) x[1]), use.names = FALSE), 
        levels = names(patterns.list))
    df.mismatch <- 
      data.frame("freq" = as.numeric(table(substitution.modDist.poisson)) * 100 / 
                   (aux.modDist * mod.len))
    pl.mismatch <- ggplot(df.mismatch, aes(x = freq.Var1, y = freq.Freq))
    pl.mismatch <- pl.mismatch + 
      geom_bar(stat = "identity", position = "identity", width = 0.7) + 
      labs(x = "Modular variants", y = "Normalized mismatch frequencies [%]") +
      ggtitle("ModSeq | Distribution of mismatches per modular variants") + 
      theme_bw() + 
      theme(plot.title = element_text(face = "bold", size = 16),
            text = element_text(size = 14))
    ggsave(paste(out.dir, out.filename.run3, 
                 '_MismatchPerVariant_normalizedFreq_poisson_graph.pdf', 
                 sep = ""), pl.mismatch,  paper = "a4r", width = 12)
    write.csv(as.data.frame(table(substitution.modDist.poisson)), 
              file = paste(out.dir, out.filename.run3, 
                           '_MismatchPerVariant_normalizedFreq_poisson_table.csv',
                           sep = ""))
    
    ### Short indels
    sourceCpp(file.path(modseq.dir, "src/indelsAlongReferenceSpace.cpp"))
    indel.modPos <- 
      mcmapply(FUN = indelsAlongReferenceSpace, cigar_ops = cigar.ops, 
               cigar_ops_lengths = cigar.ops.len, 
               mod_len_cumsum = mod.len.cumsum, 
               names_mod_len_cumsum = names.mod.len.cumsum, 
               offset = leftclipped, 
               mc.cores = num.cores)
    names(indel.modPos) <- res.sam.filt[["refID"]][
      editDistance <= editDist.thold & editDistance > 0]
    # save(indel.modPos, file = paste(out.dir, out.filename.run3, "_indelsModPos.rda", sep = ""))
    
    aux.indel.modDist <- 
      unlist(mclapply(strsplit(unlist(indel.modPos), split = ":"), 
                      function(x) x[2], mc.cores = num.cores))
    aux.indel.type <- 
      unlist(mclapply(strsplit(unlist(indel.modPos), split = ":"),
                      function(x) x[1], mc.cores = num.cores))
    indel.modDist <- table(factor(aux.indel.modDist, levels = names(patterns.list)))
    df.indel <- 
      data.frame("mod" = factor(aux.indel.modDist, levels = names(patterns.list)),
                 "type" = aux.indel.type)
    df.indel <- as.data.frame(table(df.indel))  # To normalize it
    df.indel[df.indel$Freq == 0,]$Freq <- N
    pl.indel <- ggplot(df.indel, aes(x=mod, y = as.numeric(Freq) * 100 / 
                                       (rep(aux.modDist, 2) * rep(mod.len, 2)),
                                     fill = type))
    pl.indel <- pl.indel + 
      geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
      labs(x = "Modular variants", y = "Normalized indel frequencies [%]") +
      scale_fill_discrete(breaks = c("Del", "In"), 
                          labels = c("Deletion", "Insertion"), 
                          guide = guide_legend(reverse = TRUE)) + 
      ggtitle("ModSeq | Distribution of indels per modular variants") + 
      theme_bw() + 
      theme(legend.position = "top", 
            plot.title = element_text(face = "bold", size = 16),
            text = element_text(size = 14))
    ggsave(paste(out.dir, out.filename.run3, 
                 '_IndelPerVariant_normalizedFreq_graph.pdf', sep = ""), 
           pl.indel, paper = "a4r", width = 12)
    write.csv(df.indel, file = paste(out.dir, out.filename.run3, 
                                     '_IndelPerVariant_normalizedFreq_table.csv',
                                     sep = ""))
    
    ### Indels per module combination
    indel.modComb.poisson <- 
      mcmapply(PoissonVarCall, mod.dstr[mod.dstr >= min.coverage], 
               names(mod.dstr)[mod.dstr >= min.coverage], 
               MoreArgs = list(indels = indel.modPos, error.rate = seq.error),
               mc.cores = num.cores)
    # save(indel.modComb.poisson, file = paste(out.dir, out.filename.run3, "_indelsPerModComb.rda", sep = ""))
    
    aux <- unlist(mclapply(indel.modComb.poisson, names, mc.cores = num.cores),
                  use.names = FALSE)
    aux.len <- unlist(mclapply(indel.modComb.poisson, length, 
                               mc.cores = num.cores))
    df.indel.modComb <- 
      data.frame("moduleCombination" = rep(names(aux.len), aux.len),
                 "indel" = aux, 
                 "count" = unlist(indel.modComb.poisson, use.names = FALSE))
    write.csv(df.indel.modComb, 
              file = paste(out.dir, out.filename.run3, 
                           '_IndelsPerModComb_poisson_table.csv', sep = ""))
    
    indel.modDist.poisson <- 
      factor(unlist(mclapply(strsplit(unlist(mclapply(indel.modComb.poisson, names)),
                                      split = ":"), function(x) x[2])),
             levels = names(patterns.list))
    aux.indel.type <- 
      unlist(mclapply(strsplit(unlist(mclapply(indel.modComb.poisson, names)), 
                               split = ":"), function(x) x[1], 
                      mc.cores = num.cores))
    df.indel <- data.frame("mod" = indel.modDist.poisson, 
                           "type" = aux.indel.type)
    df.indel <- as.data.frame(table(df.indel))  
    df.indel[df.indel$Freq == 0,]$Freq <- NA
    pl.indel<- ggplot(df.indel, aes(x = mod, y = as.numeric(Freq) * 100 /
                                      (rep(aux.modDist, 2)*rep(mod.len, 2)),
                                    fill = type))
    pl.indel <- pl.indel + 
      geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
      labs(x = "Modular variants", y = "Normalized indel frequencies [%]") +
      scale_fill_discrete(breaks = c("Del", "In"), 
                          labels = c("Deletion", "Insertion"), 
                          guide = guide_legend(reverse = TRUE)) + 
      ggtitle("ModSeq | Distribution of indels per modular variants") + 
      theme_bw() + 
      theme(legend.position = "top", 
            plot.title = element_text(face = "bold", size = 16),
            text = element_text(size = 14))
    ggsave(paste(out.dir, out.filename.run3, 
                 '_IndelsPerVariant_normalizedFreq_poisson_graph.pdf', 
                 sep = ""), pl.indel,  paper = "a4r", width = 12)
    write.csv(df.indel, 
              file = paste(out.dir, out.filename.run3, 
                           '_IndelPerVariant_normalizedFreq_poisson_table.csv',
                           sep = ""))
  }
  
  if (exists("patterns")) {
    rm(patterns)
  }
}


###############################################################################
cat("Runtime:", round(as.numeric(difftime(Sys.time(), Ti, units = "mins")), 
                      digits=4), "min\n")
cat("Done!\n")
sink()

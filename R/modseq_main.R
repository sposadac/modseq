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
# **TODO** add R directory to search path
# modseq.dir Path to modseq head directory. At the moment, set to current directory
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
  
  source(file.path(modseq.dir, 'R/functions/PlotQualityDistribution.R'))
  source(file.path(modseq.dir, 'R/functions/PlotReadLengthDistribution.R'))
  
  cat("Sequencing mode: \"", seq.mode, "\".\n", sep = "")
  if (seq.mode == "PE") {
    
    readF1 <- readFastq(dirPath = in.seqDir, pattern = forward.filename)
    readF2 <- readFastq(dirPath = in.seqDir, pattern = reverse.filename)
    len.rawData <- length(readF1)
    # checks
    if (length(readF2) != len.rawData) {
      stop("Number of reads in input files does not match: \n Forward set: ",
           len.rawData, " reads, \n Reverse set: ", length(readF2), " reads.\n")
    }
    cat("Number of reads: ", len.rawData, " x 2 (mode: PE).\n", sep = "")
    PlotQualityDistribution(readF1, forward.filename, out.dir)
    PlotQualityDistribution(readF2, reverse.filename, out.dir)
    PlotReadLengthDistribution(readF1, forward.filename, out.dir)
    PlotReadLengthDistribution(readF2, reverse.filename, out.dir)
    
  } else if (seq.mode == "SE") {
    
    readF1 <- readFastq(dirPath = in.seqDir, pattern = in.filename)
    len.rawData <- length(readF1)
    cat("Number of reads: ", len.rawData, ".\n", sep = "")
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
  
  source(file.path(modseq.dir, "R/functions/IlluminaStat.R"))
  
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
      wdir = wdir)
  
  } else if (seq.mode == "PE") {
    
    reads.nQtrim <- Run1_qualTrimming(
      readF1 = readF1, in.filename = forward.filename, seq.mode = seq.mode
      qtrim.3end = qtrim.3end, qtrim.thold = qtrim.thold, 
      out.filename = out.filename.run1, out.ssplot = out.ssplot,
      readF2 = readF2, reverse.filename = reverse.filename, wdir = wdir)
    readF1.nQtrim <- reads.nQtrim[[1]]
    readF2.nQtrim <- reads.nQtrim[[2]]
    
  }

  qtrim.flag <- as.integer(1) 
  cat("Object \'qtrim.flag\' is set to ", qtrim.flag,".\n", sep = "")
  
  keep <- append(keep, c("len.rawData"))
  if (run[2] == 1) {
    keep <- append(keep, c("readF1.nQtrim", "IlluminaStat"))
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
  if (!existsFunction("IlluminaStat")) {
    source(file.path(modseq.dir, "R/functions/IlluminaStat.R"))
  }

  if (mem.trace) {
    memTrace <- c(memTrace, MemTrace())
  }
  
  ## PE read assembly
  source(file.path(modseq.dir, "R/functions/Run2_peAssembly.R"))
  if (qtrim.flag == 0) {
    
    Run2_peAssembly(readF1 = readF1, readF2 = readF2, qtrim.flag = qtrim.flag, 
                    forward.file = forward.filename, reverse.file = reverse.filename, 
                    out.ssplot = out.ssplot, out.filename.run1 = out.filename.run1, 
                    out.filename.run2 = out.filename.run2, in.seqDir = in.seqDir,
                    len.rawData = len.rawData, wdir = wdir)
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
    if (!exists("len.rawData")) {
      len.rawData <- NULL
    }
    Run2_peAssembly(readF1 = readF1.nQtrim, readF2 = readF2.nQtrim, 
                    qtrim.flag = qtrim.flag, forward.file = in.tr1, 
                    reverse.file = in.tr2, out.ssplot = out.ssplot, 
                    out.filename.run1 = out.filename.run1, 
                    out.filename.run2 = out.filename.run2, 
                    in.seqDir = in.seqDir, 
                    len.rawData = len.rawData, wdir = wdir)
    rm(list = c("readF1.nQtrim", "readF2.nQtrim"))
    
  }
  paired.flag <- as.integer(1)
  cat("Object \'paired.flag\' is set to ", paired.flag, ". \n", sep = "")
  
  if (mem.trace) {
    memTrace <- c(memTrace, MemTrace())
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
  ## Loading module-table
  source(file.path(modseq.dir, "R/functions/LoadModuleTable.R"))
  patterns <- 
    LoadModuleTable(in.modulesDir = in.modDir, modules.filename = mod.filename)
  
  source((file.path(modseq.dir, "R/functions/GetInputFile.R")))
  in.file <- GetInputFile(wdir)
  if (file.exists(in.file[[1]])) {
    reads.subj <- readDNAStringSet(
      file = in.file[[1]], format = in.file[[2]], use.names = TRUE) 
    cat("Input reads file: \"", in.file[[1]], "\".\n", sep = "")
    
  } else {
    stop("File \"", in.file, "\" not found.\n")
    
  }
  reads.subj.len <- length(reads.subj)
  
  if (map.mode == "gls") {
    
    source(file.path(modseq.dir, "R/functions/Run3_gls.R"))
    retList <- Run3_gls(
      patterns=patterns, reads=reads.subj, num.reads=reads.subj.len, in.modDir, 
      mod.filename, res.listName, res.counts.filename, out.dir, 
      gls.ambiguity = gls.ambiguity, gls.direction = gls.direction, 
      gls.mma = gls.mma, mem.trace = mem.trace, memTrace = memTrace, 
      run.info = run.info, num.cores = num.cores)
    
    res.list <- retList[[1]]
    if (length(retList) > 1) {
      mod.comb <- retList[[2]]
      res.list.lengths <- retList[[3]]
      res.list.name <- retList[[4]]
      if (run[4] == 1) {
        keep <- append(keep, c("mod.comb", "res.list.lengths", "res.list.name")) 
      }
    } 
    
    if (run[4] == 1) {
      keep <- append(keep, c("reads.subj", "res.list")) #TODO read.subj outside gls?, and patterns?
    }
    
  } else if (map.mode == "bwa") {
    
    reads.file <- GetInputFile(wdir)
    reads.file <- reads.file [[1]]
    if (!file.exists(reads.file)) {
      stop ("File \"", reads.file, "\" not found.\n")
    }
    
    in.file <- file.path(in.modDir, paste(mod.filename, '_modComb.fasta', sep = ""))
    if (!file.exists(in.file)) {
      source(paste(modseq.dir, "R/functions/ModuleCombinationsGen.R", sep = ""))
      ModuleCombinationsGen(mod.filename, pattern = patterns, 
                            num.cores = num.cores)
      if (run[4] == 1) {
        keep <- append(keep, c("ModuleCombinationsGen"))
      }
    }
    # Indexing reference sequences -- module-combinations
    if (!file.exists(paste(in.modDir, mod.filename, "_modComb.fasta.bwt", 
                           sep = ""))) {
      run.index <- paste(bwa.path, 'bwa index ', in.modDir, mod.filename, 
                         "_modComb.fasta", sep = "")
      cat("Building index reference sequences ... \n")
      system(run.index)
    }
    # Mcomp: stands for Picard compatibility (option M)
    run.bwa <- paste(bwa.path, 'bwa mem -M -t ', num.cores,' -c ', bwa.cVal, 
                     ' ', in.modDir, mod.filename, '_modComb.fasta ', 
                     reads.file, ' > ', out.dir, out.filename.run2, 
                     '_bwaMEM_Mcomp.sam 2> ', out.dir,'out_', out.filename.run2, 
                     '_bwaMEM_Mcomp.txt', sep = "")
    cat("BWA - Starting read mapping ...\n")
    ti.search <- Sys.time()
    system(run.bwa)   
    cat("Runtime:", round(as.numeric(
      difftime(Sys.time(), ti.search, units = "mins")), digits = 4), "min.\n")
    cat("Mapping done!\n")
    
    cat("Local realignment - Preparing \"", in.modDir, mod.filename,
        "_modComb.fasta\" to use as reference ...\n", sep = "")
    # Dictionary - reference sequences
    if (!file.exists(paste(in.modDir, mod.filename, '_modComb.dict', 
                           sep = ""))) {
      system(paste('java -jar ../gatk/picard-tools-1.129/picard.jar ', 
                   'CreateSequenceDictionary R=', in.modDir, mod.filename, 
                   '_modComb.fasta O=', in.modDir, mod.filename, 
                   '_modComb.dict', sep = ""))
    }
    # Indexing - reference sequences
    if (!file.exists(paste(in.modDir, mod.filename, '_modComb.fasta.fai', 
                           sep=""))) {
      system(paste('samtools faidx ', in.modDir, mod.filename, '_modComb.fasta', 
                   sep = ""))
    } 
    
    cat("Local realignment - Sorting \"", out.filename.run2, 
        '_bwaMEM.bam',"\" file... \n", sep = "")
    # Convert SAM to BAM
    if (!file.exists(paste(out.dir, out.filename.run2, '_bwaMEM.bam', 
                           sep = ""))) {
      system(paste('samtools view -b -S -o ', out.dir, out.filename.run2, 
                   '_bwaMEM.bam ', out.dir, out.filename.run2, 
                   '_bwaMEM_Mcomp.sam', sep = ""))
    }
    # BAM - sorted by coordinates
    if (!file.exists(paste(out.dir, out.filename.run2, '_bwaMEM_', 
                           'sortedPicard.bam', sep = ""))) {
      system(paste('java -jar ../gatk/picard-tools-1.129/picard.jar SortSam ', 
                   'INPUT=', out.dir, out.filename.run2, '_bwaMEM.bam ', 
                   'OUTPUT=', out.dir, out.filename.run2, '_bwaMEM_', 
                   'sortedPicard.bam SORT_ORDER=coordinate', sep = ""))
    }
    # Mark read duplicates
    # Duplicated could bias variant detection by adding excessive coverage 
    # depth at a variant locus
    if (bwa.dupl) {
      cat("Local realignment - Mark read duplicates ... \n")
      system(paste('java -jar ../gatk/picard-tools-1.129/picard.jar ', 
                   'MarkDuplicates INPUT=', out.dir, out.filename.run2, 
                   '_bwaMEM_sortedPicard.bam OUTPUT=', out.dir, 
                   out.filename.run3, '.bam METRICS_FILE=', out.dir, 
                   out.filename.run3, '_metrics.txt', sep = ""))
      # Add or replace read groups
      system(paste('java -jar ../gatk/picard-tools-1.129/picard.jar ', 
                   'AddOrReplaceReadGroups INPUT=', out.dir, out.filename.run3, 
                   '.bam OUTPUT=', out.dir, out.filename.run3, 
                   '_readGroup.bam RGID=group1 RGLB=lib1 RGPL=illumina ', 
                   'RGPU=unit1 RGSM=sample1', sep = ""))
    } else {
      # Add or replace read groups
      system(paste('java -jar ../gatk/picard-tools-1.129/picard.jar ', 
                   'AddOrReplaceReadGroups INPUT=', out.dir, out.filename.run3, 
                   '_sortedPicard.bam OUTPUT=', out.dir, out.filename.run3, 
                   '_readGroup.bam RGID=group1 RGLB=lib1 RGPL=illumina ', 
                   'RGPU=unit1 RGSM=sample1', sep = ""))
    }
    
    # Index bam file
    system(paste('java -jar ../gatk/picard-tools-1.129/picard.jar ', 
                 'BuildBamIndex INPUT=', out.dir, out.filename.run3, 
                 '_readGroup.bam', sep = ""))
    cat("Local realignment - identify intervals for local realignment. \n")
    # Intervals for realignment
    system(paste('java -jar ../gatk/GenomeAnalysisTK.jar -T ', 
                 'RealignerTargetCreator -R ', in.modDir, mod.filename, 
                 '_modComb.fasta -I ', out.dir, out.filename.run3, 
                 '_readGroup.bam -o ', out.dir, out.filename.run3, 
                 '_forIndelRealigner.intervals', sep = ""))
    cat("Local realignment. \n")
    # Realignment
    system(paste('java -jar ../gatk/GenomeAnalysisTK.jar -T IndelRealigner -R ',
                 in.modDir, mod.filename, '_modComb.fasta -I ', out.dir,
                 out.filename.run3, '_readGroup.bam -targetIntervals ',
                 out.dir, out.filename.run3, '_forIndelRealigner.intervals -o ', 
                 out.dir, out.filename.run3, '_realigned.bam', 
                 sep = ""))  
    # Convert BAM to SAM
    system(paste('samtools view ', out.dir, out.filename.run3, 
                 '_realigned.bam > ', out.dir, out.filename.run3, 
                 '_realigned.sam', sep = ""))
    
    if (run[4] == 1) {
      if (bwa.dupl) {
        res.sam.realn <- scan(file = paste(out.dir, out.filename.run3, 
                                           "_realigned.sam", sep = ""), 
                              what = list(character(), integer(), character(), 
                                          integer(), integer(), character(), 
                                          character(), integer(), integer(), 
                                          character(), character(), character(),
                                          character(), character(), character(),
                                          character(), character()), skip = 0,
                              flush = TRUE, fill = TRUE, quote = "")
        names(res.sam.realn) <- c("readID", "flag", "refID", "refPOS", "mapQ", 
                                  "CIGAR", "", "", "", "readSeq", "readQuality",
                                  "mismatchPOS", "", "", "editDistance", "", "")
      } else {
        res.sam.realn <- scan(file = paste(out.dir, out.filename.run3, 
                                           "_realigned.sam", sep = ""), 
                              what = list(character(), integer(), character(), 
                                          integer(), integer(), character(), 
                                          character(), integer(), integer(), 
                                          character(), character(), character(),
                                          character(), character(), character(),
                                          character()), skip = 0, 
                              flush = TRUE, fill = TRUE, quote = "") 
        names(res.sam.realn) <- c("readID", "flag", "refID", "refPOS", "mapQ", 
                                  "CIGAR", "", "", "", "readSeq", "readQuality",
                                  "mismatchPOS", "", "editDistance", "", "")
      }
      keep <- append(keep, c("res.sam.realn"))
    } else {
      res.sam.realn <- scan(file = paste(out.dir, out.filename.run3, 
                                         "_realigned.sam", sep = ""), 
                            what = list(character(), integer()), skip = 0, 
                            flush = TRUE, fill = TRUE, quote = "") 
      names(res.sam.realn) <- c("readID", "flag")
    }
    reads.dupl <- 
      sum(res.sam.realn[["flag"]] == 1024 | res.sam.realn[["flag"]] == 1040)
    reads.map <- 
      sum(res.sam.realn[["flag"]] == 0 | res.sam.realn[["flag"]] == 16) + 
      reads.dupl
    cat("Found ", reads.map, " module-combinations in ", reads.subj.len, 
        " reads (", round(reads.map * 100 / reads.subj.len, digits = 2), 
        "%). \n", sep = "")
    cat("Found", reads.dupl, " marked reads as duplicates (", 
        round(reads.dupl * 100 / reads.subj.len, digits = 2), "%). \n", 
        sep = "") 
  }
  
  if (run[4] == 1 || run[5] == 1) {
    keep <- append(keep, c("patterns", "num.cores"))
  }
  
  if (mem.trace) {
    memTrace <- c(memTrace, MemTrace())
  }
  rm(list=ls()[!(ls() %in% keep)]) 
}

###############################################################################
cat("Runtime:", round(as.numeric(difftime(Sys.time(), Ti, units = "mins")), 
                      digits=4), "min\n")
cat("Done!")
sink()

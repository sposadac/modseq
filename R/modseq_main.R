### MAIN SCRIPT ###
# Authors: 
# Susana Posada Cespedes <susana.posada@bsse.ethz.ch>
# Steven Schmitt
####

### DELETE WORKSPACE 
rm(list=ls())

## Get current work directory
wdir <- getwd()

### LOAD USER-DEFINED OPTIONS AND PACKAGES
# **TODO** add R directory to search path
source(file.path(wdir, "R/modseq_input.R"))
source(file.path(wdir, "R/functions/InputCheck.R")) 
source(file.path(wdir, "R/functions/NameGen.R"))
names(run) <- c(1, 2, 3, 4, 5)

### TRACKING TOTAL MEMORY USAGE
# MemTrace() Memory in Megabytes
if (mem.trace) {
  source(file.path(wdir, "R/functions/MemTrace.R"))
  memTrace <- MemTrace()
}

### START RUN
Ti <- Sys.time()

# Log-file
log.file <- file.path(out.dir, 
                      paste("/out_", out.filename.run2, "_run",
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
  
  source(file.path('R/functions/plotQualityDistribution.R'))
  source(file.path('R/functions/plotReadLengthDistribution.R'))
  
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
    plotQualityDistribution(readF1, forward.filename, out.dir)
    plotQualityDistribution(readF2, reverse.filename, out.dir)
    plotReadLengthDistribution(readF1, forward.filename, out.dir)
    plotReadLengthDistribution(readF2, reverse.filename, out.dir)
  } else if (seq.mode == "SE") {
    readF1 <- readFastq(dirPath = in.seqDir, pattern = in.filename)
    len.rawData <- length(readF1)
    cat("Number of reads: ", len.rawData, ".\n", sep = "")
    plotQualityDistribution(readF1, in.filename, out.dir)
    plotReadLengthDistribution(readF1, in.filename, out.dir)
  }
}

###############################################################################
## 1: QUALITY TRIMMING
###############################################################################
if (run[1] == 1) {  
  cat("====================================================================\n")  
  cat("Quality trimming. \n")
  cat("====================================================================\n")  
  
  source(file.path(wdir, 'R/functions/IlluminaStat.R'))
  
  ### QUALITY TRIMMING: 
  # Uncalled bases (i.e. N) are removed from either both ends (qtrim.3end = 0)
  # or only the right end (qtrim.3end = 1).
  # Bases are removed if quality score is less than or equal to qtrim.thold.
  source(file.path(wdir, 'R/functions/Run1_qualTrimming.R'))
  if (seq.mode == "SE") {
    readF1.nQtrim <- 
      Run1_qualTrimming(readF1, in.filename, wdir, out.filename.run1, 
                        qtrim.3end, qtrim.thold, seq.mode, out.ssplot)
  } else if (seq.mode == "PE") {
    reads.nQtrim <- 
      Run1_qualTrimming(readF1, forward.filename, wdir, out.filename.run1,
                        qtrim.3end, qtrim.thold, seq.mode, out.ssplot,
                        readF2, reverse.filename)
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
  if (!existsFunction('IlluminaStat')) {
    source(file.path(wdir, 'R/functions/IlluminaStat.R'))
  }

  if (mem.trace) {
    memTrace <- c(memTrace, MemTrace())
  }
  
  ## PE read assembly
  source(file.path(wdir, 'R/functions/Run2_peAssembly.R'))
  if (qtrim.flag == 0) {
    Run2_peAssembly(readF1, readF2, wdir, out.filename.run2, qtrim.flag, 
                    out.ssplot, forward.filename, reverse.filename, in.seqDir,
                    len.rawData)
    rm(list = c("readF1", "readF2"))
    
  } else if (qtrim.flag == 1) {

    in.tr1 <- file.path(wdir, 'data/processed', 
                        paste(out.filename.run1, "_1_nQtrim.fastq", sep = ""))
    in.tr2 <- file.path(wdir, 'data/processed',
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
    Run2_peAssembly(readF1.nQtrim, readF2.nQtrim, wdir, out.filename.run2, 
                    qtrim.flag, out.ssplot, in.tr1, in.tr2, in.seqDir, 
                    len.rawData)
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
  
  ## Loading module-table
  cat("Module-table used: \"", in.modDir, mod.filename, ".csv\". \n", sep = "")
  patterns <- 
    data.frame(
      read.csv(file = file.path(in.modDir, paste(mod.filename, ".csv", sep = "")),
               header = TRUE, sep = ",", dec = "."), stringsAsFactors = FALSE
      )
  num.cores <- detectCores()
  
  if (paired.flag == 1) {  
    
    ## Paired reads 
    cat("Read mapping of paired reads")
    in.file <- file.path(wdir, "data/paired/fasta", 
                         paste(out.filename.run2, "_PANDAseq.fasta", sep = ""))
    if (file.exists(in.file)) {
      # Must be format = "fasta", because when a FASTQ file is given as input, all
      # the sequences must have the same length
      reads.subj <- 
        readDNAStringSet(file = in.file, format = "fasta", use.names = TRUE)
      cat("PE-reads file: ", in.file, ".\n")
    } else {
      stop("File \"", in.file, "\" not found. \n")
    }
    
  } else if (paired.flag == 0) {  
    
    ## Unpaired reads
    if (qtrim.flag == 0) {  
      
      ## Unpaired and untrimmed reads - raw reads 
      if (seq.mode == "SE") {
        
        cat("Read mapping of untrimmed and non-paired reads")
        in.file <- file.path(in.seqDir, paste(in.filename, ".fastq", sep = ""))
        cat("Sequencing-reads file: ", in.file, ".\n")
        
      } else if (seq.mode == "PE") {
        
        cat("Read mapping of untrimmed and unpaired reads")
        if (paired.file == "f") {
          in.file <- file.path(in.seqDir,  
                               paste(forward.filename, ".fastq", sep = ""))
          cat("Forward-reads file: ", in.file, ".\n")
        } else if (paired.file == "r") {
          in.file <- file.path(in.seqDir,  
                               paste(reverse.filename, ".fastq", sep = ""))
          cat("Reverse-reads file: ", in.file, ".\n")
        }
        
      }
      reads.subj <- 
        readDNAStringSet(file = in.file, format = "fastq", use.names = TRUE)
      
    } else {   
      
      ## Unpaired, but trimmed reads
      cat("Read mapping of trimmed but non-paired reads")
      if (seq.mode == "SE" || paired.file == "f") { 
        in.file <- file.path(wdir, 'data/processed/fasta', 
                            paste(out.filename.run1, "_1_nQtrim.fasta", sep = ""))
      } else if (paired.file == "r") {
        in.file <- file.path(wdir, 'data/processed/fasta', 
                             paste(out.filename.run1, "_2_nQtrim.fasta", sep = ""))
      }
      if (file.exists(in.file)) {
        reads.subj <- 
          readDNAStringSet(file = in.file, format = "fasta", use.names = TRUE) 
        cat("Input reads file: ", in.file, ".\n")
      } else {
        stop("File \"", in.file, "\" not found.\n")
      }
      
    }
  }
  reads.subj.len <- length(reads.subj)
  
  if (map.mode == "gls") {
    
    source(file.path(wdir, 'R/functions/Search_vmatchV2.R'))
    if (run[4] == 1) {
      keep <- append(keep, c("mod.tot", "reads.subj", "reads.subj.len", 
                             "res.list"))
    }
    
  } else if (map.mode == "bwa") {
    
    if (paired.flag == 1) {  # Paired reads
      
      reads.file <- 
        file.path(wdir, "data/paired", 
                  paste(out.filename.run2, "_PANDAseq.fastq", sep = ""))
      
    } else if (paired.flag == 0) { 
      
      if (qtrim.flag == 0) { # Unpaired and untrimmed files - raw files
        if (seq.mode == "SE") {
          reads.file <- 
            file.path(in.seqDir, paste(in.filename, ".fastq", sep = ""))
        } else if (seq.mode == "PE") {
          if (paired.file == "f") {
            reads.file <- 
              file.path(in.seqDir, paste(forward.filename, ".fastq", sep = ""))
          } else if (paired.file == "r") {
            reads.file <- 
              file.path(in.seqDir, paste(reverse.filename, ".fastq", sep = ""))
          }
        }
      } else { # Unpaired, but trimmed files 
        ###****HERE****
        if (seq.mode == "SE" || paired.file == "f") {
          reads.file <- paste("./DATA/PROCESSED/", out.filename.run1, 
                              "_1_nQtrim.fastq", sep = "") 
        } else if (paired.file == "r") {
          reads.file <- paste("./DATA/PROCESSED/", out.filename.run1,
                              "_2_nQtrim.fastq", sep = "")
        }
      }
      
    }
    if (!file.exists(reads.file)) {
      stop ("File \"", reads.file, "\" not found.\n")
    }
    
    if (!file.exists(paste(in.modDir, mod.filename, '_modComb.fasta', 
                           sep = ""))) {
      source(paste(wdir, '/functions/ModuleCombinationsGen.R', sep = ""))
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

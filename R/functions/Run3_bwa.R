Run3_bwa <- function(reads.file, in.modDir, mod.filename, bwa.path, bwa.cVal,
                     bwa.dupl, gatk.path, picard.path, samtools.path, 
                     out.filename.run2, out.filename.run3, patterns=NULL, 
                     modseq.dir=NULL, out.dir=NULL, num.cores=numeric(0)) {
  
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
  
  mod.file <- file.path(in.modDir, 
                        paste(mod.filename, '_modComb.fasta', sep = ""))
  
  if (!file.exists(mod.file)) {  
    if (!existsFunction("ModuleCombinationsGen")) {
      source(file.path(modseq.dir, "R/functions/ModuleCombinationsGen.R"))
    }
    
    if (!is.data.frame(patterns)) {
      if (!existsFunction("LoadModuleTable")) {
        source(file.path(modseq.dir, "R/functions/LoadModuleTable.R"))
      }
      patterns <- 
        LoadModuleTable(in.modulesDir = in.modDir, 
                        modules.filename = mod.filename) 
    }
    
    ModuleCombinationsGen(modules.filename = mod.filename, pattern = patterns, 
                          in.modDir = in.modDir, modseq.dir = modseq.dir, 
                          num.cores = num.cores)
    
  }
  
  ## Setting bwa variable
  if (length(bwa.path) == 0) {
    bwa <- "bwa"
  } else {
    bwa <- file.path(bwa.path, "bwa")
  }
  
  ## Indexing reference sequences -- module-combinations sequences
  out.file <- file.path(in.modDir, 
                        paste(mod.filename, "_modComb.fasta.bwt", sep = ""))
  if (!file.exists(out.file)) {
    run.index <- paste(bwa, "index", mod.file)
    
    cat("Building index reference sequences ... \n")
    system(run.index)
    cat("Output file: \"", out.file,"\".\n", sep = "")
    
  }
  
  ## Running bwa mem
  # Mcomp: stands for Picard compatibility (option M)
  sam.file <- 
    file.path(out.dir, 
              paste(out.filename.run2, "_bwaMEM_Mcomp.sam", sep = ""))
  out.file <- 
    file.path(out.dir, 
              paste("out_", out.filename.run2, "_bwaMEM_Mcomp.txt", sep = ""))
  
  run.bwa <- paste(bwa, "mem -M -t", num.cores, "-c", bwa.cVal, mod.file, 
                   reads.file, ">", sam.file, "2>", out.file)
  
  cat("BWA - Starting read mapping ...\n")
  ti.search <- Sys.time()
  system(run.bwa)   
  cat("Runtime:", round(as.numeric(
    difftime(Sys.time(), ti.search, units = "mins")), digits = 4), "min.\n")
  cat("Mapping done!\n")
  cat("Output files: \n")
  cat("Alignment file (sam format): \"", sam.file, "\". \n", sep = "")
  cat("Log file: \"", out.file, "\". \n", sep = "")
  
  ## Preparing dictionary - reference sequences
  ## Setting picard variable
  if (length(picard.path) == 0) {
    picard <- "picard.jar"
  } else {
    picard <- file.path(picard.path, "picard.jar")
  }
  cat("Local realignment - Preparing \"", mod.file, "\" to use as reference", 
      " ...\n", sep = "")
  out.file <- file.path(in.modDir, 
                        paste(mod.filename, '_modComb.dict', sep = ""))
  
  if (!file.exists(out.file)) {
    cat("Local realignment - Creating sequence dictionary for the reference", 
        "sequences ...\n")
    system(paste("java -jar ", picard, " CreateSequenceDictionary R=", 
                 mod.file,  " O=", out.file, sep = ""))
    cat("Output file: \"", out.file,"\".\n", sep = "")
  }
  
  ## Indexing - reference sequences
  ## Setting samtools variable
  if (length(samtools.path) == 0) {
    samtools <- "samtools"
  } else {
    samtools <- file.path(samtools.path, "samtools")
  }
  out.file <- file.path(in.modDir,
                        paste(mod.filename, "_modComb.fasta.fai", sep = ""))
  
  if (!file.exists(out.file)) {
    cat("Local realignment - Indexing reference sequences using samtools", 
        "...\n")
    system(paste(samtools, "faidx", mod.file))
    cat("Output file: \"", out.file,"\".\n", sep = "")
  } 
  
  ## Converting SAM to BAM and sorting
  bam.file <- 
    file.path(out.dir, 
              paste(out.filename.run2, "_bwaMEM_Mcomp.bam", sep = ""))
  
  if (!file.exists(bam.file)) {
    cat("Local realignment - Converting sam to bam file ... \n")
    system(paste(samtools, "view -b -S -o", bam.file, sam.file))
    cat("Alignment file (bam format): \"", bam.file,"\".\n", sep = "")
  }
  
  ## BAM - sorted by coordinates
  bamSort.file <- 
    file.path(out.dir,
              paste(out.filename.run2, "_bwaMEM_sortedPicard.bam", sep = ""))
  
  if (!file.exists(bamSort.file)) {
    cat("Local realignment - Sorting bam file by coordinates .. \n")
    system(paste("java -jar ", picard, " SortSam INPUT=", bam.file, 
                 " OUTPUT=", bamSort.file, " SORT_ORDER=coordinate", sep = ""))
    cat("Output file: \"", bamSort.file,"\".\n", sep = "")
  }
  
  ## Mark read duplicates
  # Duplicated can introduce bias in variant detection by adding excessive 
  # coverage depth at a variant locus
  if (bwa.dupl) {
    
    bam.file <- file.path(out.dir, paste(out.filename.run3, ".bam", sep = ""))
    out.file <- file.path(out.dir, 
                          paste(out.filename.run3, "_metrics.txt", sep = ""))
    cat("Local realignment - Mark read duplicates ... \n")
    system(paste("java -jar ", picard, " MarkDuplicates INPUT=", bamSort.file,
                 " OUTPUT=", bam.file, " METRICS_FILE=", out.file, sep = ""))
    cat("Output files: \n")
    cat("Alignment file (bam format): \"", bam.file, "\".\n", sep = "")
    cat("Metrics file: \"", out.file, "\".\n", sep = "")
    
  } else {
    
    bam.file <- bamSort.file
    
  }
  ## Add or replace read groups
  out.file <- file.path(out.dir, 
                        paste(out.filename.run3, "_readGroup.bam", sep = ""))
  system(paste("java -jar ", picard, " AddOrReplaceReadGroups INPUT=", 
               bam.file, " OUTPUT=", out.file, " RGID=group1 RGLB=lib1 ", 
               "RGPL=illumina RGPU=unit1 RGSM=sample1", sep = ""))
  cat("Output file: \"", out.file,"\".\n", sep = "")
  
  ## Indexing bam file
  system(
    paste("java -jar ", picard, " BuildBamIndex INPUT=", out.file, sep = ""))
  
  
  ## Intervals for realignment
  ## Setting gatk variable
  if (length(gatk.path) == 0) {
    gatk <- "GenomeAnalysisTK.jar"
  } else {
    gatk <- file.path(gatk.path, "GenomeAnalysisTK.jar")
  }
  cat("Local realignment - identifying intervals for local realignment. \n")
  bam.file <- out.file
  out.file <- 
    file.path(out.dir, 
              paste(out.filename.run3, "_forIndelRealigner.intervals", 
                    sep = ""))
  system(paste("java -jar ", gatk, " -T RealignerTargetCreator -R ", mod.file, 
               " -I ", bam.file, " -o ", out.file, sep = ""))
  cat("Output file: \"", out.file,"\".\n", sep = "")
  
  ## Realignment
  cat("Local realignment. \n")
  realign.file <- 
    file.path(out.dir, 
              paste(out.filename.run3, "_realigned.bam", sep = ""))
  system(paste("java -jar ", gatk, " -T IndelRealigner -R ", mod.file, " -I ",
               bam.file, " -targetIntervals ", out.file, " -o ", realign.file,                 
               sep = ""))  
  cat("Output file: \"", realign.file,"\".\n", sep = "")
  
  ## Converting BAM to SAM
  sam.file <- 
    file.path(out.dir, 
              paste(out.filename.run3, "_realigned.sam", sep = ""))
  cat("Local realignment - Converting bam to sam file ... \n")
  system(paste(samtools, "view", realign.file, ">", sam.file))
  cat("Alignment file (sam format): \"", sam.file,"\".\n", sep = "")
  
  return(sam.file)
}
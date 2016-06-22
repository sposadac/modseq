# Library ShortRead required
SubsamplingFastq <- function(in.seqDir, in.filename, reverse.filename=NULL, num=1e4,
                             seed=47, wdir=NULL) {
  
  ## Function arguments
  # in.filename         Name of the file containing input reads (SE-mode) or 
  #                     forward reads (PE-mode).
  # reverse.filename    (optional) Name of the file containing reverse reads, 
  #                     expected in PE-mode. 
  # num                 (optional) Number of reads to keep in the sub-sample, by
  #                     default 1e4.
  # seed                (optional) Seed to ensure reproducible sampling, by 
  #                     default 47.
  # wdir                (optional) Directory to search for input files, by 
  #                     default current working directory.
  
  ## Whenever input directory is not specified, assumed to be the current
  #  working directory
  if (is.null(wdir)) {
    wdir <- getwd()
  } 
  
  readF1 <- file.path(in.seqDir, paste(in.filename, ".fastq", sep = ""))
  readF1 <- FastqSampler(readF1, n=num)
  set.seed(47)
  readF1.subsample <- yield(readF1)
  close(readF1)
  
  out.file <- 
    file.path(wdir, 'data/processed', 
              paste(in.filename, "_subsample_n", num, ".fastq", sep = ""))
  checkFile(out.file)
  cat("writing \"", out.file, "\" ...\n", sep = "")
  writeFastq(object = readF1.subsample, file = out.file, compress = FALSE, mode = 'w')
  
  if (!is.null(reverse.filename)) {
    readF2 <- file.path(in.seqDir, paste(reverse.filename, ".fastq", sep = ""))
    readF2 <- FastqSampler(readF2, n=num)
    set.seed(47)
    readF2.subsample <- yield(readF2)
    close(readF2)
    out.file <- 
      file.path(wdir, 'data/processed', 
                paste(reverse.filename, "_subsample_n", num, ".fastq", sep = ""))
    checkFile(out.file)
    cat("writing \"", out.file, "\" ...\n", sep = "")
    writeFastq(object = readF2.subsample, file = out.file, compress = FALSE, mode = 'w')
  } 
  
}

checkFile <- function(in.file) {
  
  if (file.exists(in.file)) {
    file.remove(in.file)
    warning("Previous existing file \"", in.file, 
            "\" is about to be overwritten. \n")
  }
  
}


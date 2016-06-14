## Low-quality bases and uncalled bases (i.e. N's) are removed from either both 
#  ends (qtrim.3end = 0) or from the right end (qtrim.3end = 1).
#  Intermediate files are written to folder "data/"
Run1_qualTrimming <- function(readF1, in.filename, seq.mode, qtrim.3end=1,
                              qtrim.thold=10, out.filename, out.ssplot=FALSE,
                              readF2=NULL, reverse.filename=NULL, wdir=NULL, 
                              modseq.dir=NULL, out.dir=NULL) {   
  
  ## Function arguments
  # readF2            Expected when seq.mode = "PE"
  # reverse.filename  Expected when seq.mode = "PE"
  # wdir              (optional) Directory to search for input files, by 
  #                   default current working directory.
  # out.dir           (optional) Path to output files, by default current 
  #                   working directory.
  # modseq.dir        (optional) modseq's head directory, by default current 
  #                   working directory.
  
  ## Whenever input or output directories are not specified, assumed to be the 
  #  current working directory
  if (is.null(wdir)) {
    wdir <- getwd()
  } 
  if (is.null(modseq.dir)) {
    modseq.dir <- getwd()
    warning("Object \'modseq.dir\' not found, set to: \"", modseq.dir, "\".")
  }
  if (is.null(out.dir)) {
    out.dir <- getwd()
  }
  
  if (!existsFunction("IlluminaStat")) {
    source(file.path(modseq.dir, "R/functions/IlluminaStat.R"))
  }
  if (!existsFunction("PlotQualityDistribution")) {
    source(file.path(modseq.dir, "R/functions/PlotQualityDistribution.R"))
  }
  if (!existsFunction("PlotReadLengthDistribution")) {
    source(file.path(modseq.dir, "R/functions/PlotReadLengthDistribution.R"))
  }
  
  if (qtrim.3end == 0) {
    leftTrim <- TRUE
    rightTrim <- TRUE
    cat("Trimming low-quality bases from both ends...\n")
  } else if (qtrim.3end == 1) {
    leftTrim <- FALSE
    rightTrim <- TRUE
    cat("Trimming low-quality bases from the right end...\n")
  }
  
  ## Quality trimming: forward reads (PE mode) or sequencing reads (SE mode)
  readF1.ntrim <- trimEnds(sread(readF1), "N", left = leftTrim, 
                           right = rightTrim, relation = "==", ranges = TRUE)
  readF1.ntrim <- narrow(readF1, start(readF1.ntrim), end(readF1.ntrim))
  readF1.nQtrim <- trimEnds(readF1.ntrim, a = rawToChar(as.raw(qtrim.thold + 33)),
                            left = leftTrim, right = rightTrim)

  ## Creating directories
  cat("Trimmed reads: writing fastq/fasta files...\n")
  if (!file.exists(file.path(wdir, 'data/'))) {
    dir.create(file.path(wdir, 'data/'))
  }
  if (!file.exists(file.path(wdir, 'data/processed/'))) {
    dir.create(file.path(wdir, 'data/processed/'))
  }
  if (!file.exists(file.path(wdir, 'data/processed/fasta/'))) {
    dir.create(file.path(wdir, 'data/processed/fasta/'))
  }
  
  ## Writing fasta file
  out.file <- file.path(wdir, "data/processed/fasta", 
                       paste(out.filename, "_1_nQtrim.fasta", sep = ""))
  checkFile(out.file)
  cat("writing \"", out.file, "\" ...\n", sep = "")
  writeFasta(object = readF1.nQtrim, file = out.file, mode = 'w')
  
  ## Writing fastq file: reads with zero length after trimming are omitted
  out.file <- 
    file.path(wdir, 'data/processed', 
              paste(out.filename, "_subSet_wo0_1_nQtrim.fastq", sep = ""))
  checkFile(out.file)
  cat("writing \"", out.file, "\" ...\n", sep = "")
  writeFastq(object = readF1.nQtrim[width(readF1.nQtrim) != 0], 
             file = out.file, compress = FALSE, mode = 'w')
  
  out.file <- paste(out.filename, "_1_nQtrim", sep = "")
  PlotQualityDistribution(readF1.nQtrim, out.file, out.dir, 
                          title = "Trimmed reads")
  PlotReadLengthDistribution(readF1.nQtrim, out.file, out.dir, 
                             title = "Trimmed reads")
  
  if (seq.mode == "PE") {
    
    ## Quality trimming: reverse reads (PE mode)
    readF2.ntrim <- trimEnds(sread(readF2), "N", left = leftTrim, 
                             right = rightTrim, relation = "==", ranges = TRUE)
    readF2.ntrim <- narrow(readF2, start(readF2.ntrim), end(readF2.ntrim))
    readF2.nQtrim <- trimEnds(readF2.ntrim, a = rawToChar(as.raw(qtrim.thold + 33)),
                              left = leftTrim, right = rightTrim)
    
    ## Writing fastq files: read with zero length after trimming, in any of the 
    # files (forward or reverse), are appended at the end. 
    out.file <- file.path(wdir, 'data/processed', 
                         paste(out.filename, "_1_nQtrim.fastq", sep = ""))
    checkFile(out.file)
    cat("writing \"", out.file, "\" ...\n", sep = "")
    writeFastq(
      object = append(
        readF1.nQtrim[width(readF1.nQtrim) != 0 & width(readF2.nQtrim) != 0],
        readF1.nQtrim[width(readF1.nQtrim) == 0 | width(readF2.nQtrim) == 0]),
      file = out.file, compress = FALSE, mode = 'w'
      )
    
    out.file <- file.path(wdir, 'data/processed', 
                         paste(out.filename, "_2_nQtrim.fastq", sep = ""))
    checkFile(out.file)
    cat("writing \"", out.file, "\" ...\n", sep = "")
    writeFastq(
      object = append(
        readF2.nQtrim[width(readF1.nQtrim) != 0 & width(readF2.nQtrim) != 0], 
        readF2.nQtrim[width(readF1.nQtrim) == 0 | width(readF2.nQtrim) == 0]),
      file = out.file, compress = FALSE, mode = 'w'
      )
    
    ## Writing fasta file
    out.file <- file.path(wdir, 'data/processed/fasta', 
                          paste(out.filename, "_2_nQtrim.fasta", sep = ""))
    checkFile(out.file)
    cat("writing \"", out.file, "\" ...\n", sep = "")
    writeFasta(object = readF2.nQtrim, file = out.file, mode = 'w')
    
    ## Writing fastq file: reads with zero length after trimming are omitted
    out.file <- 
      file.path(wdir, 'data/processed',
                paste(out.filename, "_subSet_wo0_2_nQtrim.fastq", sep = ""))
    checkFile(out.file)
    cat("writing \"", out.file, "\" ...\n", sep = "")
    writeFastq(object = readF2.nQtrim[width(readF2.nQtrim) != 0], 
               file = out.file, compress = FALSE, mode = 'w')
    
    out.file <- paste(out.filename, "_2_nQtrim", sep = "")
    PlotQualityDistribution(readF2.nQtrim, out.file, out.dir, 
                            title = "Trimmed reads")
    PlotReadLengthDistribution(readF2.nQtrim, out.file, out.dir, 
                               title = "Trimmed reads")
    
    if (out.ssplot) {
      
      forward.filename <- in.filename
      cat("Trimmed reads: summary statistics.\n") 
      # Adapted from script written by S.Schmitt
      plot.input <- 
        data.frame(width(readF1), width(readF2), width(readF1.nQtrim), 
                   width(readF2.nQtrim))
      names(plot.input) <- c(forward.filename, reverse.filename,
                             paste(out.filename, "_1_nQtrim", sep = ""),
                             paste(out.filename, "_2_nQtrim", sep = ""))
      mean.readF1 <- mean(plot.input[[1]], na.rm = TRUE)
      mean.readF2 <- mean(plot.input[[2]], na.rm = TRUE)
      mean.readF1.nQtrim <- mean(plot.input[[3]], na.rm = TRUE)
      mean.readF2.nQtrim <- mean(plot.input[[4]], na.rm = TRUE)
      avQScore.readF1.nQtrim <- cbind(
        width(readF1.nQtrim), alphabetScore(readF1.nQtrim) / width(readF1.nQtrim))
      avQScore.readF2.nQtrim <- cbind(
        width(readF2.nQtrim), alphabetScore(readF2.nQtrim) / width(readF2.nQtrim))
      
      cat("Trimmed reads: summary statistics - exporting plots ...\n")
      Ti <- Sys.time()
      IlluminaStat(in.dir = file.path(wdir, 'data/processed/'), 
                   pattern = paste(out.filename, "_subSet_wo0", sep = ""), 
                   SV_plotinputframe = plot.input, SV_qcwd1 = avQScore.readF1.nQtrim,
                   SV_qcwd2 = avQScore.readF2.nQtrim, SV_meanread1 = mean.readF1,
                   SV_meanread2 = mean.readF2, SV_meanread3 = mean.readF1.nQtrim,
                   SV_meanread4 = mean.readF2.nQtrim, 
                   pdfname = paste("stat_trimmed_", out.filename, "_", 
                                   format(Ti, "%Y%m%d_%H%M"), ".pdf", sep = ""),
                   type = "unpaired", sample = FALSE, in.dirRaw = in.seqDir,
                   forward.filename = forward.filename, 
                   reverse.filename = reverse.filename, out.dir = out.dir) 
    }  
    
    cat("Number of reads after quality trimming: ", length(readF1.nQtrim), 
        " x 2 (mode: PE).\n", sep = "")
    cat("Number of pairs of reads after quality trimming with non-zero length: ",
        sum(width(readF1.nQtrim) != 0 & width(readF2.nQtrim) != 0), ".\n",
        sep = "")
    cat("Forward-reads: number of reads after quality trimming with non-zero length: ",
        sum(width(readF1.nQtrim) != 0), ".\n", sep = "")
    cat("Reverse-reads: number of reads after quality trimming with non-zero length: ", 
        sum(width(readF2.nQtrim) != 0), ".\n", sep = "")
    
    return(list("readF1" = readF1.nQtrim, "readF2" = readF2.nQtrim))
    
  } else if (seq.mode == "SE") {
    
    ## Writing fastq files: reads with length zero after trimming are appended
    # at the end 
    out.file <- file.path(wdir, 'data/processed', 
                          paste(out.filename, "_1_nQtrim.fastq", sep = ""))
    checkFile(out.file)
    cat("writing \"", out.file, "\" ...\n", sep = "")
    writeFastq(
      object = append(readF1.nQtrim[width(readF1.nQtrim) != 0], 
                      readF1.nQtrim[width(readF1.nQtrim) == 0]), 
      file = out.file, compress = FALSE, mode = 'w'
      )  
    
    if (out.ssplot) {
      cat("Trimmed reads: summary statistics.\n") 
      # Adapted from script written by S.Schmitt
      plot.input <- data.frame(width(readF1), width(readF1.nQtrim))
      names(plot.input) <- c(in.filename, 
                             paste(out.filename, "_1_nQtrim" , sep = ""))
      avQScore.readF1 <- cbind(width(readF1), alphabetScore(readF1) / width(readF1))
      avQScore.readF1.nQtrim <- cbind(
        width(readF1.nQtrim), alphabetScore(readF1.nQtrim) / width(readF1.nQtrim))
      mean.readF1 <- mean(plot.input[[1]], na.rm = TRUE)
      mean.readF1.nQtrim <- mean(plot.input[[2]], na.rm = TRUE)
      
      cat("Trimmed reads: summary statistics - exporting plots ...\n")
      Ti <- Sys.time()
      IlluminaStat(in.dir = file.path(wdir, 'data/processed/'),
                   pattern = paste(out.filename, "_subSet_wo0", sep = ""), 
                   SV_plotinputframe = plot.input, SV_qcwd1 = avQScore.readF1,
                   SV_qcwd2 = avQScore.readF1.nQtrim, SV_meanread1 = mean.readF1,
                   SV_meanread2 = mean.readF1.nQtrim, SV_meanread3 = 0, 
                   SV_meanread4 = 0,
                   pdfname = paste("stat_trimmed_", out.filename, "_",
                                   format(Ti, "%Y%m%d_%H%M"), ".pdf", sep = ""),
                   type = "single", sample = FALSE, in.dirRaw = in.seqDir,
                   forward.filename = in.filename, out.dir = out.dir)
    }
    
    cat("Number of reads after quality trimming: ", length(readF1.nQtrim), ".\n",
        sep = "")
    cat("Number of reads after quality trimming with non-zero length: ",
        sum(width(readF1.nQtrim) != 0), ".\n", sep = "")
    
    return(readF1.nQtrim)
    
  }
}

checkFile <- function(in.file) {
  
  if (file.exists(in.file)) {
    file.remove(in.file)
    warning("Previous existing file \"", in.file, 
            "\" is about to be overwritten. \n")
  }
  
}

##**TODO** Illumina-stat plots for the case of input directory provided by user 
#          Use threads (option -T)
Run2_peAssembly <- function(readF1, readF2, wdir, out.filename.run2, qtrim.flag,
                            out.ssplot, forward.file, reverse.file, 
                            in.seqDir=NULL, len.rawData=NULL) { 
  ## Creating directories
  if (!file.exists(file.path(wdir, 'data/'))) {
    dir.create(file.path(wdir, 'data/'))
  }
  if (!file.exists(file.path(wdir, 'data/paired/'))) {
    dir.create(file.path(wdir, 'data/paired/'))
  }
  if (!file.exists(file.path(wdir, 'data/paired/fasta/'))) {
    dir.create(file.path(wdir, 'data/paired/fasta/'))
  }
  
  ## Pandaseq options
  # -C filter  Load a pluggable filter module.
  # -B         Allow unbarcoded sequences.
  # -F         Output FASTQ.
  # -w file    Output seqences to a FASTA/FASTQ.
  # -g file    Output log to a text file.
  # -U file    File to write unalignable read pairs with quality scores.
  if (qtrim.flag == 0) {   
    ## Paired-end reads assembly: untrimmed reads
    cat("Paired-end read assembly of untrimmed reads. \n")
    f1 <- file.path(in.seqDir, paste(forward.filename, '.fastq', sep = ""))
    f2 <- file.path(in.seqDir, paste(reverse.filename, '.fastq', sep = ""))   
  } else if (qtrim.flag == 1) {   
    ## Paired-end reads assembly: trimmed reads
    cat("Paired-end read assembly of trimmed reads. \n")
    f1 <- forward.file
    f2 <- reverse.file   
  }
  cat("Forward-reads file: \"", f1, "\".\n", sep = "")
  cat("Reverse-reads file: \"", f2, "\".\n", sep = "")
  assign(paste("runPANDAseq_", out.filename.run2, sep = ""), 
         paste("pandaseq -f ", f1, " -r ", f2, " -C qualString:data/paired", 
               out.filename.run2, "_PANDAseq.fastq -B -F -w data/paired/fasta/",
               out.filename.run2, "_PANDAseq.fasta -g data/paired/logPANDAseq_",
               out.filename.run2, ".txt -U data/paired/", out.filename.run2,
               "_unalignedPANDAseq.fastq", sep = ""))
  ## Use when the plugin (qualString) is not available 
  #  (see below how to output fasta)
  #assign(paste("runPANDAseq_", out.filename.run2, sep = ""),
  #       paste("pandaseq -f ", f1, " -r ", f2, " -B -F -w data/paired/", 
  #             out.filename.run2, "_PANDAseq.fastq -g data/paired/logPANDAseq_",
  #             out.filename.run2, ".txt -U data/paired/", out.filename.run2,
  #             "_unalignedPANDAseq.fastq", sep = ""))
  cat("Running PANDAseq: PAired-eND Assembler for illumina sequences ...\n")
  system(get(paste("runPANDAseq_", out.filename.run2, sep = "")))
  cat("Paired-end read assembly done!\n")
  cat("Output files: \n")
  cat("Fastq-file: \"", wdir, "data/paired/", out.filename.run2, 
      "_PANDAseq.fastq\",\n", sep = "")
  cat("Fasta-file: \"", wdir, "data/paired/fasta/", out.filename.run2,
      "_PANDAseq.fasta\",\n", sep = "")
  cat("Unaligned reads: \"", wdir, "data/paired/", out.filename.run2, 
      "_unalignedPANDAseq.fastq\",\n", sep = "")
  cat("Log-file: \"", wdir, "data/paired/logPANDAseq_", out.filename.run2, 
      ".txt\".\n", sep = "")
  
  PandaseqPaired <- 
    readFastq(dirPath = file.path(wdir, 'data/paired'), 
              pattern = paste(out.filename.run2, "_PANDAseq.fastq", sep = ""))
  ## Temporary work-around (when plugin qualString is not available)
  #writeFasta(PandaseqPaired, 
  #           file = file.path(wdir, 'data/paired/fasta', paste(out.filename.run2,
  #                            "_PANDAseq.fasta", sep = ""))
  #           )
  
  if (out.ssplot) {
    cat("Paired-end assembled reads: summary statistics. \n")
    plot.input1 <- width(readF1) 
    plot.input2 <- width(readF2) 
    plot.input3 <- width(PandaseqPaired)
    avQScore.readF1 <- cbind(plot.input1, alphabetScore(readF1) / plot.input1)
    avQScore.readPE <- cbind(plot.input3, 
                             alphabetScore(PandaseqPaired) / plot.input3)
    n <- max(length(plot.input1), length(plot.input2), length(plot.input3))  
    length(plot.input1) <- n
    length(plot.input2) <- n
    length(plot.input3) <- n
    plot.input <- data.frame(plot.input1, plot.input2, plot.input3)
    mean.readF1 <- mean(plot.input[[1]], na.rm = TRUE)
    mean.readF2 <- mean(plot.input[[2]], na.rm = TRUE)
    mean.readPE <- mean(plot.input[[3]], na.rm = TRUE)
    cat("Paired reads: summary statistics - exporting plots...\n")
    if (qtrim.flag == 0) {
      names(plot.input) <- c(forward.file, reverse.file, 
                             paste(out.filename.run2, "_PANDAseq", sep = ""))
      IlluminaStat(in.dir = file.path(wdir, 'data/paired/'), 
                   pattern = out.filename.run2, SV_plotinputframe = plot.input,
                   SV_qcwd1 = avQScore.readF1, SV_qcwd2 = avQScore.readPE, 
                   SV_meanread1 = mean.readF1, SV_meanread2 = mean.readF2, 
                   SV_meanread3 = mean.readPE, SV_meanread4 = 0,
                   pdfname = paste("stat_assembled_", out.filename.run2, "_", 
                                   format(Ti, "%Y%m%d_%H%M"), ".pdf", sep = ""), 
                   type = "paired", out.dir = out.dir, sample = FALSE, 
                   in.dirRaw = in.seqDir, forward.filename = forward.file)
    } else {
      names(plot.input) <- c(paste(out.filename.run1, "_1", sep = ""),
                             paste(out.filename.run1, "_2", sep = ""),
                             paste(out.filename.run2, "_PANDAseq", sep = ""))
      if (file.exists(
        file.path(wdir, 'data/processed', 
                  paste(out.filename.run1, "_subSet_wo0_1_nQtrim.fastq",
                        sep = "")))) {
        
        IlluminaStat(in.dir = file.path(wdir, 'data/paired/'), 
                     pattern = out.filename.run2, SV_plotinputframe = plot.input,
                     SV_qcwd1 = avQScore.readF1, SV_qcwd2 = avQScore.readPE, 
                     SV_meanread1 = mean.readF1, SV_meanread2 = mean.readF2, 
                     SV_meanread3 = mean.readPE, SV_meanread4 = 0,
                     pdfname = paste("stat_assembled_", out.filename.run2, "_", 
                                     format(Ti, "%Y%m%d_%H%M"), ".pdf", sep = ""), 
                     type = "paired", out.dir = out.dir, sample = FALSE, 
                     in.dirRaw = file.path(wdir, 'data/processed/'), 
                     forward.filename = paste(out.filename.run1, "_subSet_wo0_1",
                                              "_nQtrim", sep = ""))
        
      } else if (file.exists(
        file.path(in.seqDir, paste(out.filename.run1, "_subSet_wo0_1", 
                                   "_nQtrim", sep = "")))) {
                
        IlluminaStat(in.dir = file.path(wdir, 'data/paired/'), 
                     pattern = out.filename.run2, SV_plotinputframe = plot.input,
                     SV_qcwd1 = avQScore.readF1, SV_qcwd2 = avQScore.readPE, 
                     SV_meanread1 = mean.readF1, SV_meanread2 = mean.readF2, 
                     SV_meanread3 = mean.readPE, SV_meanread4 = 0,
                     pdfname = paste("stat_assembled_", out.filename.run2, "_", 
                                     format(Ti, "%Y%m%d_%H%M"), ".pdf", sep = ""), 
                     type = "paired", out.dir = out.dir, sample = FALSE, 
                     in.dirRaw = in.seqDir, 
                     forward.filename = paste(out.filename.run1, "_subSet_wo0_1",
                                              "_nQtrim", sep = ""))
        
      } else {
        #**TODO**
        # Error when forward.file have reads with length zero 
        # Sln: Search for "_subSet_wo0_1" in the in.seqDir
        #system(paste("ls ", file.path(in.seqDir), "*_subSet_wo0_1_nQtrim", sep = ""))
        
      }
      
    }
  }
  len.peData <- length(PandaseqPaired)
  if (is.null(len.rawData)) {
    cat("Number of Paired-end-reads: ", len.peData, ".\n", sep = "")
  } else {
    cat("Number of Paired-end-reads: ", len.peData, " of ", len.rawData, " (",
        round(len.peData * 100 / len.rawData, digit = 2), "%).\n", sep = "")
  }
}
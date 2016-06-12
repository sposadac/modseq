GetReadsFile <- function(wdir) {
  
  ## Works on global variables as some of them have been removed depending on 
  #  the sequencing mode. 
  if (paired.flag == 1) {  
    
    ## Paired reads 
    cat("Read mapping of paired-end reads.\n")
    if (map.mode == "gls") {
      # Must be format = "fasta", because when a FASTQ file is given as input 
      # to readDNAStringSet, all the sequences must have the same length 
      in.file <- 
        file.path(wdir, "data/paired/fasta", 
                  paste(out.filename.run2, "_PANDAseq.fasta", sep = ""))
      format <- "fasta"
      
    } else if (map.mode == "bwa") {
      in.file <- 
        file.path(wdir, "data/paired", 
                  paste(out.filename.run2, "_PANDAseq.fastq", sep = ""))
      format <- "fastq"
      
    } else {
      stop("Unsupported mapping mode, \"", map.mode, "\". \n")
    }

    
  } else if (paired.flag == 0) {  
    
    ## Unpaired reads
    if (qtrim.flag == 0) {  
      
      ## Unpaired and untrimmed reads - raw reads 
      if (seq.mode == "SE") {
        cat("Read mapping of untrimmed and non-paired reads. \n")
        in.file <- file.path(in.seqDir, paste(in.filename, ".fastq", sep = ""))
        
      } else if (seq.mode == "PE") {      
        cat("Read mapping of untrimmed and unpaired reads.\n")
        
        if (paired.file == "f") {
          cat("Set of forward reads selected.\n")
          in.file <- file.path(in.seqDir,  
                               paste(forward.filename, ".fastq", sep = ""))
          
        } else if (paired.file == "r") {
          cat("Set of reverse reads selected.\n")
          in.file <- file.path(in.seqDir,  
                               paste(reverse.filename, ".fastq", sep = ""))
          
        }
        format <- "fastq"
      }
      
    } else {   
      
      ## Unpaired, but trimmed reads
      cat("Read mapping of trimmed but non-paired reads.\n")
      if (seq.mode == "SE" || paired.file == "f") {
        
        if (map.mode == "gls"){
          in.file <- 
            file.path(wdir, "data/processed/fasta", 
                      paste(out.filename.run1, "_1_nQtrim.fasta", sep = ""))
          format <- "fasta"
          
        } else if (map.mode == "bwa") {
          in.file <- file.path(wdir, "data/processed", 
                               paste(out.filename.run1, "_1_nQtrim.fastq", sep = "")) 
          format <- "fastq"
          
        } else {
          stop("Unsupported mapping mode, \"", map.mode, "\". \n")
        }
        
        
      } else if (paired.file == "r") {      
        
        if (map.mode == "gls"){
          in.file <- 
            file.path(wdir, 'data/processed/fasta',
                      paste(out.filename.run1, "_2_nQtrim.fasta", sep = ""))
          format <- "fasta"
          
        } else if (map.mode == "bwa") {
          in.file <- file.path(wdir, "data/processed", 
                               paste(out.filename.run1, "_2_nQtrim.fastq", sep = "")) 
          format <- "fastq"
          
        } else {
          stop("Unsupported mapping mode, \"", map.mode, "\". \n")
        }
        
      }
    }
  }
  
  return(list(in.file, format))
  
}
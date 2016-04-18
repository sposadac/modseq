### MAIN SCRIPT ###
# Author: Susana Posada Cespedes <susana.posada@bsse.ethz.ch>
# April 2016
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
if (file.exists(paste(out.dir, "/out_", out.filename.run2, "_run",
                      paste(names(run)[which(run == 1)], collapse = "_"),
                      ".txt",sep=""))) {
  output <- file(paste(out.dir, "/out_", out.filename.run2, "_run",
                       paste(names(run)[which(run == 1)], collapse = "_"),
                       ".txt", sep=""), open = "at")
} else {
  output <- file(paste(out.dir, "/out_", out.filename.run2, "_run",
                       paste(names(run)[which(run == 1)], collapse = "_"),
                       ".txt",sep=""), open = "wt")
}
sink(output, append = TRUE)

cat("Run", paste(names(run)[which(run == 1)], collapse = ","), "started...\n")
keep <- ls()

# Read fastq-files
if (run[1] == 1 || (run[2] == 1 && qtrim.flag == 0)) { 
  cat("Sequencing mode: \"", seq.mode, "\".\n", sep = "")
  if (seq.mode == "PE") {
    readF1 <- readFastq(dirPath = in.seqDir, pattern = forward.filename)
    readF2 <- readFastq(dirPath = in.seqDir, pattern = reverse.filename)
    len.rawData <- length(readF1)
    cat("Number of reads: ", len.rawData, " x 2 (mode: PE).\n", sep = "")
  } else if (seq.mode == "SE") {
    readF1 <- readFastq(dirPath = in.seqDir, pattern = in.filename)
    len.rawData <- length(readF1)
    cat("Number of reads: ", len.rawData, ".\n", sep = "")
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
  
  cat("Module-table used: \"", mod.filename, ".csv\". \n", sep = "")
  patterns <- 
    data.frame(read.csv(paste(in.modDir, mod.filename, ".csv", sep = ""),
                        header = TRUE, sep = ",", dec = "."), 
               stringsAsFactors = FALSE)
  num.cores <- detectCores()
  if (paired.flag == 1) {  # Paired reads (i.e. "PE")
    if (file.exists(paste("./DATA/PAIRED/FASTA/", out.filename.run2,
                          "_PANDAcombined.fasta", sep = ""))) {
      reads.subj <- readDNAStringSet(
        file = paste("./DATA/PAIRED/FASTA/", out.filename.run2, 
                     "_PANDAcombined.fasta", sep = ""), format = "fasta", 
        use.names = TRUE)
    } else {
      stop("File \"", paste("./DATA/PAIRED/FASTA/", out.filename.run2, 
                            "_PANDAcombined.fasta", sep = ""), "\" not found. \n")
    }
  } else if (paired.flag == 0) {  # Unpaired reads
    if (qtrim.flag == 0) {  # Unpaired and untrimmed reads - raw reads
      if (seq.mode == "SE") {
        reads.subj <- readDNAStringSet(
          file = paste(in.seqDir, in.filename, sep = ""), format = 
            "fastq", use.names = TRUE) 
      } else if (seq.mode == "PE") {
        if (paired.file == "f") {
          reads.subj <- readDNAStringSet(
            file = paste(in.seqDir, forward.filename, sep = ""), 
            format = "fastq", use.names = TRUE) 
        } else if (paired.file == "r") {
          reads.subj <- readDNAStringSet(
            file = paste(in.seqDir, reverse.filename, sep = ""), 
            format = "fastq", use.names = TRUE) 
        }
      }        
    } else {  # Unpaired, but trimmed reads 
      if (seq.mode == "SE" || paired.file == "f") {
        if (file.exists(paste("./DATA/PROCESSED/FASTA/", out.filename.run1, 
                              "_1_nQtrim.fasta", sep = ""))) {
          reads.subj <- readDNAStringSet(
            file = paste("./DATA/PROCESSED/FASTA/", out.filename.run1, 
                         "_1_nQtrim.fasta", sep = ""), format = "fasta", 
            use.names = TRUE) 
        } else {
          stop("File \"", paste("./DATA/PROCESSED/FASTA/", out.filename.run1,
                                "_1_nQtrim.fasta", sep = ""), "\" not found", 
               ".\n")
        }
      } else if (paired.file == "r") {
        if (file.exists(paste("./DATA/PROCESSED/FASTA/", out.filename.run1, 
                              "_2_nQtrim.fasta", sep = ""))) {
          reads.subj <- readDNAStringSet(
            file = paste("./DATA/PROCESSED/FASTA/", out.filename.run1, 
                         "_2_nQtrim.fasta", sep = ""), format = "fasta", 
            use.names = TRUE) 
        } else {
          stop("File \"", paste("./DATA/PROCESSED/FASTA/", out.filename.run1,
                                "_2_nQtrim.fasta", sep = ""), "\" not found",
               ".\n")
        }
      }
    }
  }
  reads.subj.len <- length(reads.subj)
  
  if (map.mode == "gls") {
    
    # TODO: verbose level (display info per search round?), 
    # TODO: map.mode pairwise alignment
    
    source(paste(wdir, '/functions/Search_vmatchV2.R', sep = ""))
    Nsearch <- "Ambiguous matches enabled"
    mod.tot <- ncol(patterns)
    
    if (gls.ambiguity == FALSE) {
      Nsearch <- "Ambiguous matches disabled"
      cat("Ambiguous matches disabled for the pattern search (i.e. N only to ",
          "N). \n")
    }
    
    ti.search <- Sys.time()
    res.list <- list()
    res.list[[1]] <-  names(reads.subj)
    names(res.list)[1] <- "" 
    res.counts <- numeric(mod.tot + 1)
    res.counts[1] <- reads.subj.len
    if (mem.trace) {
      memTrace <- c(memTrace, MemTrace())
      offset_mem <- length(memTrace)
      memList <- format(object.size(res.list), units = "MB")
    }
    cat("GLS - Starting pattern search ...\n")
    if (gls.direction == "f") {
      cat("Search mode: Forward - gls. \n")
      for (i in seq_len(mod.tot)) {
        res.list <- mclapply(
          names(res.list), FUN = SearchPatListID, hitsList = res.list, 
          pat = patterns[i], input.data = reads.subj, mm = gls.mma,  
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
          pat = patterns[i], input.data = reads.subj, mm = gls.mma, mode = "r", 
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
    save(res.list, 
         file = paste(out.dir, "resList_", res.listName, ".rda", sep = ""))
    cat("Pattern search completed!\n")
    if (mem.trace) {
      memTrace <- c(memTrace, MemTrace())
    }
    
    cat("Saving and plotting search-results ...\n")
    # Reformat counts as data.frame
    res.counts <- 
      data.frame("Round" = c("Input-data", paste("Mod. ", seq_len(mod.tot), "/",
                                                 mod.tot, sep ="")),
                 "Counts" = res.counts)  
    if (gls.direction == "r") { 
      res.counts$Round <- 
        factor(res.counts$Round, 
               levels = c(levels(res.counts$Round)[1], 
                          rev(levels(res.counts$Round)[2:(mod.tot+1)])))
      print(paste("Found ", res.counts[2, 2], " module combinations in ",
                  reads.subj.len, " reads (", 
                  round(res.counts[2, 2] * 100 / reads.subj.len, 
                        digits = 2), "%)", sep = ""))
    } else {
      print(paste("Found ", res.counts[nrow(res.counts), 2], 
                  " module combinations in ", reads.subj.len, " reads (",
                  round(res.counts[nrow(res.counts), 2] * 100 / reads.subj.len, 
                        digits = 2), "%)", sep = ""))
    }
    cat("Printing module-counts plot as \"", out.dir, res.counts.filename,
        "_graphs.pdf\".\n", sep = "")
    pl.patCounts <- ggplot(res.counts, aes(x = Round, y = Counts, 
                                           ymax = max(Counts)*1.1))
    pl.patCounts <- pl.patCounts + geom_bar(stat = "identity") + 
      labs(x = "Search round", y = "Counts") +  
      ggtitle(paste("ModSeq | Module-counts (", Nsearch,")", 
                    sep = "")) + theme_bw() + 
      theme(text = element_text(size = 14), 
            plot.title = element_text(face = "bold", size = 16)) + 
      annotate("text", x = Inf, y = -Inf, label = run.info, hjust = 1.1, 
               vjust = -0.3, cex = 3) + 
      geom_text(aes(label = Counts), size = 5, hjust = 0.5, vjust = -0.2, 
                position = "stack")
    ggsave(filename = paste(oout.dir, res.counts.filename, "_graph.pdf", 
                            sep = ""), plot = pl.patCounts, paper = "a4r", 
           width = 12)
    
    ### Distribution of modular variants per search round
    mod.number <- sapply(patterns, function(x) length(na.omit(x)))
    mod.number.cumsum <- cumsum(c(1,mod.number))
    mod.empty <- 
      unlist(mapply(function(x,i) if (sum(na.omit(x) == "") > 0) return(i),
                    x = patterns, i = seq_len(mod.tot)))
    if (gls.direction == "f") {
      mod.comb.cumprod <- cumprod(mod.number)
    } else if (gls.direction == "r") {
      mod.comb.cumprod <- rev(cumprod(rev(mod.number)))
    } 
    
    var.counts <- 
      data.frame(matrix(ncol = 3, 
                        dimnames = list(c(""), c("module","variant","freq"))))
    var.counts[1, ] <- c(0, NA, reads.subj.len)
    var.counts[2:(sum(mod.number) + 1), 1] <- rep(seq(mod.tot), mod.number)
    var.counts[2:(sum(mod.number) + 1), 2] <- 
      unlist(mclapply(seq_len(mod.tot), 
                      function(x) rownames(na.omit(patterns[x])), 
                      mc.cores = num.cores))
    for (i in seq_len(mod.tot)) {
      mod.aux <- 
        data.frame(read.csv(paste(out.dir, res.counts.filename, '_round_',
                                  i, '_of_', mod.tot, '.csv', sep = ""), 
                            nrow = mod.comb.cumprod[i], header = TRUE, 
                            row.names = 1))
      var.counts[(mod.number.cumsum[i] + 1):(mod.number.cumsum[i] + mod.number[i]), 3] <-
        unlist(mclapply(seq_len(mod.number[i]), 
                        function(j) sum(mod.aux[seq(j, mod.comb.cumprod[i], 
                                                    mod.number[i]), 1]),
                        mc.cores = num.cores))
    }
    rm(mod.aux)
    var.counts$module <- as.factor(var.counts$module)
    if (gls.direction == "r") {  
      var.counts$module <- 
        factor(var.counts$module, 
               levels = c(levels(var.counts$module)[1], 
                          rev(levels(var.counts$module)[2:(mod.tot+1)])))
      mod.number <- rev(mod.number)
    } 
    var.counts$variant[
      unlist(mclapply(names(which(table(var.counts$module) == 1)),
                      function(x) which(x == var.counts$module), 
                      mc.cores = num.cores))] <- NA
    cat("Printing modular variant counts plot as \"",out.dir,
        res.counts.filename, "_modularVariantsDistribution_graph.pdf\".\n",
        sep = "")
    pl.patCounts <- ggplot(var.counts, aes(x = module, y = freq, 
                                           ymax = max(res.counts$Counts) * 1.1))
    pl.patCounts <- pl.patCounts + 
      geom_bar(stat = "identity", aes(fill = factor(variant), colour = "gray"), 
               alpha = 0.8) + labs(x = "Search round", y = "Counts") + 
      scale_x_discrete(breaks = levels(var.counts$module), 
                       labels = res.counts$Round) + 
      scale_color_manual(guide = FALSE, values = "gray") + theme_bw() + 
      ggtitle("ModSeq | Distribution of modular variants per search round") + 
      theme(text = element_text(size = 14),
            plot.title = element_text(face = "bold", size = 16), 
            legend.title = element_blank()) +
      geom_text(data = data.frame(x = seq_len(mod.tot + 1), 
                                  num = c(NA, mod.number),
                                  y = res.counts$Counts), 
                aes(label = y, x = x, y = y, 
                    ymax = max(res.counts$Counts) * 1.1), size = 5, hjust = 0.5,
                vjust = -0.2, position = "stack") + 
      geom_text(data = data.frame(x   = seq_len(mod.tot + 1), 
                                  num = c(NA, mod.number),
                                  y   = rep(0, mod.tot + 1)),
                aes(label = num, x = x, y = y, 
                    ymax = max(res.counts$Counts) * 1.1), size = 4, hjust = 0.5, 
                vjust = 1.1, position = "stack")
    ggsave(paste(out.dir, res.counts.filename, "_modularVariants", 
                 "Distribution_graph.pdf", sep = ""), pl.patCounts, 
           paper = "a4r", width = 12) 
    
    var.counts$module <- c("0", rep(names(patterns), mod.number))
    cat('Exporting: \"', ouout.dir, res.counts.filename, '_modularVariants',
        'Distribution_table.csv. \n', sep = "")
    write.csv(var.counts, paste(out.dir, res.counts.filename, 
                                '_modularVariantsDistribution_table.csv', 
                                sep = ""))
    
    if (sum(gls.mma > 0) > 0) {
      cat("Maximum mismatch allowance for module", which(gls.mma > 0), "was", 
          gls.mma[which(gls.mma > 0)], ".\n")
    }
    
    if (mem.trace) {
      source(paste(wdir, '/functions/plotMemProfile.R', sep = ""))
      plotMemProfile((min(as.numeric(
        unlist(sapply(memList, function(x) strsplit(x, split = "Mb"), 
                      USE.NAMES = FALSE)))) %/% 10) * 10, ((max(as.numeric(
                        unlist(sapply(memList, function(x) strsplit(x, split = "Mb"), 
                                      USE.NAMES = FALSE)))) %/% 10) * 10 + 10), 
        (min(memTrace[offset_mem:(mod.tot + offset_mem)]) %/% 10) * 10, 
        ((max(memTrace[offset_mem:(mod.tot + offset_mem)]) %/% 10) * 10 + 10), 
        res.listName, path = out.dir)
    }
    
    aux <- unlist(res.list, use.names = FALSE)
    if (sum(duplicated(aux)) > 0) {
      rep.tab <- table(aux)
      pdf(paste(out.dir, res.listName, "_repeated_hits_graph.pdf", 
                sep = ""))
      plot(as.numeric(rep.tab), type = "h", xlab = "Reads", ylab = "Counts",
           main = "Read mapping")
      dev.off()
      if (!file.exists(paste(in.modDir, mod.filename, "_modComb.fasta", 
                             sep = ""))) {
        source(paste(wdir, '/functions/ModuleCombinationsGen.R', sep = ""))
        ModuleCombinationsGen(mod.filename, pattern = patterns, 
                              num.cores = num.cores)
      }
      mod.comb <- 
        readDNAStringSet(paste(in.modDir, mod.filename, "_modComb.fasta", 
                               sep = ""))
      mod.comb.len <- length(mod.comb)
      res.list.lengths <- 
        unlist(mclapply(res.list, length, mc.cores = num.cores), 
               use.names = FALSE)
      # To get the name of the module combinations. 
      if (sum(names(res.list) == mod.comb) != mod.comb.len) {
        # Only needed if the order is not preserved
        mod.comb.seq2names <- names(mod.comb)
        names(mod.comb.seq2names) <- mod.comb
        res.list.name <- 
          unlist(mcmapply(function (x,n) rep(x,n), 
                          x = mod.comb.seq2names[names(res.list)], 
                          n = res.list.lengths, USE.NAMES = FALSE,
                          mc.cores = num.cores), use.names = FALSE)       
      } else {
        res.list.name <- 
          unlist(mcmapply(function (x,n) rep(x,n), x = names(mod.comb), 
                          n = res.list.lengths, USE.NAMES = FALSE,
                          mc.cores = num.cores), use.names = FALSE)
      }
      rep.hits <- names(rep.tab[rep.tab > 1])
      aux.rep.hits <- 
        unlist(mclapply(rep.hits, function(x) res.list.name[which(aux == x)],
                        mc.cores= num.cores))
      df.rep.hits <- data.frame("readID" = rep(rep.hits, rep.tab[rep.tab > 1]),
                                "refID"  = aux.rep.hits,
                                "refSEQ" = as.character(mod.comb[aux.rep.hits]))
      
      write.csv(df.rep.hits, paste(out.dir, res.listName, 
                                   "_repeated_hits_table.csv", sep = ""),
                row.names = TRUE)
      print(paste(rep.hits, "- read was counted", rep.tab[rep.tab > 1], 
                  "times"))
      save(df.rep.hits, file = paste(out.dir, res.listName, 
                                     "_repeated_hits.rda", sep = ""))
      cat(sum(rep.tab == 1), "reads mapped exactly 1 time. \n")
      cat(sum(rep.tab > 1), " reads (", 
          round(sum(rep.tab > 1) * 100 / reads.subj.len, digit = 2), 
          "%) mapped > 1 times. \n", sep = "")
      cat("Corrected for multiple hits: Found ", length(rep.tab), 
          " module combinations in ", reads.subj.len, " reads (", 
          round(length(rep.tab) / reads.subj.len * 100, digits = 2), "%).\n",
          sep = "")
      if (run[4] == 1) {
        keep <- append(keep, c("mod.comb", "mod.comb.len", "res.list.lengths",
                               "res.list.name"))
      }
    } else {
      cat("All mapped reads are mapped exactly 1 time (i.e. mapped to 1 ",
          "module combination). \n")
    }
    
    if (run[4] == 1) {
      keep <- append(keep, c("mod.tot", "reads.subj", "reads.subj.len", 
                             "res.list"))
    }
    
  } else if (map.mode == "bwa") {
    
    if (paired.flag == 1) {  # Paired reads
      reads.file <- paste("./DATA/PAIRED/", out.filename.run2, 
                          "_PANDAseq.fastq", sep = "")
    } else if (paired.flag == 0) { 
      if (qtrim.flag == 0) { # Unpaired and untrimmed files - raw files
        if (seq.mode == "SE") {
          reads.file <- paste(in.seqDir, in.filename, sep = "")
        } else if (seq.mode == "PE") {
          if (paired.file == "f") {
            reads.file <- paste(in.seqDir, forward.filename, sep = "")
          } else if (paired.file == "r") {
            reads.file <- paste(in.seqDir, reverse.filename, sep = "")
          }
        }
      } else { # Unpaired, but trimmed files 
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
    # indexing reference sequences - module combinations
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
    cat("Found ", reads.map, " module combinations in ", reads.subj.len, 
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

Run3_gls <- function() {
  # **TODO**: verbose level (display info per search round?), 
  # **TODO**: map.mode pairwise alignment
  
  Nsearch <- "Ambiguous matches enabled"
  mod.tot <- ncol(patterns)
  
  if (gls.ambiguity == FALSE) {
    Nsearch <- "Ambiguous matches disabled"
    cat("Ambiguous matches disabled for the pattern search (i.e. N only to ",
        "N). \n")
  }
  
  ## Pattern search
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
  cat("Pattern search completed!\n")
  out.file <- 
    file.path(out.dir, paste("resList_", res.listName, ".rda", sep = ""))
  cat("Exporting list containig results of the read mapping: \"", out.file,
      "\".\n", sep = "")
  save(res.list, file = out.file)
  
  if (mem.trace) {
    memTrace <- c(memTrace, MemTrace())
  }
  
  cat("Saving and plotting results of the read mapping...\n")
  ## Reformatting counts as data.frame
  res.counts <- 
    data.frame("Round" = 
                 c("Input-data", paste("Mod. ", seq_len(mod.tot), "/", mod.tot, 
                                       sep ="")),
               "Counts" = res.counts)  
  if (gls.direction == "r") { 
    res.counts$Round <- 
      factor(res.counts$Round, 
             levels = c(levels(res.counts$Round)[1], 
                        rev(levels(res.counts$Round)[2:(mod.tot+1)])))
    print(
      paste("Found ", res.counts[2, 2], " module-combinations in ",
            reads.subj.len, " reads (", 
            round(res.counts[2, 2] * 100 / reads.subj.len, digits = 2), "%)",
            sep = "")
    )
  } else {
    print(paste("Found ", res.counts[nrow(res.counts), 2], 
                " module-combinations in ", reads.subj.len, " reads (",
                round(res.counts[nrow(res.counts), 2] * 100 / reads.subj.len, 
                      digits = 2), "%)", sep = ""))
  }
  
  ### Module-counts plot
  out.file <- 
    file.path(out.dir, paste(res.counts.filename, "_graph.pdf", sep = ""))
  cat("Printing module-counts plot as \"", out.file, "\".\n", sep = "")
  pl.patCounts <- 
    ggplot(res.counts, aes(x = Round, y = Counts, ymax = max(Counts) * 1.1))
  pl.patCounts <-
    pl.patCounts + geom_bar(stat = "identity") + 
    labs(x = "Search round", y = "Counts") +  
    ggtitle(paste("ModSeq | Module-counts (", Nsearch,")", sep = "")) +
    theme_bw() + 
    theme(text = element_text(size = 14), 
          plot.title = element_text(face = "bold", size = 16)) + 
    annotate("text", x = Inf, y = -Inf, label = run.info, hjust = 1.1, 
             vjust = -0.3, cex = 3) + 
    geom_text(aes(label = Counts), size = 5, hjust = 0.5, vjust = -0.2, 
              position = "stack")
  ggsave(filename = out.file, plot = pl.patCounts, paper = "a4r", width = 12)
  
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
                      dimnames = list(c(""), c("module", "variant", "freq"))))
  var.counts[1, ] <- c(0, NA, reads.subj.len)
  var.counts[2:(sum(mod.number) + 1), 1] <- rep(seq(mod.tot), mod.number)
  var.counts[2:(sum(mod.number) + 1), 2] <- 
    unlist(mclapply(seq_len(mod.tot), 
                    function(x) rownames(na.omit(patterns[x])), 
                    mc.cores = num.cores))
  for (i in seq_len(mod.tot)) {
    mod.aux <- 
      data.frame(
        read.csv(
          file.path(out.dir, paste(res.counts.filename, '_round_', i, '_of_',
                                   mod.tot, '.csv', sep = "")), 
          nrow = mod.comb.cumprod[i], header = TRUE, row.names = 1)
      )
    var.counts[(mod.number.cumsum[i] + 1):(mod.number.cumsum[i] + mod.number[i]), 3] <-
      unlist(mclapply(
        seq_len(mod.number[i]), 
        function(j) sum(mod.aux[seq(j, mod.comb.cumprod[i], mod.number[i]), 1]),
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
  out.file <- 
    file.path(out.dir, 
              paste(res.counts.filename, "_modularVariantsDistribution", 
                    "_graph.pdf", sep = ""))
  cat("Printing modular-variant counts plot as: \"", out.file, "\".\n", 
      sep = "")
  pl.patCounts <- 
    ggplot(var.counts, aes(x = module, y = freq, 
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
    geom_text(data = 
                data.frame(x = seq_len(mod.tot + 1), num = c(NA, mod.number),
                           y = res.counts$Counts), 
              aes(label = y, x = x, y = y, 
                  ymax = max(res.counts$Counts) * 1.1), size = 5, hjust = 0.5,
              vjust = -0.2, position = "stack") + 
    geom_text(data = 
                data.frame(x = seq_len(mod.tot + 1), num = c(NA, mod.number),
                           y = rep(0, mod.tot + 1)),
              aes(label = num, x = x, y = y, 
                  ymax = max(res.counts$Counts) * 1.1), size = 4, hjust = 0.5, 
              vjust = 1.1, position = "stack")
  ggsave(out.file, pl.patCounts, paper = "a4r", width = 12) 
  
  var.counts$module <- c("0", rep(names(patterns), mod.number))
  out.file <- 
    file.path(out.dir, paste(res.counts.filename, "_modularVariants", 
                             "Distribution_table.csv", sep = ""))
  cat("Exporting distribution of modular variants per search round as table:",
      " \"", out.file, ".\" \n", sep = "")
  write.csv(var.counts, out.file)
  
  if (sum(gls.mma > 0) > 0) {
    cat("Maximum mismatch allowance for module", which(gls.mma > 0), "was", 
        gls.mma[which(gls.mma > 0)], ".\n")
  }
  
  ### Plot of memory consumption
  if (mem.trace) {
    source(file.path(wdir, 'R/functions/plotMemProfile.R'))
    plotMemProfile(
      (min(
        as.numeric(unlist(
          sapply(memList, function(x) strsplit(x, split = "Mb"), 
                 USE.NAMES = FALSE)))) %/% 10) * 10, 
      ((max(as.numeric(unlist(
        sapply(memList, function(x) strsplit(x, split = "Mb"), 
               USE.NAMES = FALSE)))) %/% 10) * 10 + 10), 
      (min(memTrace[offset_mem:(mod.tot + offset_mem)]) %/% 10) * 10, 
      ((max(memTrace[offset_mem:(mod.tot + offset_mem)]) %/% 10) * 10 + 10), 
      res.listName, path = out.dir)
  }
  
  ### Plot of repeated hits (reads that were mapped to several 
  #   module-combination)
  aux <- unlist(res.list, use.names = FALSE)
  if (sum(duplicated(aux)) > 0) {
    
    out.file <- 
      file.path(out.dir, 
                paste(res.listName, "_repeated_hits_graph.pdf", sep = ""))
    cat("Printing plot of repeated-hits in: \"", out.file, "\".\n", sep = "")
    rep.tab <- table(aux)
    pdf(out.file)
    plot(as.numeric(rep.tab), type = "h", xlab = "Reads", ylab = "Counts",
         main = "Read mapping")
    dev.off()
    
    mod.file <- 
      file.path(in.modDir, paste(mod.filename, "_modComb.fasta", sep = ""))
    if (!file.exists(mod.file)) {
      source(file.path(wdir, 'R/functions/ModuleCombinationsGen.R'))
      ModuleCombinationsGen(mod.filename, pattern = patterns, 
                            num.cores = num.cores, in.modDir)
    }
    mod.comb <- readDNAStringSet(mod.file)
    mod.comb.len <- length(mod.comb)
    res.list.lengths <- 
      unlist(mclapply(res.list, length, mc.cores = num.cores), 
             use.names = FALSE)
    
    ## To get the name of the module-combinations. 
    if (sum(names(res.list) == mod.comb) != mod.comb.len) {
      ## Only needed if the order is not preserved
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
    out.file <- 
      file.path(out.dir, 
                paste(res.listName, "_repeated_hits_table.csv", sep = ""))
    cat("Exporting repeated-hits as table: \"", out.file, ".\" \n", sep = "")
    write.csv(df.rep.hits, out.file, row.names = TRUE)
    print(
      paste(rep.hits, "- read was counted", rep.tab[rep.tab > 1], "times"))
    
    out.file <- 
      file.path(out.dir, paste(res.listName, "_repeated_hits.rda", sep = ""))
    cat("Exporting list containig repeated-hits: \"", out.file, "\".\n",
        sep = "")
    save(df.rep.hits, file = out.file)
    cat(sum(rep.tab == 1), "reads mapped exactly 1 time. \n")
    cat(sum(rep.tab > 1), " reads (", 
        round(sum(rep.tab > 1) * 100 / reads.subj.len, digit = 2), 
        "%) mapped > 1 times. \n", sep = "")
    cat("Corrected for multiple hits: Found ", length(rep.tab), 
        " module-combinations in ", reads.subj.len, " reads (", 
        round(length(rep.tab) / reads.subj.len * 100, digits = 2), "%).\n",
        sep = "")
    if (run[4] == 1) {
      keep <- append(keep, c("mod.comb", "mod.comb.len", "res.list.lengths",
                             "res.list.name"))
    }
  } else {
    cat("All mapped reads are mapped exactly 1 time (i.e. mapped to 1 ",
        "module-combination). \n")
  }
  
}
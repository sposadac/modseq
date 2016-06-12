# **TODO**: Exclude duplicates for the distribution (res.sam.real: 1024, 1040)?

Run4_bwa <- function(mod.comb, res.sam.realn=NULL, bwa.dupl=TRUE, coverage.left=0, coverage.right=0, mapQ.thold=8, editDist.thold=8, 
                     num.reads=NULL, out.filename) {
  
  source(file.path(modseq.dir, "R/functions/PlotModuleCombinationsDistribution.R"))
  
  if (is.character("res.sam.realn")) {
    source(file.path(modseq.dir, "R/functions/LoadAlignmentFile.R"))
    res.sam.realn <- 
      LoadAlignmentFile(sam.file = res.sam.realn, bwa.dupl = bwa.dupl)
  }

  ## Mapped reads (excluding secondary alignments)
  #  0 and 16      -> mapped read (forward or reverse strand, respectively)
  #  1024 and 1040 -> reads marked as duplicates (forward or reverse strand, 
  #                   respectively)

  ind.sam.realn <- 
    which(res.sam.realn[["flag"]] == 0 | res.sam.realn[["flag"]] == 16 |
            res.sam.realn[["flag"]] == 1024 | res.sam.realn[["flag"]] == 1040)
  
  ### FILTERING 
  # 1. By coverage (first nucleotide consensus 1 and last of consensus 2)
  leftclipped <- as.integer(res.sam.realn[["refPOS"]][ind.sam.realn] - 1) 
  rightclipped <- 
    width(mod.comb[res.sam.realn[["refID"]][ind.sam.realn]]) - 
    cigarWidthAlongReferenceSpace(res.sam.realn[["CIGAR"]][ind.sam.realn]) -
    leftclipped
  ind.coverage <- 
    which(leftclipped <= coverage.left & rightclipped <= coverage.right) 
  
  # 2. By mapping quality
  ind.mapQ <- 
    which(as.numeric(res.sam.realn[["mapQ"]][ind.sam.realn][ind.coverage]) >=
            mapQ.thold)
  
  # 3. By edit distance
  md.aux <- unlist(mclapply(strsplit(res.sam.realn[["mismatchPOS"]], 
                                     split = ":Z:"), function(x) x[[1]]))
  md.aux2 <- unlist(mclapply(strsplit(res.sam.realn[[13]], split = ":Z:"),
                             function(x) x[[1]]))
  
  res.sam.realn[[12]][ind.sam.realn][which(md.aux[ind.sam.realn] == "XA")] <- 
    res.sam.realn[[13]][ind.sam.realn][which(md.aux[ind.sam.realn] == "XA")]
  res.sam.realn[[12]][ind.sam.realn][
    which(md.aux[ind.sam.realn] == "SA" & md.aux2[ind.sam.realn] == "XA")] <- 
    res.sam.realn[[14]][ind.sam.realn][
      which(md.aux[ind.sam.realn] == "SA" & md.aux2[ind.sam.realn] == "XA")]
  res.sam.realn[[12]][ind.sam.realn][
    which(md.aux[ind.sam.realn] == "SA" & md.aux2[ind.sam.realn] == "MD")] <- 
    res.sam.realn[[13]][ind.sam.realn][
      which(md.aux[ind.sam.realn] == "SA" & md.aux2[ind.sam.realn] == "MD")]
  
  if (bwa.dupl) {
    res.sam.realn[[15]][ind.sam.realn][which(md.aux[ind.sam.realn] == "XA")] <- 
      res.sam.realn[[16]][ind.sam.realn][which(md.aux[ind.sam.realn] == "XA")]   
    res.sam.realn[[15]][ind.sam.realn][
      which(md.aux[ind.sam.realn] == "SA" & md.aux2[ind.sam.realn] == "XA")] <-
      res.sam.realn[[17]][ind.sam.realn][
        which(md.aux[ind.sam.realn] == "SA" & md.aux2[ind.sam.realn] == "XA")]
    res.sam.realn[[15]][ind.sam.realn][
      which(md.aux[ind.sam.realn] == "SA" & md.aux2[ind.sam.realn] == "MD")] <- 
      res.sam.realn[[16]][ind.sam.realn][
        which(md.aux[ind.sam.realn] == "SA" & md.aux2[ind.sam.realn] == "MD")]
  } else {
    res.sam.realn[[14]][ind.sam.realn][which(md.aux[ind.sam.realn] == "XA")] <- 
      res.sam.realn[[15]][ind.sam.realn][which(md.aux[ind.sam.realn] == "XA")]   
    res.sam.realn[[14]][ind.sam.realn][
      which(md.aux[ind.sam.realn] == "SA" & md.aux2[ind.sam.realn] == "XA")] <-
      res.sam.realn[[16]][ind.sam.realn][
        which(md.aux[ind.sam.realn] == "SA" & md.aux2[ind.sam.realn] == "XA")]
    res.sam.realn[[14]][ind.sam.realn][
      which(md.aux[ind.sam.realn] == "SA" & md.aux2[ind.sam.realn] == "MD")] <- 
      res.sam.realn[[15]][ind.sam.realn][
        which(md.aux[ind.sam.realn] == "SA" & md.aux2[ind.sam.realn] == "MD")]  
  }
  
  editDistance <- 
    as.integer(unlist(mclapply(strsplit(
      res.sam.realn[["editDistance"]][ind.sam.realn][ind.coverage][ind.mapQ],
      split = "NM:i:"), function(x) x[2], mc.cores = num.cores)))
  ind.exactMatch <- which(editDistance == 0)
  ind.edit <- which(editDistance <= editDist.thold)
  
  if (exists("num.reads")) {
    cat("After filtering: found ", length(ind.edit), " module combinations in ",
        num.reads, " reads (", 
        round(length(ind.edit) * 100 / num.reads, digits = 2), "%)", sep = "")
  } else {
    cat("After filtering: found ", length(ind.edit), " module combinations.", 
        sep = "")
  }
  
  
  res.sam.filt <- 
    mclapply(res.sam.realn, function(x) x[ind.sam.realn][ind.coverage][
      ind.mapQ], mc.cores = num.cores)
  
  ###### 1: Distribution of modular variants
  aux.modDist <- table(factor(unlist(strsplit(res.sam.filt[["refID"]][
    ind.edit], split = ":")), levels = names(patterns.list)))
  mod.id2names <- 
    unlist(mclapply(seq_len(mod.tot), 
                    function(x) row.names(na.omit(patterns[x])), 
                    mc.cores = num.cores))
  names(mod.id2names) <- names(patterns.list)
  
  df.modDist <- 
    data.frame("module"  = factor(rep(names(patterns), mod.number), 
                                  levels = names(patterns)),
               "variant" = mod.id2names[names(aux.modDist)], 
               "counts"  = as.numeric(aux.modDist), row.names = NULL)
  
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, "_FreqVariantsPerModule_table.csv", sep = ""))
  cat("Exporting distribution of modular variants: \"", out.file,"\". \n", sep = "")
  write.csv(df.modDist, file = out.file)
  
  df.modDist <- 
    data.frame("module"  = factor(rep(names(patterns), mod.number), 
                                  levels = names(patterns)), 
               "variant" = mod.id2names[names(aux.modDist)],
               "vnumber" = unlist(mclapply(mod.number, seq_len, 
                                           mc.cores = num.cores), 
                                  use.names = FALSE), 
               "freq"    = as.numeric(aux.modDist) * 100 / length(ind.edit), 
               row.names = NULL)
  
  cat('Printing distribution of modular variants plot as \"', out.dir,
      out.filename, '_FreqVariantsPerModule_graph.pdf\". \n', sep = "")
  pl.modDist <- ggplot(df.modDist[df.modDist$freq != 100, ],
                       aes(x = vnumber, y = freq, fill = variant))
  pl.modDist <- pl.modDist +  geom_bar(stat = "identity", width = 0.7) + 
    facet_wrap(~module) + labs(x = "Modular variants", y = "Frequencies [%]") +
    ggtitle("ModSeq | Distribution of modular variants") + theme_bw() + 
    theme(axis.title = element_text(size = 14), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.text.y = element_text(size = 14),
          plot.title = element_text(face = "bold", size = 16),
          legend.title = element_blank())
  ggsave(paste(out.dir, out.filename, '_FreqVariantsPerModule_', 
               'graph.pdf', sep = ""), pl.modDist,  paper = "a4r", width = 12)

  
  ##### 2: Distribution of module combinations
  mod.dstr <- table(res.sam.filt[["refID"]][ind.edit])
  mod.dstr <- sort(mod.dstr, decreasing = TRUE)
  length(mod.dstr) <- mod.comb.len
  df.mod.dstr <- data.frame("x" = seq_len(mod.comb.len) * 100 / mod.comb.len, 
                            "y" = mod.dstr)
  
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, "_ModComb_distribution_graph.pdf", sep = ""))
  cat("Plotting distribution of module combinations: \"", out.file, "\". \n", 
      sep = "")
  
  PlotModuleCombinationsDistribution(data = df.mod.dstr, x.var = "x", 
                                     y.var = "y", out.file = out.file)
  
}
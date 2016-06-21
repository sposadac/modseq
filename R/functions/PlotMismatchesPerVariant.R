PlotMismatchesPerVariant <- function(substitution.modPos, ind, mod.len, var,  
                                     tot, patterns, filename, plotParams=NULL, 
                                     out.dir=NULL, num.cores=numeric(0)) {
  # plotParams
  # 1 ymax
  # 2 y-ticks
  # 3 width
  
  ## Whenever output directory is not specified, assumed to be the current 
  #  working directory
  if (is.null(out.dir)) {
    out.dir <- getwd()
  } 
  
  if (length(num.cores) == 0) {
    num.cores <- detectCores()
  }
  
  aux <- mclapply(
    strsplit(unlist(substitution.modPos, use.names = FALSE)[ind], split = ":"),
    function(x) c(x[2], x[3]), mc.cores = num.cores
    )
  
  if (length(aux) > 0) {
    
    aux <- unlist(aux)
    df <- data.frame(
      "pos"   = factor(as.numeric(aux[seq(1, length(aux), 2)]),
                       levels = seq(1, mod.len[var], 1)),
      "subst" = factor(aux[seq(2, length(aux), 2)], 
                       levels = c("A", "C", "G", "T", "N"))
      )
    
    if (is.null(plotParams)) {
      
      aux <- table(df$pos)
      
      plotParams <- numeric(3)
      plotParams[1] <- roundUp(max(aux))
      plotParams[2] <- plotParams[1] / 5 
      plotParams[3] <- max(mod.len[var] / 2, 6)
      
    }
    
    out.file <- 
      file.path(out.dir, 
                paste(filename, "_", var, "_mismatchPerPos.pdf", sep = ""))
    
    pl.base <- ggplot(df, aes(x = pos, fill = subst, colour = subst))
    pl <- pl.base + 
      stat_count(position = "stack", origin = -0.5, width = 0.6) + 
      scale_x_discrete(breaks = seq(1, mod.len[var], 1), 
                       labels = s2c(patterns[var]), 
                       name   = "Bases per position\n") +
      #scale_y_continuous(name   = "Mismatch counts", 
      #                   limits = c(0, as.numeric(plotParams[1]))) + 
      labs(y = "Counts\n") + 
      scale_fill_manual(breaks = c("A", "C", "G", "T", "N"), 
                        values = c("A" = "#F67771", "C" = "#FFB165", 
                                   "G" = "#1EBE7F", "T" = "#4B93FF", 
                                   "N" = "#7519FF"),
                        labels = c("A" = "A", "C" = "C", "G" = "G", "T" = "T", 
                                   "N" = "N")) + 
      scale_colour_manual(values = c("A" = "#F6746E", "C" = "#FFB165", 
                                     "G" = "#1EBE7F", "T" = "#4B93FF", 
                                     "N" = "#7519FF") , guide = FALSE) + 
      ggtitle(paste("ModSeq | Distribution of mismatches for variant", var)) + 
      theme_bw() +
      theme(legend.position = "top", legend.title = element_blank(), 
            plot.title = element_text(face = "bold", size = 16),
            text = element_text(size = 14))
    
    ggsave(out.file, pl,  paper = "a4r", width = as.numeric(plotParams[3]), 
           height = 3.5)
    
    # Normalized by the number of total reads that are mapped to a module 
    # combination containing the corresponding variant
    out.file <- 
      file.path(out.dir, 
                paste(filename, "_", var, "_mismatchPerPos_frequencies.pdf",
                      sep = ""))
                
    pl.norm <- pl.base + 
      stat_count(position = "stack", origin = -0.5, width = 0.6) +
      scale_x_discrete(breaks = seq(1, mod.len[var], 1), 
                       labels = s2c(patterns[var]), 
                       name = "\nBases per position") +
      scale_y_continuous(breaks = seq(0, tot[var], as.numeric(plotParams[2])),
                         labels = round(
                           seq(0, 1, as.numeric(plotParams[2]) / tot[var]) * 100, 
                           digit = 3), 
                         name = "Frequency [%]\n", 
                         limits = c(0, as.numeric(plotParams[1]))) +   
      scale_fill_manual(breaks = c("A", "C", "G", "T", "N"), 
                        values = c("A" = "#F67771", "C" = "#FFB165", 
                                   "G" = "#1EBE7F", "T" = "#4B93FF", 
                                   "N" = "#7519FF"), 
                        labels = c("A" = "A","C" = "C","G" = "G","T" = "T",
                                   "N" = "N")) + 
      scale_colour_manual(values = c("A" = "#F6746E", "C" = "#FFB165", 
                                     "G" = "#1EBE7F", "T" = "#4B93FF", 
                                     "N" = "#7519FF") , guide = FALSE) + 
      ggtitle(paste("ModSeq | Frequecies of mismatches for variant", var)) +
      theme_bw() +
      theme(legend.position = "top", legend.title = element_blank(), 
            plot.title = element_text(face = "bold", size = 16),
            text = element_text(size = 14))
    
    ggsave(out.file, pl.norm, paper = "a4r", width = as.numeric(plotParams[3]),
           height = 3.5)
    
    return(pl.norm)
    
    }

} 

roundUp <- function(x, accuracy=10) {
  return(ceiling(x/accuracy) * accuracy)
}
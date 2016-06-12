# library(seqinr)
PlotMismatchesPerVariant <- function(substitution.modPos, ind, mod.len, var, tot, 
                                     patterns, filename, subdirectory, 
                                     plotParams = NULL) {
  # plotParams
  # 1 ymax
  # 2 y-ticks
  # 3 width
  aux <- mclapply(strsplit(unlist(substitution.modPos)[ind], split = ":"), 
                  function(x) c(x[2], x[3]), mc.cores = num.cores)
  if (length(aux) > 0) {
    df <- 
      data.frame("pos" = factor(as.numeric(unlist(aux)[seq(1, length(unlist(aux)), 2)]), 
                                levels = seq(1, mod.len[var], 1)),
                 "subst" = factor(unlist(aux)[seq(2, length(unlist(aux)), 2)], 
                                  levels = c("A", "C", "G", "T", "N")))
    if (length(plotParams) == 0) {
      plotParams <- numeric(3)
      plotParams[1] <- signif(max(table(df$pos)), digit = 1)
      plotParams[2] <- max(50, signif(min(table(df$pos)), digit = 1))
      plotParams[3] <- max(mod.len[var] / 2, 2.5)
    }
    
    # normalized by the number of total reads that are mapped to a module 
    # combination containing the corresponding variant
    pl.base <- ggplot(df, aes(x = pos, fill = subst, colour = subst))
    pl <- pl.base + stat_bin(position = "stack", binwidth = 1, origin = -0.5, 
                        drop = FALSE, width = 0.6) + 
      scale_x_discrete(breaks = seq(1, mod.len[var], 1),
                       labels = s2c(patterns[var]), 
                       name   = "Bases per position") +
      scale_y_continuous(name   = "Mismatch counts", 
                         limits = c(0, as.numeric(plotParams[1]))) + 
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
    ggsave(paste('./OUTPUT/', subdirectory, "/", filename, "_", var, 
                 "_mismatchPerPos.pdf", sep = ""), pl, 
           width = as.numeric(plotParams[3]), height = 3.5)
    
    pl.norm <- pl.base + stat_bin(position = "stack", binwidth = 1, origin = -0.5, 
                             drop = FALSE, width = 0.6) +
      scale_x_discrete(breaks = seq(1, mod.len[var], 1), 
                       labels = s2c(patterns[var]), 
                       name = "Bases per position") +
      scale_y_continuous(breaks = seq(0, tot[var], as.numeric(plotParams[2])),
                         labels = round(seq(0, 1, as.numeric(plotParams[2]) / 
                                              tot[var]) * 100, digit = 2), 
                         name = "Mismatch frequencies [%]", 
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
    ggsave(paste('./OUTPUT/', subdirectory,"/", filename, "_", var, 
                 "_mismatchPerPos_frequencies.pdf", sep = ""), pl.norm, 
           width = as.numeric(plotParams[3]), height = 3.5)
    }
} 
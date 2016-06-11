PlotVariantCounts <- function(data, x.var, y.var, z.var, counts, labels,
                              modules.tot, modules.num, out.file) {
  
  ## Function arguments
  # data        Data frame containing number of hits per variant per search 
  #             round.
  # x.var       Identifier of the x-axis (i.e. different modules).
  # y.var       Identifier of the y-axis (i.e. frequencies).
  # z.var       Identifier of modular variants.
  # counts      Total number of hits per module per search round.
  # labels      Names of the different modules.
  # modules.tot Total nuumber of modules.
  # modules.num Number of variants per module.
  # out.file    Path to the output plot.
  
  ymax <- max(counts)
  breaks <- data[ , x.var]
  aux.data <- data.frame(
    x = seq_len(modules.tot + 1), y1 = counts, y2 = rep(0, modules.tot + 1),
    num = c(NA, modules.num))
  cols <- rep("gray", length(levels(data[ , z.var])))
    
  pl.patCounts <- 
    ggplot(data, aes_string(x = x.var, y = y.var, ymax = ymax * 1.1)) 
  
  pl.patCounts <- pl.patCounts + 
    geom_bar(stat = "identity", aes_string(fill = z.var, colour = z.var), #aes(color = "gray"), 
             alpha = 0.8) + labs(x = "\nSearch round", y = "Counts\n") + 
    scale_x_discrete(breaks = levels(breaks), labels = labels) + 
    scale_colour_manual(guide = FALSE, values = cols) +
    theme_bw() +
    ggtitle("ModSeq | Distribution of modular variants per search round") + 
    theme(text = element_text(size = 14),
          plot.title = element_text(face = "bold", size = 16), 
          legend.title = element_blank()) +
    geom_text(data = aux.data, 
              aes(label = y1, x = x, y = y1, ymax = max(y1) * 1.1), size = 5, 
              hjust = 0.5, vjust = -0.2, position = "stack") + 
    geom_text(data = aux.data,
              aes(label = num, x = x, y = y2, 
                  ymax = max(y1) * 1.1), size = 4, hjust = 0.5, 
              vjust = 1.1, position = "stack")
  
  ggsave(out.file, pl.patCounts, paper = "a4r", width = 12, height = 7) 
  
  return(pl.patCounts)
}
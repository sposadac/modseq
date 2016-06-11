PlotModuleCounts <- function(data, x.var, y.var, ymax, out.file, 
                             plot.label = "ModSeq", gls.ambiguity = TRUE) {
  
  ## Function arguments
  # data        data frame with two columns containing name of search rounds
  #             and counts
  # out.file    path and name of the file to save the plot
  # plot.label  optional argument
  # ambiguity   optional argument
  
  Nsearch <- "Ambiguous matches enabled"
  if (gls.ambiguity == FALSE) {
    Nsearch <- "Ambiguous matches disabled"
    cat("Ambiguous matches disabled for the pattern search (i.e., N only to ",
        "N). \n")
  }
  
  pl.patCounts <- 
    ggplot(data, aes_string(x = x.var, y = y.var, ymax = ymax * 1.1))

  pl.patCounts <-
    pl.patCounts + geom_bar(stat = "identity") + 
    labs(x = "\nSearch round", y = "Counts\n") +  
    ggtitle(paste("ModSeq | Module counts (", Nsearch,")", sep = "")) +
    theme_bw() + 
    theme(text = element_text(size = 14), 
          plot.title = element_text(face = "bold", size = 16)) + 
    annotate("text", x = Inf, y = -Inf, label = plot.label, hjust = 1.1, 
             vjust = -0.6, cex = 3) + 
    geom_text(aes_string(label = y.var), size = 5, hjust = 0.5, vjust = -0.2, 
              position = "stack")
  
  ggsave(filename = out.file, plot = pl.patCounts, paper = "a4r", width = 12,
         height = 7)
  
  return(pl.patCounts)
}
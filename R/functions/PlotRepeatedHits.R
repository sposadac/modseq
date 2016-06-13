library(scales)

PlotRepeatedHits <- function(data, out.file) {
  
  ## Function arguments
  # data      Table containing how many time each mapped read has been 
  #           assigned to different module combinations.
  # out.file  path and file name to output the plot
  
#   pdf(out.file)
#   plot(as.numeric(data), type = "h", xlab = "Reads", ylab = "Counts",
#        main = "Read mapping")
#   dev.off()
  
  df.rep <- as.data.frame(data)
  pl.rep <- ggplot(df.rep, aes(x = factor(Freq)))
  pl.rep <- pl.rep + geom_bar(width = 0.5) + 
    labs(x = "\nNumber of hits", y = "Number of reads\n") + 
    scale_y_continuous(trans = "log10",
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))) +
    ggtitle(paste("ModSeq | Repeated hits", sep = "")) + theme_bw() + 
    theme(text = element_text(size = 14), 
          plot.title = element_text(face = "bold", size = 16))
  
  ggsave(filename = out.file, plot = pl.rep, paper = "a4r", width = 6,
         height = 3.5)
  
  return(pl.rep)
}


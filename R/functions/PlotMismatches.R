PlotMismatches <- function(data, x.var, y.var, out.file) {
  
  ## Function arguments
  # data        data frame with two columns containing ids of modular variants
  #             and frequencies of mismatches normalized by their length
  
  pl.mismatch <- ggplot(data, aes_string(x = x.var, y = y.var))
  
  pl.mismatch <- pl.mismatch + 
    geom_bar(stat = "identity", position = "identity", width = 0.7) +
    labs(x = "\nModular variants", y = "Normalized mismatch frequency [%] \n") +
    ggtitle("ModSeq | Distribution of mismatches per modular variants") + 
    theme_bw() + 
    theme(plot.title = element_text(face = "bold", size = 16),
          text = element_text(size = 14))
  
  ggsave(out.file, pl.mismatch, paper = "a4r", width = 12) 
  
  return(pl.mismatch)
}
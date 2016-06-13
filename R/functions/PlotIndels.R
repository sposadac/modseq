PlotIndels <- function(data, x.var, y.var, z.var, out.file) {
  
  ## Function arguments
  # data        data frame with three columns containing ids of modular variants
  #             type (insertion or deletion) and frequencies of indels normalized 
  #             by their length
  
  pl.indel <- ggplot(data, aes_string(x = x.var, y = y.var, fill = z.var))
  
  pl.indel <- pl.indel + 
    geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
    labs(x = "\nModular variants", y = "Normalized indel frequency [%] \n") +
    scale_fill_discrete(breaks = c("Del", "In"), 
                        labels = c("Deletion", "Insertion"), 
                        guide = guide_legend(reverse = TRUE)) + 
    ggtitle("ModSeq | Distribution of indels per modular variants") + 
    theme_bw() + 
    theme(legend.position = "top", 
          plot.title = element_text(face = "bold", size = 16),
          text = element_text(size = 14))

  ggsave(out.file, pl.indel, paper = "a4r", width = 12)
  
  return(pl.indel)
}





PlotModuleCombinationsDistribution <- function(data, x.var, y.var, out.file) {
  
  ## Function arguments
  # data        Data frame containing number of hits per module combination.
  # x.var       Identifier of the x-axis.
  # y.var       Identifier of the y-axis.
  # out.file    Path to the output plot.
  
  pl.dstr <- ggplot(data, aes_string(x = x.var, y = y.var))
  
  pl.dstr <- pl.dstr + geom_point(na.rm = TRUE) +
    labs (x = "Module combinations [%]", y = "Frequency") + 
    ggtitle("ModSeq | Distribution of module combinations") + theme_bw() + 
    theme(plot.title = element_text(face = "bold", size = 16), 
          text = element_text(size = 14))
  
  ggsave(filename = out.file, plot = pl.dstr, paper = "a4r", width = 6)
  
  return(pl.dstr)
}

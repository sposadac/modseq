PlotVariantFrequencies <- function(data, x.var, y.var, z.var, facets,
                                            out.file) {
  
  ## Function arguments
  # data        Data frame containing number of hits per variant per search 
  #             round.
  # x.var       Identifier of the x-axis (i.e. different variants per module).
  # y.var       Identifier of the y-axis (i.e. frequencies).
  # z.var       Identifier of modular variants.
  # facets      Identifier of the facets (i.e. different modules).
  # out.file    Path to the output plot.
  
  pl.modDist <- ggplot(data, aes_string(x = x.var, y = y.var, fill = z.var))
  
  pl.modDist <- pl.modDist +  geom_bar(stat = "identity", width = 0.7) + 
    facet_wrap(facets) + labs(x = "Modular variants", y = "Frequencies [%]") +
    ggtitle("ModSeq | Distribution of modular variants") + 
    theme_bw() + 
    theme(axis.title = element_text(size = 14), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.text.y = element_text(size = 14),
          plot.title = element_text(face = "bold", size = 16),
          legend.title = element_blank())
  
  ggsave(out.file, pl.modDist,  paper = "a4r", width = 12)
  
  return(pl.modDist)
  
}
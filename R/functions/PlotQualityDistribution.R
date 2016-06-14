PlotQualityDistribution <- function(reads, out.filename, out.dir, title=NULL) {
  
  # Average quality per read
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, "_averageQualityDistribution_graph.pdf",
                    sep = ""))
  cat("Plotting the distribution of average quality scores per reads as: \"",
      out.file, "\".\n", sep = "")
  avQual <- alphabetScore(reads) / width(reads)
  avQual <- as.data.frame("x" = avQual)
  pl.quality <- ggplot(avQual, aes(x = avQual))
  
  pl.quality <- pl.quality + 
    geom_histogram(aes(y=..density..), binwidth=1, colour="black", 
                   fill="gray") +
    geom_density(alpha=.5, fill="#33CCFF", colour="#33CCFF") + 
    labs(x = "\nAverage quality score", y = "density\n") +
    ggtitle(paste("ModSeq | Distribution of the average quality scores per read |",
            title)) +
    theme_bw() + 
    theme(text = element_text(size = 14), 
          plot.title = element_text(face = "bold", size = 16))
  ggsave(out.file, pl.quality, paper = "a4r", width = 12) 
  
  return(pl.quality)
}
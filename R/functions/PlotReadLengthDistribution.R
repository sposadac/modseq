PlotReadLengthDistribution <- function(reads, out.filename, out.dir, title=NULL) {
  
  # Average quality per read
  read.len <- width(reads)
  
  if (length(unique(read.len)) != 1) {
    out.file <- 
      file.path(out.dir, 
                paste(out.filename, "_readLengthDistribution_graph.pdf",
                      sep = ""))
    cat("Plotting the read-length distribution as: \"", out.file, "\".\n", sep = "")
    
    read.len <- as.data.frame("x" = read.len)
    
    pl.readLen <- ggplot(read.len, aes(x = read.len))
    
    pl.readLen <- pl.readLen + 
      geom_histogram(aes(y=..density..), binwidth=1, colour="black", 
                     fill="gray") +
      labs(x = "\nRead length", y = "density\n") +
      ggtitle(paste("ModSeq | Read-length distribution |", title)) +
      theme_bw() + 
      theme(text = element_text(size = 14), 
            plot.title = element_text(face = "bold", size = 16))
    ggsave(out.file, pl.readLen, paper = "a4r", width = 12) 
    
    return(pl.readLen)
  } else {
    cat("Read-length: ", unique(read.len), " bp.\n", sep = "")
  }

}
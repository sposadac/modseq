plotReadLengthDistribution <- function(reads, out.filename, out.dir) {
  
  # Average quality per read
  out.file <- 
    file.path(out.dir, 
              paste(out.filename, "_readLengthDistribution_graph.pdf",
                    sep = ""))
  cat("Plotting the read-length distribution as: \"", out.file, "\".\n", sep = "")
  read.len <- width(reads)
  read.len <- as.data.frame("x" = read.len)
  pl.readLen <- ggplot(read.len, aes(x = read.len))
  
  pl.readLen <- pl.readLen + 
    geom_histogram(aes(y=..density..), binwidth=1, colour="black", 
                   fill="gray") +
    labs(x = "Read length", y = "density") +
    ggtitle("ModSeq | Read-length distribution") +
    theme_bw() + 
    theme(text = element_text(size = 14), 
          plot.title = element_text(face = "bold", size = 16))
  ggsave(out.file, pl.readLen, paper = "a4r", width = 12) 
  
  return(pl.readLen)
}
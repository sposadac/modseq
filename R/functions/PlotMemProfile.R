PlotMemProfile <- function(memList, memTrace, offsetMem, num.mod, out.file, 
                           out.dir) {
  
  out.file <- file.path(out.dir, paste("MemProf_", out.file, "_V2.pdf", sep=""))
  cat("Printing plot of memory consumption in: \"", out.file, "\".\n", sep = "")
  pdf(out.file, paper = "a4r", width = 12)
   
  ## Adding extra space to right margin of plot within frame
  par(mar=c(5, 4, 4, 6) + 0.1)
  
  ## Producing first plot and drawing its axis
  x <- seq(0, num.mod)
  y <- as.numeric(unlist(
    sapply(memList, function(x) strsplit(x, split = "Mb"), USE.NAMES=FALSE)))
  ylim1 <- min(y) %/% 10 * 10
  ylim2 <- max(y) %/% 10 * 10 + 10
  plot(x, y, pch=21, axes=FALSE, ylim=c(ylim1, ylim2), xlab="", type="p", 
       ylab="", col="black", main="Memory profile")
  lines(x, y, type="s", col="black")
  axis(2, ylim=c(ylim1,ylim2), col="black", las=1)  ## las=1 makes horizontal labels
  box()
  mtext("Size of the list (MB)", side=2, line=2.5)
  
  ## Adding Legend
  legend(x=1, y=ylim2, legend=c("object.size", "Ncells/Vcells"), 
         text.col=c("black","red"), pch=c(21,15), col=c("black","red"), 
         bty = "n")
  
  ## Allow a second plot on the same graph
  par(new=TRUE)
  
  ## Producing second plot and setting axis scale on right
  y <- memTrace[offsetMem:(num.mod + offsetMem)]
  ylim1 <- min(y) %/% 10 * 10
  ylim2 <- max(y) %/% 10 * 10 + 10
  plot(x, y, pch=15, xlab="", ylab="", ylim=c(ylim1, ylim2), axes=FALSE, 
       type="p", col="red")
  lines(x, y, col="red", type="s")
  
  ## A little farther (line=4) to make room for labels
  mtext("Memory consumption (MB)", side=4, col="red", line=4) 
  axis(4, ylim=c(ylim1, ylim2), col="red", col.axis="red", las=1)
  
  ## Draw the x axis
  axis(1, pretty(range(seq(0, 9)), 10))
  mtext("Search round", side=1, col="black", line=2.5) 
  
  dev.off()
}
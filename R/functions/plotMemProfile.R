plotMemProfile <- function(r1, r2, r3, r4, outname, path) {
  
  out.file <- paste(path, "MemProf_", outname, "_V2.pdf", sep="")
  cat("Printing plot of memory consumption in: \"", out.file, "\".\n", sep = "")
  pdf(out.file, paper = "a4r", width = 12)
   
  ## add extra space to right margin of plot within frame
  par(mar=c(5, 4, 4, 6) + 0.1)
  
  ## Plot first set of data and draw its axis
  #plot(seq(0,9), as.numeric(substr(memList,1,5)), pch=21, axes=FALSE,
  #     ylim=c(40,140), xlab="", type="p", ylab="", col="black", 
  #     main="Memory profile")
  plot(seq(0, num.mod),
       as.numeric(unlist(
         sapply(memList, function(x) strsplit(x, split = "Mb"), 
                USE.NAMES=FALSE))), pch=21, axes=FALSE, ylim=c(r1,r2), xlab="",
       type="p", ylab="",col="black", main="Memory profile")
  lines(seq(0,num.mod), 
        as.numeric(unlist(
          sapply(memList, function(x) strsplit(x, split = "Mb"), 
                 USE.NAMES=FALSE))), type="s", col = "black")
  axis(2, ylim=c(r1,r2),col="black",las=1)  ## las=1 makes horizontal labels
  box()
  mtext("Size of the list (MB)",side=2,line=2.5)
  
  ## Add Legend
  legend(x=1, y=r2, legend=c("object.size", "Ncells/Vcells"), 
         text.col=c("black","red"), pch=c(21,15), col=c("black","red"), 
         bty = "n")
  
  ## Allow a second plot on the same graph
  par(new=TRUE)
  
  ## Plot the second plot and put axis scale on right
  plot(seq(0,num.mod), Mem_trace[offset_mem:(num.mod+offset_mem)], pch=15,  
       xlab="", ylab="", ylim=c(r3,r4), axes=FALSE, type="p", col="red")
  lines(seq(0, num.mod), Mem_trace[offset_mem:(num.mod + offset_mem)], 
        col="red", type="s")
  
  ## A little farther (line=4) to make room for labels
  mtext("Memory consumption (MB)", side=4, col="red", line=4) 
  axis(4, ylim=c(r3, r4), col="red", col.axis="red", las=1)
  
  ## Draw the x axis
  axis(1, pretty(range(seq(0, 9)), 10))
  mtext("Search round", side=1, col="black", line=2.5) 
  
  dev.off()
}
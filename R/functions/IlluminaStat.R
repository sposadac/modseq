# Author: S.Schmitt. Modified by S.Posada.
library(grid)
library(gridExtra)
IlluminaStat <- function(in.dir, pattern, SV_plotinputframe, SV_qcwd1, SV_qcwd2, 
                         SV_meanread1, SV_meanread2, SV_meanread3, SV_meanread4,
                         pdfname, type, out.dir, sample=FALSE, in.dirRaw=NULL,
                         forward.filename=NULL, reverse.filename=NULL) {
  
  if (type == "unpaired") {
    
    # Read-in Illumina-files
    SV_illuminaqcF1raw <- 
      qa(dirPath=in.dirRaw, pattern=paste(forward.filename, ".fastq", sep=""),
         type="fastq")
    SV_illuminaqcF2raw <- 
      qa(dirPath=in.dirRaw, pattern=paste(reverse.filename, ".fastq", sep=""),
         type="fastq")
    
    SV_illuminaqcF1 <- 
      qa(dirPath=in.dir, pattern=paste(pattern, "_1_nQtrim.fastq", sep=""), 
         type="fastq")
    SV_illuminaqcF2 <- 
      qa(dirPath=in.dir, pattern=paste(pattern, "_2_nQtrim.fastq", sep=""),
         type="fastq")
    
    # Make graphs and print as PDF
    pdf(file.path(out.dir, pdfname), paper = "a4r", width = 12)
    print("... plot counts")
    plot(
      ShortRead:::.plotReadCount(
        rbind(SV_illuminaqcF1raw, SV_illuminaqcF2raw, SV_illuminaqcF1,
              SV_illuminaqcF2))
      )
    plot(
      ShortRead:::.plotNucleotideCount(
        rbind(SV_illuminaqcF1raw, SV_illuminaqcF2raw, SV_illuminaqcF1, 
              SV_illuminaqcF2)))
    print("... plot quality score density")
    SV_qcsubF1raw <- SV_illuminaqcF1raw[["readQualityScore"]];
    SV_qcsubF2raw <- SV_illuminaqcF2raw[["readQualityScore"]];
    SV_qcsubF1 <- SV_illuminaqcF1[["readQualityScore"]];
    SV_qcsubF2 <- SV_illuminaqcF2[["readQualityScore"]];
    par(mfrow = c(2,1));
    plot(
      ShortRead:::.plotReadQuality(
        rbind(SV_qcsubF1raw[SV_qcsubF1raw$type=="read", ], 
              SV_qcsubF1[SV_qcsubF1$type=="read", ]))
      )
    plot(
      ShortRead:::.plotReadQuality(
        rbind(SV_qcsubF2raw[SV_qcsubF2raw$type=="read", ],
              SV_qcsubF2[SV_qcsubF2$type=="read", ]))
      )
    print("... plot cycle vs. basecalls and quality")
    SV_qcsubF1raw <- SV_illuminaqcF1raw[["perCycle"]]
    SV_qcsubF2raw <- SV_illuminaqcF2raw[["perCycle"]]
    SV_qcsubF1 <- SV_illuminaqcF1[["perCycle"]]
    SV_qcsubF2 <- SV_illuminaqcF2[["perCycle"]]
    par(mfrow = c(2,1))
    plot(
      ShortRead:::.plotCycleBaseCall(
        rbind(SV_qcsubF1raw$baseCall, SV_qcsubF1$baseCall))
      )
    plot(
      ShortRead:::.plotCycleBaseCall(rbind(
        SV_qcsubF2raw$baseCall,SV_qcsubF2$baseCall))
      )
    par(mfrow = c(2,1));
    plot(
      ShortRead:::.plotCycleQuality(
        rbind(SV_qcsubF1raw$quality, SV_qcsubF1$quality))
      )
    plot(ShortRead:::.plotCycleQuality(
      rbind(SV_qcsubF2raw$quality, SV_qcsubF2$quality))
      )
    print("... plot read length distributions")
    df <- data.frame(x = SV_plotinputframe[, 3], 
                     y = SV_plotinputframe[, 4])
    SG_hist <- ggplot(df, aes(x = x))
    SG_hist <- SG_hist + 
      geom_density(binwidth=1, fill="red", alpha=.1, size=.3, color="red") + 
      geom_density(binwidth=1, fill="blue", alpha=.1, size=.3, aes(x = y),
                   color="blue") + theme_bw() +
      scale_x_continuous("Read length") + scale_y_continuous("Counts")
    SG_hist <- SG_hist + geom_vline(xintercept=SV_meanread1, color="black",
                                    linetype="dashed", size=.2)
    SG_hist <- SG_hist + geom_vline(xintercept=SV_meanread2, color="black",
                                    linetype="dashed", size=.2) 
    SG_hist <- SG_hist + geom_vline(xintercept=SV_meanread3, color="red", 
                                    linetype="dashed", size=.2) 
    SG_hist <- SG_hist + geom_vline(xintercept=SV_meanread4, color="blue",
                                    linetype="dashed", size=.2)
    plot(SG_hist)
    print("... plot read length vs. quality")
    par(mfrow = c(1,2))
    smoothScatter(
      SV_qcwd1, xlim = c(0, max(SV_qcwd1, na.rm=TRUE)), xlab="Read width", 
      ylab="Quality", main="Forward trimmed-reads: Width vs. Mean Quality", 
      colramp=colorRampPalette(c("white", "darkblue")), nrpoints = 0)
    smoothScatter(
      SV_qcwd2, xlim = c(0, max(SV_qcwd2, na.rm=TRUE)), xlab="Read width", 
      ylab="Quality", main="Reverse trimmed-reads: Width vs. Mean Quality", 
      colramp=colorRampPalette(c("white", "darkblue")), nrpoints = 0)
    print("... plot sequence distribution")
    SV_qcsubF1 <- SV_illuminaqcF1[["sequenceDistribution"]];
    SV_qcsubF2 <- SV_illuminaqcF2[["sequenceDistribution"]];
    plot(
      ShortRead:::.plotReadOccurrences(
        rbind(SV_qcsubF1[SV_qcsubF1$type=="read", ], 
              SV_qcsubF2[SV_qcsubF2$type=="read", ]), cex=.5))
    #SV_mostread <- ShortRead:::.freqSequences(
    #  rbind(SV_illuminaqcF1,SV_illuminaqcF2), "read")
    
  } else if (type == "single") {
    
    # Read in Illumina-files
    SV_illuminaqcF1 <- 
      qa(dirPath=in.dirRaw, pattern=paste(forward.filename, ".fastq", sep=""),
         type="fastq")
    SV_illuminaqcF1trim <- 
      qa(dirPath=in.dir, pattern=paste(pattern, "_1_nQtrim.fastq", sep=""), 
         type="fastq")
    
    # Make graphs and print as PDF
    pdf(file.path(out.dir, pdfname), paper = "a4r", width = 12)
    print("... plot counts")
    plot(
      ShortRead:::.plotReadCount(
        rbind(SV_illuminaqcF1, SV_illuminaqcF1trim))
      )
    plot(
      ShortRead:::.plotNucleotideCount(
        rbind(SV_illuminaqcF1, SV_illuminaqcF1trim))
      )
    print("... plot quality score density")
    SV_qcsubF1 <- SV_illuminaqcF1[["readQualityScore"]]
    SV_qcsubF1trim <- SV_illuminaqcF1trim[["readQualityScore"]]
    plot(
      ShortRead:::.plotReadQuality(
        rbind(SV_qcsubF1[SV_qcsubF1$type=="read", ], 
              SV_qcsubF1trim[SV_qcsubF1trim$type=="read", ]))
      )
    print("... plot cycle vs. basecalls and quality");
    SV_qcsubF1 <- SV_illuminaqcF1[["perCycle"]];
    SV_qcsubF1trim <- SV_illuminaqcF1trim[["perCycle"]];
    plot(
      ShortRead:::.plotCycleBaseCall(
        rbind(SV_qcsubF1$baseCall, SV_qcsubF1trim$baseCall))
      )
    plot(
      ShortRead:::.plotCycleQuality(
        rbind(SV_qcsubF1$quality, SV_qcsubF1trim$quality))
      )
    print("... plot read length distributions")
    df <- data.frame(x = SV_plotinputframe[, 1], 
                     y = SV_plotinputframe[, 2])
    SG_hist <- ggplot(df, aes(x = x))
    SG_hist <- SG_hist + 
      geom_density(binwidth=1, fill="blue", alpha=.1, size=.3, color="red") + 
      geom_density(binwidth=1, fill="red", alpha=.1, size=.3, aes(x = y), 
                   color="blue") + theme_bw() + 
      scale_x_continuous("Read length") + scale_y_continuous("Counts")
    SG_hist <- SG_hist + geom_vline(xintercept=SV_meanread1, color="blue",
                                    linetype="dashed", size=.2)
    SG_hist <- SG_hist + geom_vline(xintercept=SV_meanread2, color="red",
                                    linetype="dashed", size=.2)
    plot(SG_hist)
    par(mfrow = c(1,2));
    print("... plot read length vs. quality");
    smoothScatter(
      SV_qcwd1, xlim = c(0, max(SV_qcwd1, na.rm=TRUE)), xlab="Read width",
      ylab="Quality", main="Raw reads: Width vs. Mean Quality", 
      colramp=colorRampPalette(c("white", "darkblue")), nrpoints = 0)
    smoothScatter(
      SV_qcwd2, xlim = c(0, max(SV_qcwd2, na.rm=TRUE)), xlab="Read width", 
      ylab="Quality", main="Trimmed reads: Width vs. Mean Quality", 
      colramp=colorRampPalette(c("white", "darkblue")), nrpoints = 0)
    print("... plot sequence distribution")
    SV_qcsubF1 <- SV_illuminaqcF1[["sequenceDistribution"]];
    SV_qcsubF1trim <- SV_illuminaqcF1trim[["sequenceDistribution"]];
    plot(
      ShortRead:::.plotReadOccurrences(
        rbind(SV_qcsubF1[SV_qcsubF1$type=="read", ], 
              SV_qcsubF1trim[SV_qcsubF1trim$type=="read", ]), cex=.5)
      )
    #SV_mostread <- ShortRead:::.freqSequences(SV_illuminaqcF1trim, "read")
    
  } else if(type == "paired") {
    
    # Read in Illumina-files
    SV_illuminaqcF1 <- 
      qa(dirPath=in.dirRaw, pattern=paste(forward.filename, ".fastq", sep=""),
         type="fastq", sample = sample)
    SV_illuminaqc <- 
      qa(dirPath= in.dir, pattern=paste(pattern, "_PANDAseq.fastq", sep= ""), 
         type="fastq", sample = sample)
    
    # Make graphs and print as PDF
    pdf(file.path(out.dir, pdfname), paper = "a4r", width = 12)
    print("... plot counts")
    plot(ShortRead:::.plotReadCount(rbind(SV_illuminaqcF1, SV_illuminaqc)))
    plot(ShortRead:::.plotNucleotideCount(rbind(SV_illuminaqcF1, SV_illuminaqc)))
    print("... plot quality score density")
    SV_qcsubF1 <- SV_illuminaqcF1[["readQualityScore"]]
    SV_qcsub <- SV_illuminaqc[["readQualityScore"]]
    plot(
      ShortRead:::.plotReadQuality(
        rbind(SV_qcsubF1[SV_qcsubF1$type=="read", ], 
              SV_qcsub[SV_qcsub$type=="read", ]))
      )
    print("... plot cycle vs. basecalls and quality");
    SV_qcsubF1 <- SV_illuminaqcF1[["perCycle"]];
    SV_qcsub <- SV_illuminaqc[["perCycle"]];
    plot(
      ShortRead:::.plotCycleBaseCall(
        rbind(SV_qcsubF1$baseCall, SV_qcsub$baseCall))
      )
    plot(
      ShortRead:::.plotCycleQuality(
        rbind(SV_qcsubF1$quality, SV_qcsub$quality))
      )
    print("... plot read length distributions")
    df <- data.frame(x = SV_plotinputframe[, 1], 
                     y = SV_plotinputframe[, 2],
                     z = SV_plotinputframe[, 3])
    SG_hist <- ggplot(df, aes(x = x))
    SG_hist <- SG_hist +
      geom_density(binwidth=1, fill="blue", alpha=.1, size=.3, color="blue") + 
      geom_density(binwidth=1, fill="darkgrey", alpha=.1, size=.3, aes(x = y), 
                   color="darkgrey") + 
      geom_density(binwidth=1, fill="red", alpha=.1, size=.3, aes(x = z), 
                   color="red") + theme_bw() +
       scale_x_continuous("Read length") + scale_y_continuous("Counts")
    SG_hist <- SG_hist + geom_vline(xintercept=SV_meanread1, color="blue",
                                    linetype="dashed", size=.2)
    SG_hist <- SG_hist + geom_vline(xintercept=SV_meanread2, color="darkgrey",
                                    linetype="dashed", size=.2)
    SG_hist <- SG_hist + geom_vline(xintercept=SV_meanread3, color="red",
                                    linetype="dashed", size=.2)
    plot(SG_hist)
    par(mfrow = c(1,2))
    print("... plot read length vs. quality");
    smoothScatter(
      SV_qcwd1, xlim = c(0, max(SV_qcwd1, na.rm=TRUE)), xlab="Read width",
      ylab="Quality", main="Forward trimmed reads: Width vs. Mean Quality", 
      colramp=colorRampPalette(c("white", "darkblue")), nrpoints = 0)
    smoothScatter(
      SV_qcwd2, xlim = c(0, max(SV_qcwd2, na.rm=TRUE)), xlab="Read width", 
      ylab="Quality", main="Assembled trimmed reads: Width vs. Mean Quality", 
      colramp=colorRampPalette(c("white", "darkblue")), nrpoints = 0)
    print("... plot sequence distribution")
    SV_qcsubF1 <- SV_illuminaqcF1[["sequenceDistribution"]]
    SV_qcsub <- SV_illuminaqc[["sequenceDistribution"]]
    plot(
      ShortRead:::.plotReadOccurrences(
        rbind(SV_qcsubF1[SV_qcsubF1$type=="read", ], 
              SV_qcsub[SV_qcsub$type=="read",]), cex=.5))
    #SV_mostread <- ShortRead:::.freqSequences(SV_illuminaqc, "read")
  }
  
  #par(mfrow = c(1,1));
  #colnames(SV_mostread)[1] <- "TOP 20 Sequences in Illumina-Data";
  #grid.arrange(tableGrob(SV_mostread, gpar.coretext = gpar(fontsize=4)));
  print("... write to PDF");
  dev.off();
  print("... done!");
}
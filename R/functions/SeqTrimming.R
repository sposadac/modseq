SeqTrimming <- function(cigarOps, cigarOpLengths, seq){
  # cigarOps - list
  # cigarOpLengths - list
  # seq - character
  len <- length(cigarOps)
  if (cigarOps[1] == "S" & cigarOps[len] == "S") {
    if (sum(cigarOps == "D") > 0) {
      end <- sum(head(cigarOpLengths[(cigarOps != "D")], n = -1))
    } else {
      end <- sum(head(cigarOpLengths, n = -1))
    }
    seq.trim <- subseq(seq, start = (cigarOpLengths[1] + 1), end = end)
  } else if (cigarOps[1] == "S") {
    if (sum(cigarOps == "D") > 0) {
      end <- sum(cigarOpLengths[(cigarOps != "D")])
    } else {
      end <- sum(cigarOpLengths)
    }
    seq.trim <- subseq(seq, start = (cigarOpLengths[1] + 1), end = end)
  } else if (cigarOps[len] == "S") {
    if (sum(cigarOps == "D") > 0) {
      end <- sum(head(cigarOpLengths[(cigarOps != "D")], n = -1))
    } else {
      end <- sum(head(cigarOpLengths, n = -1))
    }
    seq.trim <- subseq(seq, start = 1, end = end)
  } else{
    seq.trim <- seq
  }
  return(seq.trim)
}

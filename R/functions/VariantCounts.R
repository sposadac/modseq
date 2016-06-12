VariantCounts <- function(patterns, num.reads, counts, labels, file.prefix, 
                          out.dir, gls.direction = "f", modseq.dir = NULL, 
                          num.cores = numeric(0)) {
  
  ## Function arguments
  # patterns      Data frame containing the sequences of the various variants
  #               per module. Columns correspond to different modules and rows 
  #               correspond to variants identifiers. Alternatively, it can 
  #               represent the path to the csv file.
  # num.reads     Number of input reads.
  # counts        Number of hits per search round.
  # labels        Identifiers of the search rounds.
  # file.prefix   Prefix for input as well as output files. Number of hits per 
  #               search round are read from output files containing this prefix.
  # out.dir       Path to output files.
  # gls.direction (optional) direction in wich the gls was performed ("f" or "r").
  # modseq.dir    (optional) modseq's head directory, by default current working
  #               directory.
  # num.cores     (optional) number of cores available for performing parallel 
  #               tasks.
  
  ## Whenever modseq path is not specified, assumed to be the current working
  #  directory
  if (is.null(modseq.dir)) {
    modseq.dir <- getwd()
    warning("Object \'modseq.dir\' not found, set to: \"", modseq.dir, "\".")
  } 
  
  if (!exists("PlotVariantCountsPerRound", mode = "function")) {
    source(file.path(modseq.dir, "R/functions/PlotVariantCountsPerRound.R"))
  }
  source(file.path(modseq.dir, "R/functions/GetVariantsNames.R"))
  
  if (length(num.cores) == 0) {
    num.cores <- detectCores()
  }
  
  if (!is.data.frame(patterns)) {
    source(file.path(modseq.dir, "R/functions/LoadModuleTable.R"))
    if (is.character(patterns)) {
      patterns <- LoadModuleTable(in.modulesDir = patterns) 
    } else {
      stop("Invalid data type for \'patterns\', \"", class(patterns), "\".\n")
    }
  }
  mod.tot <- ncol(patterns)  
  mod.number <- sapply(patterns, function(x) length(na.omit(x)))
  mod.number.cumsum <- cumsum(c(1, mod.number))
  mod.empty <- 
    unlist(mapply(function(x,i) if (sum(na.omit(x) == "") > 0) return(i),
                  x = patterns, i = seq_len(mod.tot)))
  
  if (gls.direction == "f") {
    mod.comb.cumprod <- cumprod(mod.number)
  } else if (gls.direction == "r") {
    mod.comb.cumprod <- rev(cumprod(rev(mod.number)))
  } 
  
  var.counts <- data.frame(
    matrix(ncol = 3, dimnames = list(c(""), c("module", "variant", "freq")))
  )
  var.counts[1, ] <- c(0, NA, num.reads)
  var.counts[2:(sum(mod.number) + 1), 1] <- rep(seq(mod.tot), mod.number)
  var.counts[2:(sum(mod.number) + 1), 2] <- 
    GetVariantsNames(patterns = patterns, modules.tot = mod.tot, 
                     modseq.dir = modseq.dir, num.cores = num.cores)
  
  # Reading number of hits per modular variants from csv files
  for (i in seq_len(mod.tot)) {
    in.file <- file.path(
      out.dir, paste(file.prefix, '_round_', i, '_of_', mod.tot,
                     '.csv', sep = ""))
    mod.aux <- data.frame(
      read.csv(in.file, nrow = mod.comb.cumprod[i], header = TRUE, row.names = 1))
    var.counts[(mod.number.cumsum[i] + 1):(mod.number.cumsum[i] + mod.number[i]), 3] <-
      unlist(mclapply(
        seq_len(mod.number[i]), 
        function(j) sum(mod.aux[seq(j, mod.comb.cumprod[i], mod.number[i]), 1]),
        mc.cores = num.cores))
  }
  rm(mod.aux)
  
  # Defining modules as factors 
  var.counts$module <- as.factor(var.counts$module)
  # Reversing order of modules 
  if (gls.direction == "r") {  
    var.counts$module <- 
      factor(var.counts$module, 
             levels = c(levels(var.counts$module)[1], 
                        rev(levels(var.counts$module)[2:(mod.tot+1)])))
    mod.number <- rev(mod.number)
  } 
  
  # Defining variants as factors
  var.counts$variant[
    unlist(mclapply(names(which(table(var.counts$module) == 1)),
                    function(x) which(x == var.counts$module), 
                    mc.cores = num.cores))] <- NA
  var.counts$variant <- as.factor(var.counts$variant)
  
  # Plotting number of hits per modular variant per search round
  out.file <- file.path(
    out.dir, paste(file.prefix, "_modularVariantsDistribution_graph.pdf", 
                   sep = ""))
  cat("Printing modular-variant counts plot as: \"", out.file, "\".\n", 
      sep = "")
  PlotVariantCountsPerRound(
    data=var.counts, x.var="module", y.var="freq", z.var="variant", 
    counts=counts, labels=labels, modules.tot=mod.tot, modules.num=mod.number,
    out.file=out.file
    )
    
  # Writing csv file containing number of hits per modular variants per search round
  var.counts$module <- c("0", rep(names(patterns), mod.number))
  out.file <- 
    file.path(out.dir, paste(file.prefix, "_modularVariants", 
                             "Distribution_table.csv", sep = ""))
  cat("Exporting distribution of modular variants per search round as table:",
      " \"", out.file, "\". \n", sep = "")
  write.csv(var.counts, out.file)
  
}
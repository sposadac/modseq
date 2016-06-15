ModuleCombinationsGen <- function(modules.filename, pattern=character(0), 
                                  out.file=NULL, in.modDir=NULL, modseq.dir=NULL, 
                                  num.cores=numeric(0)) {
  
  ## Function arguments
  # modules.filename  Character. Name of the file containing table of modules.
  # patterns          (optional) Data frame containing the sequences of the 
  #                   various variants per module. Columns correspond to 
  #                   different modules and rows correspond to variants 
  #                   identifiers. 
  # in.modDir         (optional) Expected when patterns or out.file are not 
  #                   provided.
  # modseq.dir        (optional) Expected when patterns are not provided.
  # num.cores         (optional) number of cores available for performing 
  #                   parallel tasks.
  
  ## Whenever input or output directories are not specified, assumed to be 
  #  the current working directory
  if (is.null(modseq.dir)) {
    modseq.dir <- getwd()
    warning("Object \'modseq.dir\' not found, set to: \"", modseq.dir, "\".")
  } 
  if (is.null(in.modDir)) {
    in.modDir <- getwd()
    warning("Object \'in.modDir\' not found, set to: \"", in.modDir, "\".")
  }
  
  if (length(num.cores) == 0) {
    num.cores <- detectCores()
  }
  
  ## Loading module table
  if (!is.data.frame(patterns)) {
    if (!existsFunction("LoadModuleTable")) {
      source(file.path(modseq.dir, "R/functions/LoadModuleTable.R"))
    }
    pattern <- LoadModuleTable(in.modulesDir = in.modDir,
                                modules.filename = modules.filename) 
  }
  
  Modules.modComb <- expand.grid(mclapply(as.list(pattern), na.omit))
  Modules.modComb <- 
    by(data = Modules.modComb, INDICES = c(1:nrow(Modules.modComb)), 
       FUN = function(x) paste(unlist(x), collapse = ""))
  Modules.modComb <- DNAStringSet(as.character(Modules.modComb))
  aux.len <- lapply(pattern, function(x) length(na.omit(x)))
  
  if (sum(unlist(aux.len)) != nrow(pattern)) {   
    rownames <- tolower(names(unlist(lapply(aux.len, function(x) seq(1:x)))))
    rownames <- 
      expand.grid(split(rownames, factor(rep(names(pattern), aux.len), 
                                          levels = names(pattern))))
  } else {
    rownames <- 
      expand.grid(split(row.names(pattern), factor(rep(names(pattern), aux.len), 
                                                   levels = names(pattern))))
  }
  rownames <- as.matrix(rownames)
  rownames <- 
    unlist(mclapply(c(1:nrow(rownames)), 
                    function(x) paste(rownames[x,], collapse = ":"), 
                    mc.cores = num.cores))  
  names(Modules.modComb) <- rownames
  
  if (is.null(out.file)) {
    if (is.null(in.modDir)) {
      in.modDir <- getwd() 
    }
    out.file <- file.path(in.modDir, 
                          paste(modules.filename, "_modComb.fasta", sep = ""))
  }
  
  cat("Writing module combinations as fasta file: \"", out.file, "\" ...\n",
      sep = "")
  writeFasta(object = Modules.modComb, file = out.file)
}

# moduleCombinationsGen <- function(modules.filename, pattern = character(0), mode = c('RData','fasta')){
#   mode <- match.arg(mode)
#   rownames <- as.data.frame(t(sapply(row.names(pattern), FUN=function(x) sapply(pattern[x,], FUN=function(cell) {if(is.na(cell) == FALSE) {return(x)} else {return(NA)}}))))
#   rownames <- expand.grid(mclapply(as.list(rownames), na.omit))
#   rownames <- by(data=rownames, INDICES=c(1:nrow(rownames)), FUN=function(x) paste(unlist(x), collapse=":"))
#   if (mode == 'fasta') {
#     writeFasta(Modules.modComb, paste('./MODULES/', modules.filename,'_modComb.fasta', sep=""))
#   }else{
#     save(Modules.modComb, file=paste('./MODULES/', modules.filename,'_modComb.RData', sep=""))
#   } 


ModuleCombinationsGen <- function(modules.filename, pattern = character(0), 
                                  in.modDir = NULL, modseq.dir = "", 
                                  num.cores = numeric(0)) {
  
  if (length(num.cores) == 0) {
    num.cores <- detectCores()
  }
  if (length(pattern) == 0) {
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
  
  out.file <- file.path(in.modDir, 
                        paste(modules.filename, '_modComb.fasta', sep = ""))
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


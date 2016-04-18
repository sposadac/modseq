ModuleCombinationsGen <- function(module.filename, pattern = character(0), 
                                  num.cores = numeric(0), in.modDir = NULL) {
  if (length(num.cores) == 0) {
    num.cores <- detectCores()
  }
  if (length(pattern) == 0) { 
    pattern <- 
      data.frame(
        read.csv(file.path(in.modDir, paste(module.filename, '.csv', sep = "")),
                 header = TRUE, sep = ",", dec = ".")
        )
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
  writeFasta(object = Modules.modComb, 
             file = file.path(in.modDir,
                              paste(module.filename, '_modComb.fasta', sep = "")))
}

# moduleCombinationsGen <- function(module.filename, pattern = character(0), mode = c('RData','fasta')){
#   mode <- match.arg(mode)
#   rownames <- as.data.frame(t(sapply(row.names(pattern), FUN=function(x) sapply(pattern[x,], FUN=function(cell) {if(is.na(cell) == FALSE) {return(x)} else {return(NA)}}))))
#   rownames <- expand.grid(mclapply(as.list(rownames), na.omit))
#   rownames <- by(data=rownames, INDICES=c(1:nrow(rownames)), FUN=function(x) paste(unlist(x), collapse=":"))
#   if (mode == 'fasta') {
#     writeFasta(Modules.modComb, paste('./MODULES/',module.filename,'_modComb.fasta', sep=""))
#   }else{
#     save(Modules.modComb, file=paste('./MODULES/',module.filename,'_modComb.RData', sep=""))
#   } 


### Running mode
if (!exists("run")) {
  warning("Unspecified running mode. \nObject \'run\' is set to vector of ",
          "zeros.")
  run <- rep(0, 5)
}
if (length(run) > 5) {
  warning("Too many options for running mode, object \'run\' has ", 
          length(run), " entries.\nObject \'run\' must have 5 entries only.")
  run <- run[1:5]
}
if (length(run) <= 0 || length(run) < 5) { 
  warning("Too few options for running mode, object \'run\' has ", 
          length(run), " entries.\nObject \'run\' must have 5 entries.")
  aux.run <- run
  run <- rep(0, 5)
  run[seq_along(aux.run)] <- aux.run
}
if (run[1] != 0 && run[1] != 1) {
  stop("Invalid entry for run[1], \"", run[1], "\".\nQuality trimming mode:", 
       " either 0-off or 1-on. \n")
}
if (run[2] != 0 && run[2] != 1) { 
  stop("invalid entry for run[2], \"", run[2], "\".\nPaired-end read assembly", 
       " mode: either 0-off or 1-on. \n")
}
if (run[3] != 0 && run[3] != 1) {
  stop("invalid entry for run[3], \"", run[3], "\".\nRead mapping mode: ", 
       "either 0-off or 1-on. \n")
}
if (run[4] != 0 && run[4] != 1) {
  stop("invalid entry for run[4], \"", run[4], "\".\nPost-analysis mode: ", 
       "either 0-off or 1-on. \n")
}
if (run[5] != 0 && run[5] != 1) {
  stop("invalid entry for run[5], \"", run[5], "\".\nVariant calling mode: ", 
       "either 0-off or 1-on. \n")
}

### Sequencing mode
if (!exists("seq.mode")) {
  warning("Object \'seq.mode\' not found, set to default value: \"PE\".")
  seq.mode <- "PE"
} else if (!is.character(seq.mode)) {
  warning("Invalid data type for \'seq.mode\', \"", class(seq.mode), "\".", 
          "\nObject \'seq.mode\' is set to \"PE\".")
  seq.mode <- "PE"
} else if (seq.mode != "PE" && seq.mode != "SE") {
  warning("Invalid entry for sequencing mode.\nObject \'seq.mode\' should be", 
          " either \"PE\" or \"SE\", \"", seq.mode, "\".\nObject \'seq.mode\'", 
          " is set to \"PE\".")
  seq.mode <- "PE"
}
if (seq.mode == "SE") {
  rm(forward.filename)
  rm(reverse.filename)
  paired.file <- character(0)
  paired.flag <- as.integer(0)
  run[2] <- as.integer(0)
} else if (seq.mode == "PE") {
  rm(in.filename)
}

#### Input files
## Raw sequencing data - input directory
if (!exists("in.seqDir")) {
  warning("Object \'in.seqDir\' not found, set to default value: \"",
          wdir, "/data/raw/\".")
  in.seqDir <- file.path(wdir, "data/raw/")
} else if (!is.character(in.seqDir)) {
  stop("Invalid data type for \'in.seqDir\', \"", class(in.seqDir), "\".")
} else if (in.seqDir == "" || in.seqDir == " ") {
  in.seqDir <- file.path(wdir)
} 
if (!file.exists(in.seqDir)) {
  stop ("Invalid directory for raw sequencing data files, \"", in.seqDir, 
        "\". \n")
}

## Raw sequencing data - file names
if (seq.mode == "SE"){
  if (exists("in.filename")) {
    aux <- file.path(in.seqDir, paste(in.filename, '.fastq', sep = ""))
    if (!is.character(in.filename)) {
      stop("Invalid data type for \'in.filename\', \"",  class(in.filename), 
           "\".\nObject \'in.filename\' should be of type \"character\".\n")
    } else if (!file.exists(aux)) {
      stop("File \"", aux, "\" not found.\n")
    }
  } else {
    stop("Object \'in.filename\' not found. \n")
  }
} else if (seq.mode == "PE") {
  if (exists("forward.filename") && exists("reverse.filename")) {
    if (!is.character(forward.filename) || !is.character(reverse.filename)) {
      stop("Sequencing data input: file name should be of type character. \n")
    } else {
      aux <- file.path(in.seqDir, paste(forward.filename, '.fastq', sep = ""))
      if (!file.exists(aux)) {
        stop("File \"", aux, "\" not found. \n")
      }
      aux <- file.path(in.seqDir, paste(reverse.filename, '.fastq', sep = ""))
      if (!file.exists(aux)) {
        stop("File \"", aux, "\" not found. \n")
      }
    }
  } else {
    if (!exists("forward.filename")) {
      stop("Object \'forward.filename\' not found. \n")
    }
    if (!exists("reverse.filename")) {
      stop("Object \'reverse.filename\' not found. \n")
    }
  }
}

# Module-table - input directory
if (!exists("in.modDir")) {
  warning("Object \'in.modDir\' not found, set to default value: \"",
          wdir, "/data/modules/\".")
  in.modDir <- file.path(wdir, "data/modules/")
} else if (!is.character(in.modDir)) {
  stop("Invalid data type for \'in.modDir\', \"", class(in.modDir), "\".")
} else if (in.modDir == "" || in.modDir == " ") {
  in.modDir <- file.path(wdir)
}
if (!file.exists(in.modDir)) {
  stop("Invalid directory for module-table file, \"", in.modDir, "\". \n")
}

# Module-table - file name
if (!exists("mod.filename")) {
  stop("Object \'mod.filename\' not found. \n")
} else {
  aux <- file.path(in.modDir, paste(mod.filename, ".csv", sep = ""))
  if (!is.character(mod.filename)) {
    stop("Module-table input: file name should be of type character, \"", 
         class(mod.filename), "\". \n")
  } else if (!file.exists(aux)) {
    stop("File \"", aux, "\" not found. \n")
  } 
}

### Output files
# Output directory
if (!exists("out.dir")) {
  warning("Object \'out.dir\' not found, set to default value: \"",
          wdir, "/output/\".")
  out.dir <- file.path(wdir, "output/")
} else if (!is.character(out.dir)) {
  stop("Invalid data type for \'out.dir\', \"", class(out.dir), "\".")
} else if (out.dir == "" || out.dir == " ") {
  out.dir <- file.path(wdir)
}
if (!file.exists(out.dir)) {
  dir.create(file.path(out.dir))
} 

if (!exists("out.filename")) {
  warning("Object \'out.filename\' not found, set to default value: \"",
          "out\".")
  out.filename <- "out"
} else if (!is.character(out.filename)) {
  warning("Invalid data type for \'out.filename\', \"", class(out.filename),
          "\".\nObject \'out.filename\' is set to \"out\".")
  out.filename <- "out"
}

## used in run 1 and run 2 only
if (!exists("out.ssplot")) {
  warning("Object \'out.ssplot\' not found, set to default value: \"",
          "FALSE\".")
  out.ssplot <- FALSE
} else if (!is.logical(out.ssplot)) {
  aux <- as.logical(out.ssplot)
  if (is.na(aux)) { # e.g. when out.ssplot set to  "" or " "
    warning("Invalid data type for \'out.ssplot\', \"", class(out.ssplot),
            "\".\nObject \'out.ssplot\' is set to FALSE.")
    out.ssplot <- FALSE
  } else {
    warning("Invalid data type for \'out.ssplot\', \"", class(out.ssplot), 
            "\".\nObject \'out.ssplot\' is set to ", aux, ".")
    out.ssplot <- aux
  }
}

if (!exists("out.varFiles")) {
  warning("Object \'out.varFiles\' not found, set to default value: \"",
          "FALSE\".")
  out.varFiles <- FALSE
} else if (!is.logical(out.varFiles)) {
  aux <- as.logical(out.varFiles)
  if (is.na(aux)) { # e.g. when out.varFiles set to  "" or " "
    warning("Invalid data type for \'out.varFiles\', \"", class(out.varFiles),
            "\".\nObject \'out.varFiles\' is set to FALSE.")
    out.varFiles <- FALSE
  } else {
    warning("Invalid data type for \'out.varFiles\', \"", class(out.varFiles), 
            "\".\nObject \'out.varFiles\' is set to ", aux, ".")
    out.varFiles <- aux
  }
}

### Quality trimming options
if (run[1] == 1) {
  if (!exists("qtrim.thold")) {
    warning("Object \'qtrim.thold\' not found, set to default value: 10.")
    qtrim.thold <- as.integer(10)
  } else {
    aux <- as.integer(qtrim.thold)
    if (!is.numeric(qtrim.thold)) {
      if (!is.na(aux)) {
        warning("Invalid data type for \'qtrim.thold\', \"", class(qtrim.thold), 
                "\".\nObject \'qtrim.thold\' is set to ", aux, ".")
      } else {
        warning("Invalid data type for \'qtrim.thold\', \"", class(qtrim.thold), 
                "\".\nObject \'qtrim.thold\' is set to \"10\".")
        qtrim.thold <- as.integer(10)
      }
    } else if (aux != qtrim.thold) {
      warning("Invalid value for \'qtrim.thold\', ", qtrim.thold, ".\n", 
              "Object \'qtrim.thold\' is set to ", aux, ".")
    }
    rm(aux)
    qtrim.thold <- as.integer(qtrim.thold)
    if (qtrim.thold < 0) {
      warning("Invalid value for \'qtrim.thold\', ", qtrim.thold, ".\n", 
              "Object \'qtrim.thold\' is set to 10.")
      qtrim.thold <- as.integer(10)
    } else if (qtrim.thold > 41) {
      warning("Check quality threshold: too high value, \"", qtrim.thold, 
              "\".")
    }
  }
}

# qtrim.3end and qtrim.flag needed for NameGen.R
if (!exists("qtrim.3end")) {
  warning("Object \'qtrim.3end\' not found, set to default value: 1.")
  qtrim.3end <- as.integer(1)
} else {
  aux <- as.integer(qtrim.3end)
  if (!is.numeric(qtrim.3end)) {
    if (!is.na(qtrim.3end)) {
      warning("Invalid data type for \'qtrim.3end\', \"", class(qtrim.3end), 
              "\".\nObject \'qtrim.end\' is set to ", aux,".")
    } else {
      warning("Invalid data type for \'qtrim.3end\', \"", class(qtrim.3end), 
              "\".\nObject \'qtrim.end\' is set to 1.")
      qtrim.3end <- as.integer(1)
    }
  } else if (aux != qtrim.3end) {
    warning("Invalid value for \'qtrim.3end\', ", qtrim.3end, ".\nObject ", 
            "\'qtrim.end\' is set to ", aux,".")
  }
  rm(aux)
  qtrim.3end <- as.integer(qtrim.3end)
  if (qtrim.3end != 0 && qtrim.3end != 1) {
    warning("Invalid value for \'qtrim.3end\', ", qtrim.3end, ".\nObject ", 
            "\'qtrim.end\' should be either 0 or 1.\nObject \'qtrim.end\' ",
            "is set to 1.")
    qtrim.3end <- as.integer(1)
  }
}

if (!exists("qtrim.flag")) {
  warning("Object \'qtrim.flag\' not found, set to default value: 1.")
  qtrim.flag <- as.integer(1)
} else {
  aux <- as.integer(qtrim.flag)
  if (!is.numeric(qtrim.flag)) { 
    if (!is.na(aux)) {
      warning("Invalid data type for \'qtrim.flag\', \"", class(qtrim.flag), 
              "\".\nObject \'qtrim.flag\' is set to ", aux, ".")
      qtrim.flag <- aux
    } else {
      warning("Invalid data type for \'qtrim.flag\', \"", class(qtrim.flag), 
              "\".\nObject \'qtrim.flag\' is set to 1.")
      qtrim.flag <- as.integer(1)
    }    
  } else if(aux != qtrim.flag) {
    warning("Invalid value for \'qtrim.flag\', ", qtrim.flag, ".\nObject", 
            " \'qtrim.flag\' is set to ", aux,".")
    qtrim.flag <- aux
  }
  rm(aux)
  qtrim.flag <- as.integer(qtrim.flag)
  if (qtrim.flag != 0 && qtrim.flag != 1) {
    warning("Invalid value for \'qtrim.flag\', ", qtrim.flag, ".\nObject", 
            " \'qtrim.flag\' should be either 0 or 1.\nObject \'",
            "qtrim.flag\' is set to 1.")
    qtrim.flag <- as.integer(1)
  } 
}

### Paired-end read assembly options
if (run[2] == 1) { 
  if (!exists("pandaseq.path")) {
    warning("Object \'pandaseq.path\' not found, set to default value: ", 
            "character(0).")
    pandaseq.path <- character(0)
  } else {
    if (!is.character(pandaseq.path)) {
      warning("Invalid data type for \'pandaseq.path\', \"", 
              class(pandaseq.path), "\".\nObject \'pandaseq.path\' is set to", 
              " character(0).")
      pandaseq.path <- character(0)
    } else if (length(pandaseq.path) > 0 && (pandaseq.path == "" || pandaseq.path == " ")) {
      pandaseq.path <- character(0)
    }    
  }
  if (length(pandaseq.path) > 0) {
    if (!file.exists(pandaseq.path)) {
      stop("Invalid path for pandaseq binaries, \"", pandaseq.path, "\". \n")
    }
    if (!file.exists(file.path(pandaseq.path, "pandaseq"))) {
      stop("pandaseq binary file not found, \"", 
           file.path(pandaseq.path, "pandaseq"), "\". \n")
    }
  }
} 

if (run[3] == 1 || run[4] == 1) {
  if (!exists("paired.flag")) {
    warning("Object \'paired.flag\' not found, set to default value: 1.")
    paired.flag <- as.integer(1)
  } else {
    aux <- as.integer(paired.flag)
    if (!is.numeric(paired.flag)) {
      if (!is.na(aux)) {
        warning("Invalid data type for \'paired.flag\', \"", class(paired.flag), 
                "\".\nObject \'paired.flag\' is set to ", aux, ".")
        paired.flag <- aux
      } else {
        warning("Invalid data type for \'paired.flag\', \"", class(paired.flag), 
                "\".\nObject \'paired.flag\' is set to 1.")
        paired.flag <- as.integer(1)
      }
    } else if(aux != paired.flag) {
      warning("Invalid value for \'paired.flag\', ", paired.flag, ".\nObject",
              " \'paired.flag\' is set to ", aux, ".")
      paired.flag <- aux
    }
    rm(aux)
    paired.flag <- as.integer(paired.flag)
    if (paired.flag != 0 && paired.flag != 1) {
      warning("Invalid value for \'paired.flag\', ", paired.flag, ".\nObject",
              " \'paired.flag\' should be either 0 or 1.\nObject \'",
              "paired.flag\' is set to 1.")
      paired.flag <- as.integer(1)
    }
  }
  
  if (seq.mode == "PE") {
    if (!exists("paired.file")) {
      warning("Object \'paired.file\' not found, set to default option: \"f\".")
      paired.file <- "f"
    } else {
      if (!is.character(paired.file)) {
        warning("Invalid data type for \'paired.file\', \"", class(paired.file), 
                "\".\nObject \'paired.file\' is set to \"f\".")
        paired.file <- "f"
      }
      if (paired.file != "f" && paired.file != "r") {
        warning("Invalid entry for \'paired.file\', \"", paired.file, "\".\n",
                "Object \'paired.file\' should be either \"f\" or \"r\".\n", 
                "Object \'paired.file\' is set to \"f\".")
        paired.file <- "f"
      }
    }
  }
}

### Read mapping options
if (run[3] == 1 || run[4] == 1 || run[5] == 1) {
  if (!exists("map.mode")) {
    warning("Object \'map.mode\' not found, set to default option: \"bwa\".")
    map.mode <- "bwa"
  } else {
    if (!is.character(map.mode)) {
      warning("Invalid data type for \'map.mode\', \"", class(map.mode), 
              "\".\nObject \'map.mode\' is set to \"bwa\".")
      map.mode <- "bwa"
    }
    if (map.mode != "gls" && map.mode != "grPA" && map.mode != "bwa") {
      warning("Invalid entry for \'map.mode\', \"", map.mode, "\".\nObject \'", 
              "map.mode\' should be either \"gls\", \"grPA\" or \"bwa\".\n", 
              "Object \'map.mode\' is set to \"bwa\".")
      map.mode <- "bwa"
    }
  }
  
  if (map.mode == "gls") {
    
    if (!exists("gls.ambiguity")) {
      warning("Object \'gls.ambiguity\' not found, set to default option: \"", 
              "TRUE\".")
      gls.ambiguity <- TRUE
    } else if (!is.logical(gls.ambiguity)) {
      aux <- as.logical(gls.ambiguity)
      if (is.logical(aux)) {
        warning("Invalid data type for \'gls.ambiguity\', \"", 
                class(gls.ambiguity), "\".\nObject \'gls.ambiguity\' is set ",
                "to ", aux, ".")
        gls.ambiguity <- aux
      } else {
        warning("Invalid data type for \'gls.ambiguity\', \"", 
                class(gls.ambiguity), "\".\nObject \'gls.ambiguity\' is set ",
                "to TRUE.")
        gls.ambiguity <- TRUE
      }
      rm(aux)
    }
    
    if (!exists("gls.direction")) {
      warning("Object \'gls.direction\' not found, set to default option: \"",
              "r\".")
      gls.direction <- "r"
    } else {
      if (!is.character(gls.direction)) {
        warning("Invalid data type for \'gls.direction\', \"", 
                class(gls.direction), "\".\nObject \'gls.direction\' is set",
                " to \"r\".")
        gls.direction <- "r"
      } else if (gls.direction != "f" && gls.direction != "r" ) {
        warning("Invalid entry for \'gls.direction\', \"", gls.direction, "\"",
              ".\nObject \'gls.direction\' should be either \"f\" or \"r\".\n", 
              "Object \'gls.direction\' is set to \"r\".")
        gls.direction <- "r"
      }
    }
    
    if (!exists("gls.mma")) {
      warning("Object \'gls.mma\' not found, set to default value: 0.")
      gls.mma <- as.integer(0)
    } else {
      aux <- as.integer(gls.mma)
      if (!is.numeric(gls.mma)) {
        if (!is.na(aux)) {
          warning("Invalid data type for \'gls.mma\', \"", class(gls.mma), "\"",
                  ".\nObject \'gls.mma\' is set to ", aux, ".")
          gls.mma <- aux
        } else {
          warning("Invalid data type for \'gls.mma\', \"", class(gls.mma), "\"",
                  ".\nObject \'gls.mma\' is set to 0.")
          gls.mma <- as.integer(0)
        }
      } else if (aux != gls.mma) {
        warning("Invalid value for \'gls.mma\', \"", gls.mma, "\".\nObject \'",
                "gls.mma\' is set to ", aux,".")
        gls.mma <- aux
      }
      rm(aux)
      gls.mma <- as.integer(gls.mma)
    }
    
  } else if (map.mode == "bwa") {
    
    if (run[3] == 1) {
      if (!exists("bwa.path")) {
        warning("Object \'bwa.path\' not found, set to default value: ", 
                "character(0).")
        bwa.path <- character(0)
      } else {
        if (!is.character(bwa.path)) {
          warning("Invalid data type for \'bwa.path\', \"", class(bwa.path), "\"",
                  ".\nObject \'bwa.path\' is set to character(0).")
          bwa.path <- character(0)
        } else if (length(bwa.path) > 0 && (bwa.path == "" || bwa.path == " ")) {
          bwa.path <- character(0)
        }    
      }
      if (length(bwa.path) > 0) {
        if (!file.exists(bwa.path)) {
          stop("Invalid path for bwa binaries, \"", bwa.path, "\". \n")
        }
        if (!file.exists(file.path(bwa.path, "bwa"))) {
          stop("bwa binaries not found, \"", file.path(bwa.path, "bwa"), 
               "\". \n")
        }
      }   
      
      if (!exists("bwa.cVal")) {
        warning("Object \'bwa.cVal\' not found, set to default value: 20000.")
        bwa.cVal <- as.integer(20000)
      } else {
        aux <- as.integer(bwa.cVal)
        if (!is.numeric(bwa.cVal)) {      
          if (!is.na(aux)) {
            warning("Invalid data type for \'bwa.cVal\', \"", class(bwa.cVal), 
                    "\".\nObject \'bwa.cVal\' is set to ", aux, ".")
          } else {
            warning("Invalid data type for \'bwa.cVal\', \"", class(bwa.cVal), 
                    "\".\nObject \'bwa.cVal\' is set to 20000.")
            bwa.cVal <- as.integer(20000)
          }
        } else if (aux != bwa.cVal) {
          warning("Invalid value for \'bwa.cVal\', ", bwa.cVal, ".\nObject \'", 
                  "bwa.cVal\' is set to ", aux, ".")
        }
        rm(aux)
        bwa.cVal <- as.integer(bwa.cVal)
        if (bwa.cVal < 10000) {
          warning("Invalid value for \'bwa.cVal\', ", bwa.cVal, ".\nObject \'", 
                  "bwa.cVal\' is set to 20000.")
          bwa.cVal <- as.integer(20000)
        } 
      }
      
      if (!exists("gatk.path")) {
        warning("Object \'gatk.path\' not found, set to default value: ", 
                "character(0).")
        gatk.path <- character(0)
      } else {
        if (!is.character(gatk.path)) {
          warning("Invalid data type for \'gatk.path\', \"", class(gatk.path), 
                  "\".\nObject \'gatk.path\' is set to character(0).")
          gatk.path <- character(0)
        } else if (length(gatk.path) > 0 && 
                     (gatk.path == "" || gatk.path == " ")) {
          gatk.path <- character(0)
        }    
      }
      if (length(gatk.path) > 0) {
        if (!file.exists(gatk.path)) {
          stop("Invalid path for gatk's jar file, \"", gatk.path, "\". \n")
        }
        if (!file.exists(file.path(gatk.path, "GenomeAnalysisTK.jar"))) {
          stop("gatk's jar file not found, \"", 
               file.path(gatk.path, "GenomeAnalysisTK.jar"), "\". \n")
        }
      }  
      
      if (!exists("picard.path")) {
        warning("Object \'picard.path\' not found, set to default value: ", 
                "character(0).")
        picard.path <- character(0)
      } else {
        if (!is.character(picard.path)) {
          warning("Invalid data type for \'picard.path\', \"", class(picard.path),
                  "\".\nObject \'picard.path\' is set to character(0).")
          picard.path <- character(0)
        } else if (length(picard.path) > 0 && 
                     (picard.path == "" || picard.path == " ")) {
          picard.path <- character(0)
        }    
      }
      if (length(picard.path) > 0) {
        if (!file.exists(picard.path)) {
          stop("Invalid path for picard's jar file, \"", picard.path, "\". \n")
        }
        if (!file.exists(file.path(picard.path, "picard.jar"))) {
          stop("picard jar file not found, \"", 
               file.path(picard.path, "picard.jar"), "\". \n")
        }
      } 
      
      if (!exists("samtools.path")) {
        warning("Object \'samtools.path\' not found, set to default value: ", 
                "character(0).")
        samtools.path <- character(0)
      } else {
        if (!is.character(samtools.path)) {
          warning("Invalid data type for \'samtools.path\', \"", class(samtools.path),
                  "\".\nObject \'samtools.path\' is set to character(0).")
          samtools.path <- character(0)
        } else if (length(samtools.path) > 0 && 
                     (samtools.path == "" || samtools.path == " ")) {
          samtools.path <- character(0)
        }    
      }
      if (length(samtools.path) > 0) {
        if (!file.exists(samtools.path)) {
          stop("Invalid path for samtools binary, \"", samtools.path, "\". \n")
        }
        if (!file.exists(file.path(samtools.path, "samtools"))) {
          stop("samtools binary not found, \"", 
               file.path(samtools.path, "samtools"), "\". \n")
        }
      }
      
    }
    
    if (run[4] == 1 || run[5] == 1) {
      if (!exists("bwa.dupl")) {
        warning("Object \'bwa.dupl\' not found, set to default option: FALSE.")
        bwa.dupl <- FALSE
      } else if (!is.logical(bwa.dupl)) {
        aux <- as.logical(bwa.dupl)
        if (is.logical(aux)) {
          warning("Invalid data type for \'bwa.dupl\', \"", class(bwa.dupl), 
                  "\".\nObject \'bwa.dupl\' is set to ", aux, ".")
          bwa.dupl <- aux
        } else {
          warning("Invalid data type for \'bwa.dupl\', \"", class(bwa.dupl), 
                  "\".\nObject \'bwa.dupl\' is set to FALSE.")
          bwa.dupl <- FALSE
        }
        rm(aux)
      }
      
      if (!exists("coverage.left")) {
        warning("Object \'coverage.left\' not found, set to default value: 0.")
        coverage.left <- as.integer(0)
      } else {
        aux <- as.integer(coverage.left)
        if (!is.numeric(coverage.left)) {
          if (!is.na(aux)) {
            warning("Invalid data type for \'coverage.left\', \"", 
                    class(coverage.left), "\".\nObject \'coverage.left\' is set",
                    " to ", aux, ".")
          } else {
            warning("Invalid data type for \'coverage.left\', \"", 
                    class(coverage.left), "\".\nObject \'coverage.left\' is set", 
                    " to 0.")
            coverage.left <- as.integer(0)
          }
        } else if (aux != coverage.left) {
          warning("Invalid value for \'coverage.left\', ", coverage.left, ".\n", 
                  "Object \'coverage.left\' is set to ", aux, ".")
        }
        rm(aux)
        coverage.left <- as.integer(coverage.left)
        if (coverage.left < 0) {
          warning("Invalid value for \'coverage.left\', ", coverage.left, ".\n", 
                  "Object \'coverage.left\' is set to 0.")
          coverage.left <- as.integer(0)
        } 
      }
      
      if (!exists("coverage.right")) {
        warning("Object \'coverage.right\' not found, set to default value: 0.")
        coverage.right <- as.integer(0)
      } else {
        aux <- as.integer(coverage.right)
        if (!is.numeric(coverage.right)) {
          if (!is.na(aux)) {
            warning("Invalid data type for \'coverage.right\', \"", 
                    class(coverage.right), "\".\nObject \'coverage.right\' is set",
                    " to ", aux, ".")
          } else {
            warning("Invalid data type for \'coverage.right\', \"", 
                    class(coverage.right), "\".\nObject \'coverage.right\' is set", 
                    " to 0.")
            coverage.right <- as.integer(0)
          }
        } else if (aux != coverage.right) {
          warning("Invalid value for \'coverage.right\', ", coverage.right, ".\n", 
                  "Object \'coverage.right\' is set to ", aux, ".")
        }
        rm(aux)
        coverage.right <- as.integer(coverage.right)
        if (coverage.right < 0) {
          warning("Invalid value for \'coverage.right\', ", coverage.right, ".\n", 
                  "Object \'coverage.right\' is set to 0.")
          coverage.right <- as.integer(0)
        } 
      }
      
      if (!exists("mapQ.thold")) {
        warning("Object \'mapQ.thold\' not found, set to default value: 8.")
        mapQ.thold <- as.integer(8)
      } else {
        aux <- as.integer(mapQ.thold)
        if (!is.numeric(mapQ.thold)) {
          if (!is.na(aux)) {
            warning("Invalid data type for \'mapQ.thold\', \"", class(mapQ.thold),
                    "\".\nObject \'mapQ.thold\' is set to ", aux, ".")
          } else {
            warning("Invalid data type for \'mapQ.thold\', \"", class(mapQ.thold), 
                    "\".\nObject \'mapQ.thold\' is set to 8.")
            mapQ.thold <- as.integer(8)
          }
        } else if (aux != mapQ.thold) {
          warning("Invalid value for \'mapQ.thold\', ", mapQ.thold, ".\nObject ", 
                  "\'mapQ.thold\' is set to ", aux, ".")
        }
        rm(aux)
        mapQ.thold <- as.integer(mapQ.thold)
        if (mapQ.thold < 0) {
          warning("Invalid value for \'mapQ.thold\', ", mapQ.thold, ".\nObject ", 
                  "\'mapQ.thold\' is set to 8.")
          mapQ.thold <- as.integer(8)
        } 
      }
      
      if (!exists("editDist.thold")) {
        warning("Object \'editDist.thold\' not found, set to default value: 8.")
        editDist.thold <- as.integer(8)
      } else {
        aux <- as.integer(editDist.thold)
        if (!is.numeric(editDist.thold)) {
          if (!is.na(aux)) {
            warning("Invalid data type for \'editDist.thold\', \"", 
                    class(editDist.thold), "\".\nObject \'editDist.thold\' is set",
                    " to ", aux, ".")
          } else {
            warning("Invalid data type for \'editDist.thold\', \"", 
                    class(editDist.thold), "\".\nObject \'editDist.thold\' is set",
                    " to 8.")
            editDist.thold <- as.integer(8)
          }
        } else if (aux != editDist.thold) {
          warning("Invalid value for \'editDist.thold\', ", editDist.thold, ".\n", 
                  "Object \'editDist.thold\' is set to ", aux, ".")
        }
        rm(aux)
        editDist.thold <- as.integer(editDist.thold)
        if (editDist.thold < 0) {
          warning("Invalid value for \'editDist.thold\', ", editDist.thold, ".\n", 
                  "Object \'editDist.thold\' is set to 8.")
          editDist.thold <- as.integer(8)
        } 
      }
    }
    
    if (run[5] == 1) {
      
      if (!exists("mismatch.filter")) {
        warning("Object \'mismatch.filter\' not found, set to default option:", 
                " TRUE.")
        mismatch.filter <- TRUE
      } else if (!is.logical(mismatch.filter)) {
        aux <- as.logical(mismatch.filter)
        if (is.logical(aux)) {
          warning("Invalid data type for \'mismatch.filter\', \"", 
                  class(mismatch.filter), "\".\nObject \'mismatch.filter\' is",
                  " set to ", aux, ".")
          mismatch.filter <- aux
        } else {
          warning("Invalid data type for \'mismatch.filter\', \"", 
                  class(mismatch.filter), "\".\nObject \'mismatch.filter\' is",
                  " set to TRUE.")
          mismatch.filter <- TRUE
        }
        rm(aux)
      }
      
      if (mismatch.filter) {
        if (!exists("mismatch.qthold")) {
          warning("Object \'mismatch.qthold\' not found, set to default value:", 
                  " 15.")
          mismatch.qthold <- as.integer(15)
        } else {
          aux <- as.integer(mismatch.qthold)
          if (!is.numeric(mismatch.qthold)) {
            if (!is.na(aux)) {
              warning("Invalid data type for \'mismatch.qthold\', \"", 
                      class(mismatch.qthold), "\".\nObject \'mismatch.qthold\'", 
                      " is set to ", aux, ".")
            } else {
              warning("Invalid data type for \'mismatch.qthold\', \"", 
                      class(mismatch.qthold), "\".\nObject \'mismatch.qthold\'", 
                      " is set to 15.")
              mismatch.qthold <- as.integer(15)
            }
          } else if (aux != mismatch.qthold) {
            warning("Invalid value for \'mismatch.qthold\', ", mismatch.qthold, 
                    ".\nObject \'mismatch.qthold\' is set to ", aux, ".")
          }
          rm(aux)
          mismatch.qthold <- as.integer(mismatch.qthold)
          if (mismatch.qthold < 0) {
            warning("Invalid value for \'mismatch.qthold\', ", mismatch.qthold, 
                    ".\nObject \'mismatch.qthold\' is set to 15.")
            mismatch.qthold <- as.integer(15)
          }
        }
      }
      
      if (!exists("min.coverage")) {
        warning("Object \'min.coverage\' not found, set to default value: 20.")
        min.coverage <- as.integer(20)
      } else {
        aux <- as.integer(min.coverage)
        if (!is.numeric(min.coverage)) {
          if (!is.na(aux)) {
            warning("Invalid data type for \'min.coverage\', \"", 
                    class(min.coverage), "\".\nObject \'min.coverage\' is set", 
                    " to ", aux, ".")
          } else {
            warning("Invalid data type for \'min.coverage\', \"", 
                    class(min.coverage), "\".\nObject \'min.coverage\' is set", 
                    "to 20.")
            min.coverage <- as.integer(20)
          }
        } else if (aux != min.coverage) {
          warning("Invalid value for \'min.coverage\', ", min.coverage,
                  ".\nObject \'min.coverage\' is set to ", aux, ".")
        }
        rm(aux)
        min.coverage <- as.integer(min.coverage)
        if (min.coverage < 0) {
          warning("Invalid value for \'min.coverage\', ", min.coverage, 
                  ".\nObject \'min.coverage\' is set to 20.")
          min.coverage <- as.integer(20)
        }
      }
      
      if (!exists("seq.error")) {
        stop("Object \'seq.error\' not found.")
      } else {
        aux <- as.numeric(seq.error)
        if (!is.numeric(seq.error)) {
          if (!is.na(aux)) {
            warning("Invalid data type for \'seq.error\', \"", class(seq.error), 
                    "\".\nObject \'seq.error\' is set to ", aux, ".")
            seq.error <- aux
          } else {
            stop("Invalid data type for \'seq.error\', \"", class(seq.error), 
                 "\".")
          }
        }
        rm(aux)
      }
  
    }
  }
}

if (!exists("mem.trace")) {
  warning("Object \'mem.trace\' not found, set to default option:", 
          " FALSE.")
  mem.trace <- FALSE
} else if (!is.logical(mem.trace)) {
  aux <- as.logical(mem.trace)
  if (is.logical(aux)) {
    warning("Invalid data type for \'mem.trace\', \"", class(mem.trace), 
            "\".\nObject \'mem.trace\' is set to ", aux, ".")
    mem.trace <- aux
  } else {
    warning("Invalid data type for \'mem.trace\', \"", 
            class(mem.trace), "\".\nObject \'mem.trace\' is",
            " set to FALSE.")
    mem.trace <- FALSE
  }
  rm(aux)
}
